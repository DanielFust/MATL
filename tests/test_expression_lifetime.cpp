/**************************************************************************************************
 * File: test_expression_lifetime.cpp
 *
 * Description: Expressions are safe to store. An expression object (e.g. 'auto e = A + B + C;')
 *              must remain valid after the statement that created it:
 *
 *                - sub-expressions are stored by value
 *                - named (lvalue) matrices are stored by reference, so a stored expression is a
 *                  live formula: it reads the matrices' current values when it is evaluated
 *                - temporary (rvalue) matrices are moved into the expression, which owns them
 *
 *              Several of these cases read freed stack memory if sub-expressions are held by
 *              reference; run 'make check-sanitize' to detect that reliably.
 **************************************************************************************************/
#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    RowMajorMatrix<double> make_matrix(index_t m, index_t n, double start)
    {
        RowMajorMatrix<double> result(m, n);
        matl_test::fill_sequence(result, start);
        return result;
    }

    FixedSizeRowMajorMatrix<double, 2, 2> make_fixed(double start)
    {
        FixedSizeRowMajorMatrix<double, 2, 2> result;
        matl_test::fill_sequence(result, start);
        return result;
    }

    // Builds an expression from function parameters and returns it. The named matrices outlive
    // the call, so the returned expression is valid.
    auto weighted_sum(RowMajorMatrix<double> const& X, RowMajorMatrix<double> const& Y)
    {return 2.0 * X + Y - 1.0;}
}

/************************* Sub-expressions are stored by value ****************************/

TEST_CASE("Lifetime: stored chain of elementwise operations")
{
    RowMajorMatrix<double> A = make_matrix(3, 3, 1.0), B = make_matrix(3, 3, 5.0), C = make_matrix(3, 3, -2.0);
    Ref a = Ref::from(A), b = Ref::from(B), c = Ref::from(C);

    auto e = A + B + C;                       // (A + B) is a temporary that dies here
    auto f = (A - B) * 2.0 - transpose(C) / 4.0 + 1.0;

    RowMajorMatrix<double> R(3, 3);
    R = e;  CHECK_MATRIX_EQ(R, a + b + c);
    R = f;  CHECK_MATRIX_EQ(R, (a - b) * 2.0 - transpose(c) / 4.0 + 1.0);
}

TEST_CASE("Lifetime: stored expressions can be combined into new expressions")
{
    RowMajorMatrix<double> A = make_matrix(2, 3, 1.0), B = make_matrix(2, 3, 4.0);
    Ref a = Ref::from(A), b = Ref::from(B);

    auto e1 = A + B;
    auto e2 = e1 * 3.0;
    auto e3 = e2 - e1 + transpose(transpose(A));

    RowMajorMatrix<double> R(2, 3);
    R = e3;  CHECK_MATRIX_EQ(R, (a + b) * 3.0 - (a + b) + a);
}

TEST_CASE("Lifetime: expressions are independent copies of their sub-expressions")
{
    RowMajorMatrix<double> A = make_matrix(2, 2, 1.0);
    Ref a = Ref::from(A);
    RowMajorMatrix<double> R(2, 2);

    auto e = A * 2.0 + 1.0;
    {
        auto copy = e;                        // copy, then let the copy go out of scope
        R = copy;
        CHECK_MATRIX_EQ(R, a * 2.0 + 1.0);
    }
    R = e;  CHECK_MATRIX_EQ(R, a * 2.0 + 1.0);
}

TEST_CASE("Lifetime: an expression can be returned from a function")
{
    RowMajorMatrix<double> X = make_matrix(3, 2, 1.0), Y = make_matrix(3, 2, 10.0);
    auto e = weighted_sum(X, Y);
    RowMajorMatrix<double> R(3, 2);
    R = e;  CHECK_MATRIX_EQ(R, 2.0 * Ref::from(X) + Ref::from(Y) - 1.0);
}

TEST_CASE("Lifetime: sub-expressions are held by value, named matrices by reference")
{
    using M   = RowMajorMatrix<double>;
    using Sum = expr::Elementwise<expr::Add, M, M>;
    static_assert(std::is_reference_v<detail::stored_t<M>>, "named matrices are stored by reference");
    static_assert(!std::is_reference_v<detail::stored_t<Sum>>, "sub-expressions are stored by value");
    static_assert(!std::is_reference_v<detail::stored_t<expr::Elementwise<expr::Multiply, M, double>>>);
    static_assert(!std::is_reference_v<detail::stored_t<expr::Transpose<M>>>);
    static_assert(!std::is_reference_v<detail::stored_t<double>>, "scalars are stored by value");

    // The operators build these expression types
    M A(2, 2), B(2, 2);
    static_assert(std::is_same_v<decltype(A + B), Sum>);
    static_assert(std::is_same_v<decltype(A * 2.0), expr::Elementwise<expr::Multiply, M, double>>);
    static_assert(std::is_same_v<decltype(2.0 - A), expr::Elementwise<expr::Subtract, double, M>>);
    static_assert(std::is_same_v<decltype(A / 2), expr::Elementwise<expr::Divide, M, int>>);
    static_assert(std::is_same_v<decltype(A / 2.0), expr::Elementwise<expr::Multiply, M, double>>);
    static_assert(std::is_same_v<decltype((A + B) * A), expr::Product<Sum, M>>);
    static_assert(std::is_same_v<decltype(transpose(A)), expr::Transpose<M>>);
    CHECK(true);
}

/********************** Named matrices are referenced (live formula) ************************/

TEST_CASE("Lifetime: a stored expression reads the current values of named matrices")
{
    RowMajorMatrix<double> A(2, 2), B(2, 2), R(2, 2);
    A = 1.0;
    B = 2.0;
    auto e = A + B * 10.0;

    R = e;   CHECK_MATRIX_EQ(R, Ref(2, 2, 21.0));
    A = 5.0;
    R = e;   CHECK_MATRIX_EQ(R, Ref(2, 2, 25.0));
    B = 0.0;
    R = e;   CHECK_MATRIX_EQ(R, Ref(2, 2, 5.0));
}

TEST_CASE("Lifetime: a typed declaration evaluates immediately (value semantics)")
{
    RowMajorMatrix<double> A(2, 3), B(2, 3);
    A = 1.0;
    B = 2.0;
    Matrix<double> m = A + B;                 // value: evaluated now
    Matrix<double> p(transpose(A) * B);       // direct initialization from a product
    A = 100.0;
    CHECK_MATRIX_EQ(m, Ref(2, 3, 3.0));
    CHECK_MATRIX_EQ(p, Ref(3, 3, 4.0));

    Matrix<double, 2, 2> F, G;
    F = 1.0;
    G = 2.0;
    Matrix<double, 2, 2> f = F * G + 1.0;     // fixed-size value
    F = 0.0;
    CHECK_MATRIX_EQ(f, Ref(2, 2, 5.0));
}

TEST_CASE("Lifetime: a stored expression can initialize a matrix")
{
    RowMajorMatrix<double> A(2, 2);
    A = 3.0;
    auto e = A * 2.0 - 1.0;
    Matrix<double> m = e;
    CHECK_MATRIX_EQ(m, Ref(2, 2, 5.0));

    Matrix<double> s = e + A;
    CHECK_MATRIX_EQ(s, Ref(2, 2, 8.0));
}

/************************* Temporary matrices are moved in and owned ************************/

TEST_CASE("Lifetime: temporary matrices are owned by the expression")
{
    RowMajorMatrix<double> A = make_matrix(2, 2, 1.0);
    Ref a = Ref::from(A);
    Ref t = Ref::from(make_matrix(2, 2, 10.0));
    RowMajorMatrix<double> R(2, 2);

    auto e1 = A + make_matrix(2, 2, 10.0);    R = e1;  CHECK_MATRIX_EQ(R, a + t);
    auto e2 = make_matrix(2, 2, 10.0) + A;    R = e2;  CHECK_MATRIX_EQ(R, t + a);
    auto e3 = make_matrix(2, 2, 10.0) - A;    R = e3;  CHECK_MATRIX_EQ(R, t - a);
    auto e4 = A - make_matrix(2, 2, 10.0);    R = e4;  CHECK_MATRIX_EQ(R, a - t);
    auto e5 = make_matrix(2, 2, 10.0) + make_matrix(2, 2, 10.0);
    R = e5;  CHECK_MATRIX_EQ(R, t + t);
}

TEST_CASE("Lifetime: temporary matrices with scalars and transposes")
{
    Ref t = Ref::from(make_matrix(2, 3, 1.0));
    RowMajorMatrix<double> R(2, 3), T(3, 2);

    auto e1 = make_matrix(2, 3, 1.0) * 2.0;   R = e1;  CHECK_MATRIX_EQ(R, t * 2.0);
    auto e2 = 2.0 * make_matrix(2, 3, 1.0);   R = e2;  CHECK_MATRIX_EQ(R, 2.0 * t);
    auto e3 = make_matrix(2, 3, 1.0) + 1.0;   R = e3;  CHECK_MATRIX_EQ(R, t + 1.0);
    auto e4 = 1.0 - make_matrix(2, 3, 1.0);   R = e4;  CHECK_MATRIX_EQ(R, 1.0 - t);
    auto e5 = make_matrix(2, 3, 1.0) / 2.0;   R = e5;  CHECK_MATRIX_EQ(R, t / 2.0);
    auto e6 = make_matrix(2, 3, 1.0) / 2;     R = e6;  CHECK_MATRIX_EQ(R, t / 2.0);
    auto e7 = transpose(make_matrix(2, 3, 1.0)); T = e7; CHECK_MATRIX_EQ(T, transpose(t));
}

TEST_CASE("Lifetime: temporary matrices in products")
{
    RowMajorMatrix<double> A = make_matrix(2, 3, 1.0);
    Ref a = Ref::from(A), b = Ref::from(make_matrix(3, 2, -1.0));
    RowMajorMatrix<double> R(2, 2);

    auto e1 = A * make_matrix(3, 2, -1.0);                      R = e1;  CHECK_MATRIX_EQ(R, a * b);
    auto e2 = make_matrix(2, 3, 1.0) * make_matrix(3, 2, -1.0); R = e2;  CHECK_MATRIX_EQ(R, a * b);
    auto e3 = A * make_matrix(3, 2, -1.0) + 1.0;                R = e3;  CHECK_MATRIX_EQ(R, a * b + 1.0);
}

TEST_CASE("Lifetime: temporary identity matrix (eye) in a stored expression")
{
    RowMajorMatrix<double> A = make_matrix(3, 3, 1.0);
    Ref expected = Ref::from(A);
    for (index_t k = 1; k <= 3; ++k) expected(k, k) += 1.0;
    auto e = A + eye<double>(3);
    RowMajorMatrix<double> R(3, 3);
    R = e;  CHECK_MATRIX_EQ(R, expected);
}

TEST_CASE("Lifetime: temporary fixed-size matrices are owned by the expression")
{
    FixedSizeRowMajorMatrix<double, 2, 2> A = make_fixed(1.0), R;
    Ref a = Ref::from(A), t = Ref::from(make_fixed(5.0));
    auto e = A + make_fixed(5.0) * 2.0;
    R = e;  CHECK_MATRIX_EQ(R, a + t * 2.0);
}

TEST_CASE("Lifetime: a temporary matrix over external memory keeps referencing the buffer")
{
    RowMajorMatrix<double> A = make_matrix(2, 2, 1.0);
    double buffer[4] = {10, 20, 30, 40};
    auto e = A - RowMajorMatrix<double>(2, 2, buffer);   // the expression owns the view, not the data
    RowMajorMatrix<double> R(2, 2);
    R = e;  CHECK_MATRIX_EQ(R, Ref::from(A) - Ref::from(RowMajorMatrix<double>(2, 2, buffer)));
    buffer[0] = 0.0;                                     // the expression reads the current buffer
    R = e;  CHECK_EQ(R(1,1), 1.0);
}

/******************************** Mixed element types ***************************************/

TEST_CASE("Lifetime: nested product of mixed element types keeps the promoted type")
{
    RowMajorMatrix<int>    I(2, 2);
    RowMajorMatrix<double> D(2, 2), Z(2, 2), C(2, 2);
    I = 1;
    D = 0.25;
    Z = 0.0;
    C = I * D;        CHECK_MATRIX_EQ(C, Ref(2, 2, 0.5));
    C = I * D + Z;    CHECK_MATRIX_EQ(C, Ref(2, 2, 0.5));   // product temporary must be double

    // Elementwise operands of a product are evaluated into a temporary of the promoted type
    RowMajorMatrix<double> E(2, 2);
    E = 1.0;
    C = (I + D) * E;  CHECK_MATRIX_EQ(C, Ref(2, 2, 2.5));   // (1 + 0.25) * 2
    C = (I - D) * E;  CHECK_MATRIX_EQ(C, Ref(2, 2, 1.5));   // (1 - 0.25) * 2
}
