/**************************************************************************************************
 * File: test_functions.cpp
 *
 * Description: Elementwise functions: negation, the <cmath>-style functions (abs, sqrt, exp, ...),
 *              and times, rdivide and pow, compared with the same functions applied element by
 *              element to the reference matrices.
 **************************************************************************************************/
#include <cmath>
#include <type_traits>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using M = RowMajorMatrix<double>;

    // Reference: f applied to every element
    template<typename F>
    Ref map(Ref r, F f) {for (double& x : r.data) x = f(x); return r;}

    // Reference: f applied to corresponding elements
    template<typename F>
    Ref zip(Ref a, Ref const& b, F f) {for (std::size_t k = 0; k < a.data.size(); ++k) a.data[k] = f(a.data[k], b.data[k]); return a;}

    // A 3x4 matrix with values in (lo, hi)
    M values(double lo, double hi)
    {
        M A(3, 4);
        for (index_t i = 1; i <= A.size(); ++i) A(i) = lo + (hi - lo) * double(i) / double(A.size() + 1);
        return A;
    }

    // times(a, b) is defined
    template<typename A, typename B>
    concept Timesable = requires(A const& a, B const& b) { times(a, b); };
}

/*********************************************** Negation ****************************************/

TEST_CASE("Functions: negation")
{
    M A = values(-2, 3), B = values(1, 5), R(3, 4);
    Ref a = Ref::from(A), b = Ref::from(B);
    static_assert(std::is_same_v<decltype(-A), expr::Unary<expr::Negate, M>>);

    R = -A;             CHECK_MATRIX_EQ(R, -1.0 * a);
    R = -(A + B);       CHECK_MATRIX_EQ(R, -1.0 * (a + b));
    R = -(-A) + B;      CHECK_MATRIX_EQ(R, a + b);
    R = A - -B;         CHECK_MATRIX_EQ(R, a + b);

    Matrix<int> I = {{1, -2}}, J(1, 2);
    J = -I;
    CHECK_EQ(J(1,1), -1);
    CHECK_EQ(J(1,2), 2);

    Matrix<double, 2, 2> F = {{1, 2}, {3, 4}}, G;
    G = -F;
    CHECK_EQ(G(2,1), -3.0);
}

/******************************** Elementwise math functions ********************************/

TEST_CASE("Functions: elementwise math functions match <cmath>")
{
    M P = values(0.1, 9.0), S = values(-0.9, 0.9), A = values(-5.5, 5.5), R(3, 4);
    Ref p = Ref::from(P), s = Ref::from(S), a = Ref::from(A);

    R = abs(A);    CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::abs(x);}));
    R = sqrt(P);   CHECK_MATRIX_EQ(R, map(p, [](double x) {return std::sqrt(x);}));
    R = cbrt(A);   CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::cbrt(x);}));
    R = exp(S);    CHECK_MATRIX_EQ(R, map(s, [](double x) {return std::exp(x);}));
    R = log(P);    CHECK_MATRIX_EQ(R, map(p, [](double x) {return std::log(x);}));
    R = log2(P);   CHECK_MATRIX_EQ(R, map(p, [](double x) {return std::log2(x);}));
    R = log10(P);  CHECK_MATRIX_EQ(R, map(p, [](double x) {return std::log10(x);}));
    R = sin(A);    CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::sin(x);}));
    R = cos(A);    CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::cos(x);}));
    R = tan(S);    CHECK_MATRIX_EQ(R, map(s, [](double x) {return std::tan(x);}));
    R = asin(S);   CHECK_MATRIX_EQ(R, map(s, [](double x) {return std::asin(x);}));
    R = acos(S);   CHECK_MATRIX_EQ(R, map(s, [](double x) {return std::acos(x);}));
    R = atan(A);   CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::atan(x);}));
    R = sinh(S);   CHECK_MATRIX_EQ(R, map(s, [](double x) {return std::sinh(x);}));
    R = cosh(S);   CHECK_MATRIX_EQ(R, map(s, [](double x) {return std::cosh(x);}));
    R = tanh(A);   CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::tanh(x);}));
    R = floor(A);  CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::floor(x);}));
    R = ceil(A);   CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::ceil(x);}));
    R = round(A);  CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::round(x);}));

    // Composed with other operations, and as product operands
    R = sqrt(P) * 2.0 + exp(-abs(A));
    CHECK_MATRIX_EQ(R, map(p, [](double x) {return std::sqrt(x);}) * 2.0 + map(a, [](double x) {return std::exp(-std::abs(x));}));
    M T(3, 3);
    T = abs(A) * transpose(sin(A));
    CHECK_MATRIX_EQ(T, map(a, [](double x) {return std::abs(x);}) * transpose(map(a, [](double x) {return std::sin(x);})));
}

TEST_CASE("Functions: element types follow the function's result")
{
    Matrix<int> I = {{4, 9}, {16, -25}};
    static_assert(std::is_same_v<detail::dtype_t<decltype(sqrt(I))>, double>);
    static_assert(std::is_same_v<detail::dtype_t<decltype(abs(I))>, int>);
    static_assert(std::is_same_v<detail::dtype_t<decltype(pow(I, 2))>, double>);
    static_assert(std::is_same_v<detail::dtype_t<decltype(times(I, I))>, int>);

    Matrix<double> R(2, 2);
    R = sqrt(abs(I));
    CHECK_MATRIX_EQ(R, Ref::from(Matrix<double>{{2, 3}, {4, 5}}));

    Matrix<float, 1, 2> F = {{4.0f, 9.0f}};
    static_assert(std::is_same_v<detail::dtype_t<decltype(sqrt(F))>, float>);
}

/******************************* times, rdivide and pow *************************************/

TEST_CASE("Functions: times, rdivide and pow")
{
    M A = values(-3, 4), B = values(0.5, 2.5), R(3, 4);
    Ref a = Ref::from(A), b = Ref::from(B);

    R = times(A, B);       CHECK_MATRIX_EQ(R, zip(a, b, [](double x, double y) {return x * y;}));
    R = rdivide(A, B);     CHECK_MATRIX_EQ(R, zip(a, b, [](double x, double y) {return x / y;}));
    R = pow(B, A);         CHECK_MATRIX_EQ(R, zip(b, a, [](double x, double y) {return std::pow(x, y);}));

    // A scalar on either side
    R = pow(A, 2);         CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::pow(x, 2);}));
    R = pow(2.0, A);       CHECK_MATRIX_EQ(R, map(a, [](double x) {return std::pow(2.0, x);}));
    R = rdivide(1.0, B);   CHECK_MATRIX_EQ(R, map(b, [](double x) {return 1.0 / x;}));
    R = times(3, A);       CHECK_MATRIX_EQ(R, 3.0 * a);

    // Expression operands
    R = times(A + B, A - B);
    CHECK_MATRIX_EQ(R, zip(a + b, a - b, [](double x, double y) {return x * y;}));
}

TEST_CASE("Functions: binary functions check dimensions")
{
    using F23 = FixedSizeRowMajorMatrix<double, 2, 3>;
    using F32 = FixedSizeRowMajorMatrix<double, 3, 2>;
    static_assert(Timesable<F23, F23> && Timesable<F23, double> && Timesable<double, F23>);
    static_assert(!Timesable<F23, F32>,       "compile-time dimensions must match");
    static_assert(!Timesable<double, double>, "at least one operand is a matrix or expression");
    static_assert(Timesable<F23, M>,          "dynamic dimensions are checked at run time");

    #ifdef MATL_DIMENSION_CHECK
    M A(2, 2), B(3, 3), R(2, 2);
    A = 1.0; B = 1.0;
    CHECK_MATL_ERROR(dimension_error, R = times(A, B));
    CHECK_MATL_ERROR(dimension_error, R = pow(A, B));
    #endif
    CHECK(true);
}

/******************************** Lifetimes and aliasing ************************************/

TEST_CASE("Functions: stored expressions are live formulas that own temporaries")
{
    M A = values(1, 2), R(3, 4);
    auto e = sqrt(A) + abs(ones<double>(3, 4));     // ones(...) is moved into the expression
    A = 4.0;                                        // e reads A's current values
    R = e;
    CHECK_MATRIX_EQ(R, Ref(3, 4, 3.0));

    auto f = pow(zeros<double>(3, 4) + 2.0, 3);     // a stored expression of temporaries
    R = f;
    CHECK_MATRIX_EQ(R, Ref(3, 4, 8.0));
}

TEST_CASE("Functions: the target may appear in the expression")
{
    M S(3, 3), T(3, 3);
    for (index_t i = 1; i <= 9; ++i) S(i) = double(i) - 5.0;
    Ref s = Ref::from(S);

    T = S; T = -T;                  CHECK_MATRIX_EQ(T, -1.0 * s);
    T = S; T = abs(T) + T;          CHECK_MATRIX_EQ(T, map(s, [](double x) {return std::abs(x) + x;}));
    T = S; T = -transpose(T);       CHECK_MATRIX_EQ(T, -1.0 * transpose(s));
    T = S; T = times(T, transpose(T));
    CHECK_MATRIX_EQ(T, zip(s, transpose(s), [](double x, double y) {return x * y;}));
    T = S; T = abs(T) * T;          CHECK_MATRIX_EQ(T, map(s, [](double x) {return std::abs(x);}) * s);
}
