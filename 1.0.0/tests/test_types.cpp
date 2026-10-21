/**************************************************************************************************
 * File: test_types.cpp
 *
 * Description: Concepts, expression types and traits: which types the operators accept, the
 *              expression types they build, compile-time dimensions and storage order, and the
 *              assignment forms the constrained signatures allow.
 **************************************************************************************************/
#include <string>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using M  = RowMajorMatrix<double>;
    using F  = FixedSizeRowMajorMatrix<double, 2, 3>;
    using FT = FixedSizeRowMajorMatrix<double, 3, 2>;

    RowMajorMatrix<double> make_matrix(index_t m, index_t n, double value)
    {
        RowMajorMatrix<double> result(m, n);
        result = value;
        return result;
    }

    // Detects whether 'a + b' is a valid expression
    template<typename A, typename B>
    concept Addable = requires(A const& a, B const& b) { a + b; };
}

/*********************************************** Concepts ****************************************/

TEST_CASE("Types: concepts classify matrices, expressions and scalars")
{
    static_assert(MatrixType<M> && MatrixType<F> && MatrixType<M const&>);
    static_assert(!MatrixType<expr::Transpose<M>> && !MatrixType<double>);

    static_assert(ExpressionType<expr::Elementwise<expr::Add, M, M>>);
    static_assert(ExpressionType<expr::Product<M, F>> && ExpressionType<expr::Owned<M>>);
    static_assert(!ExpressionType<M>);

    static_assert(MatrixAlgebraObject<M> && MatrixAlgebraObject<expr::Transpose<M>>);
    static_assert(!MatrixAlgebraObject<double> && !MatrixAlgebraObject<std::string>);

    static_assert(Scalar<double> && Scalar<int> && Scalar<float>);
    static_assert(!Scalar<double*> && !Scalar<M>);

    static_assert(FixedSize<F> && FixedSize<expr::Transpose<F>> && !FixedSize<M>);
    CHECK(true);
}

TEST_CASE("Types: operators accept only MATL objects and arithmetic scalars")
{
    static_assert(Addable<M, M> && Addable<M, F> && Addable<M, double> && Addable<int, M>);
    static_assert(Addable<expr::Transpose<M>, M>);
    static_assert(!Addable<M, std::string>, "unrelated types are not operands");
    static_assert(!Addable<M, double*>,     "pointers are not scalars");
    static_assert(!Addable<double, double*>);
    CHECK(true);
}

/******************************************* Expression types ************************************/

TEST_CASE("Types: one elementwise expression type covers every elementwise operation")
{
    M A(2, 2), B(2, 2);
    A = 6.0; B = 2.0;
    RowMajorMatrix<double> R(2, 2);

    static_assert(std::is_same_v<decltype(A - B),   expr::Elementwise<expr::Subtract, M, M>>);
    static_assert(std::is_same_v<decltype(A + 1.0), expr::Elementwise<expr::Add, M, double>>);
    static_assert(std::is_same_v<decltype(1.0 + A), expr::Elementwise<expr::Add, double, M>>);
    static_assert(std::is_same_v<decltype(3 * A),   expr::Elementwise<expr::Multiply, int, M>>);

    R = A - B;    CHECK_MATRIX_EQ(R, Ref(2, 2, 4.0));
    R = A + 1.0;  CHECK_MATRIX_EQ(R, Ref(2, 2, 7.0));
    R = 1.0 + A;  CHECK_MATRIX_EQ(R, Ref(2, 2, 7.0));
    R = 1.0 - A;  CHECK_MATRIX_EQ(R, Ref(2, 2, -5.0));
    R = 3 * A;    CHECK_MATRIX_EQ(R, Ref(2, 2, 18.0));
    R = A / 4.0;  CHECK_MATRIX_EQ(R, Ref(2, 2, 1.5));
    R = A / 4;    CHECK_MATRIX_EQ(R, Ref(2, 2, 1.5));    // double / int divides exactly
}

TEST_CASE("Types: integer matrices divide by integers with integer division")
{
    RowMajorMatrix<int> I(1, 3), R(1, 3);
    I(1) = 7; I(2) = -7; I(3) = 9;
    R = I / 2;
    CHECK_EQ(R(1), 3);
    CHECK_EQ(R(2), -3);
    CHECK_EQ(R(3), 4);
}

TEST_CASE("Types: element types are promoted")
{
    using I = RowMajorMatrix<int>;
    static_assert(std::is_same_v<detail::dtype_t<decltype(I(1,1) + 1.5)>, double>);
    static_assert(std::is_same_v<detail::dtype_t<decltype(I(1,1) * I(1,1))>, int>);
    static_assert(std::is_same_v<detail::dtype_t<decltype(I(1,1) * M(1,1))>, double>);
    static_assert(std::is_same_v<detail::dtype_t<decltype(I(1,1) / 2.0f)>, float>);
    CHECK(true);
}

/********************************************* Traits ********************************************/

TEST_CASE("Types: fixed sizes and storage order propagate through expressions")
{
    using detail::traits;
    using detail::StorageOrder;

    // Fixed-size operands give fixed-size expressions with compile-time dimensions
    static_assert(traits<expr::Elementwise<expr::Add, F, F>>::fixed_size);
    static_assert(traits<expr::Elementwise<expr::Multiply, double, F>>::rows == 2);
    static_assert(traits<expr::Transpose<F>>::rows == 3 && traits<expr::Transpose<F>>::cols == 2);
    static_assert(traits<expr::Product<F, FT>>::rows == 2 && traits<expr::Product<F, FT>>::cols == 2);

    // Any dynamic operand makes the expression dynamic
    static_assert(!traits<expr::Elementwise<expr::Add, F, M>>::fixed_size);
    static_assert(!traits<expr::Product<F, M>>::fixed_size);

    // Temporaries hold the result in the matching matrix type
    static_assert(std::is_same_v<detail::result_matrix_t<expr::Product<F, FT>>, FixedSizeRowMajorMatrix<double, 2, 2>>);
    static_assert(std::is_same_v<detail::result_matrix_t<expr::Product<F, M>>, M>);

    // Transposing flips the storage order; mixing storage orders gives no common order
    static_assert(traits<M>::storage_order == StorageOrder::RowMajor);
    static_assert(traits<expr::Transpose<M>>::storage_order == StorageOrder::ColMajor);
    static_assert(traits<expr::Elementwise<expr::Add, M, double>>::storage_order == StorageOrder::RowMajor);
    static_assert(traits<expr::Elementwise<expr::Add, M, expr::Transpose<M>>>::storage_order == StorageOrder::Mixed);
    CHECK(true);
}

/********************************** Assignment forms ****************************************/

TEST_CASE("Types: compound assignment accepts temporary matrices")
{
    RowMajorMatrix<double> A(2, 2);
    A = 1.0;
    A += make_matrix(2, 2, 3.0);   CHECK_MATRIX_EQ(A, Ref(2, 2, 4.0));
    A -= make_matrix(2, 2, 1.0);   CHECK_MATRIX_EQ(A, Ref(2, 2, 3.0));
    A *= make_matrix(2, 2, 0.5);   CHECK_MATRIX_EQ(A, Ref(2, 2, 3.0));   // 2x2 product
}

TEST_CASE("Types: scalar assignment returns the matrix")
{
    RowMajorMatrix<double> A(2, 2);
    (A = 1.0) += 2.0;
    CHECK_MATRIX_EQ(A, Ref(2, 2, 3.0));

    FixedSizeRowMajorMatrix<double, 2, 2> Fx;
    (Fx = 2.0) = 5.0;
    CHECK_MATRIX_EQ(Fx, Ref(2, 2, 5.0));
}

TEST_CASE("Types: matrices can be assigned from const matrices")
{
    RowMajorMatrix<double> A(2, 3);
    matl_test::fill_sequence(A);
    F Fx;
    matl_test::fill_sequence(Fx, 10.0);

    RowMajorMatrix<double> const& a  = A;
    F const&                      fx = Fx;

    F G;
    G = fx;  CHECK_MATRIX_EQ(G, Ref::from(Fx));    // fixed-size from const fixed-size
    G = a;   CHECK_MATRIX_EQ(G, Ref::from(A));     // fixed-size from const dynamic
    RowMajorMatrix<double> R(2, 3);
    R = fx;  CHECK_MATRIX_EQ(R, Ref::from(Fx));    // dynamic from const fixed-size

    static_assert(std::is_trivially_copyable_v<F>, "fixed-size matrices are plain arrays");
}
