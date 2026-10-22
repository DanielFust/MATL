/**************************************************************************************************
 * File: test_column_major.cpp
 *
 * Description: Column-major matrices (ColMajorMatrix, FixedSizeColMajorMatrix): storage layout,
 *              external (e.g. Fortran) memory, and every kind of operation with operands and
 *              targets of either storage order, compared with the reference implementation.
 **************************************************************************************************/
#include <sstream>
#include <type_traits>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using R = RowMajorMatrix<double>;
    using C = ColMajorMatrix<double>;

    template<typename MAT>
    MAT random_matrix(index_t m, index_t n)
    {
        MAT A(m, n);
        matl_test::fill_random(A);
        return A;
    }
}

/************************************* Types and layout **************************************/

TEST_CASE("Column-major: types and aliases")
{
    static_assert(std::is_same_v<RowMajorMatrix<double>, DynamicMatrix<double, StorageOrder::RowMajor>>);
    static_assert(std::is_same_v<ColMajorMatrix<double>, DynamicMatrix<double, StorageOrder::ColMajor>>);
    static_assert(std::is_same_v<Matrix<double>, RowMajorMatrix<double>>, "the default stays row-major");
    static_assert(std::is_same_v<FixedSizeColMajorMatrix<double, 2, 3>, FixedSizeMatrix<double, 2, 3, StorageOrder::ColMajor>>);
    static_assert(detail::traits<C>::storage_order == StorageOrder::ColMajor);
    static_assert(detail::traits<expr::Transpose<C>>::storage_order == StorageOrder::RowMajor);
    static_assert(MatrixType<C> && MatrixType<FixedSizeColMajorMatrix<float, 2, 2>>);
    static_assert(std::is_trivially_copyable_v<FixedSizeColMajorMatrix<double, 2, 3>>);
    CHECK(true);
}

TEST_CASE("Column-major: elements are stored column by column")
{
    C A = {{1, 2, 3},
           {4, 5, 6}};                      // initializer lists are rows in either storage order
    CHECK_EQ(A(1,2), 2.0);
    CHECK_EQ(A(2,1), 4.0);
    // The single index runs through the storage: down the first column, then the second, ...
    CHECK_EQ(A(1), 1.0); CHECK_EQ(A(2), 4.0); CHECK_EQ(A(3), 2.0);
    CHECK_EQ(A(4), 5.0); CHECK_EQ(A(5), 3.0); CHECK_EQ(A(6), 6.0);

    FixedSizeColMajorMatrix<double, 2, 3> F = {{1, 2, 3}, {4, 5, 6}};
    CHECK_EQ(F(2), 4.0);
    CHECK_EQ(F(2,3), 6.0);
}

TEST_CASE("Column-major: external memory in Fortran order")
{
    // A 2x3 array as Fortran (and LAPACK) store it: column by column
    double fortran[6] = {1, 4, 2, 5, 3, 6};
    C A(2, 3, fortran);
    CHECK_EQ(A(1,1), 1.0); CHECK_EQ(A(1,3), 3.0); CHECK_EQ(A(2,2), 5.0);

    A = A * 10.0;                           // results are written into the caller's array
    CHECK_EQ(fortran[1], 40.0);
    CHECK_EQ(fortran[5], 60.0);
}

/******************************* Operations across storage orders ***************************/

TEST_CASE("Column-major: elementwise operations with either storage order")
{
    C A = random_matrix<C>(4, 5);
    R B = random_matrix<R>(4, 5);
    Ref a = Ref::from(A), b = Ref::from(B);

    C X(4, 5);
    R Y(4, 5);
    X = A + B;                     CHECK_MATRIX_EQ(X, a + b);
    Y = A + B;                     CHECK_MATRIX_EQ(Y, a + b);
    X = A * 2.0 - transpose(transpose(B));   CHECK_MATRIX_EQ(X, a * 2.0 - b);
    Y = abs(A) + times(A, B);      CHECK_MATRIX_EQ(Y, Ref::from(R(abs(A))) + Ref::from(R(times(A, B))));

    C T(5, 4);
    T = transpose(B);              CHECK_MATRIX_EQ(T, transpose(b));      // storage-order copy
    T = transpose(A) - transpose(B);         CHECK_MATRIX_EQ(T, transpose(a) - transpose(b));
}

TEST_CASE("Column-major: copies and conversions between storage orders")
{
    C A = random_matrix<C>(3, 4);
    Ref a = Ref::from(A);

    R B(A);                        CHECK_MATRIX_EQ(B, a);       // construct from the other order
    C D(B);                        CHECK_MATRIX_EQ(D, a);
    R E(3, 4);
    E = A;                         CHECK_MATRIX_EQ(E, a);       // assign from the other order
    C F;
    F = E;                         CHECK_MATRIX_EQ(F, a);       // empty target is allocated
    C G = std::move(D);            CHECK_MATRIX_EQ(G, a);       // moves within one order

    FixedSizeColMajorMatrix<double, 3, 4> X;
    X = A;                         CHECK_MATRIX_EQ(X, a);
    X = B;                         CHECK_MATRIX_EQ(X, a);
    FixedSizeRowMajorMatrix<double, 3, 4> Y;
    Y = X;                         CHECK_MATRIX_EQ(Y, a);
}

TEST_CASE("Column-major: products with either storage order, small and blocked")
{
    for (index_t n : {3, 70})
    {
        C A = random_matrix<C>(n, n + 1);
        R B = random_matrix<R>(n + 1, n);
        Ref a = Ref::from(A), b = Ref::from(B);

        C P(n, n);
        R Q(n, n);
        P = A * B;                 CHECK_MATRIX_EQ(P, a * b);
        Q = A * B;                 CHECK_MATRIX_EQ(Q, a * b);
        P = transpose(B) * transpose(A);     CHECK_MATRIX_EQ(P, transpose(b) * transpose(a));
        Q = (A + 1.0) * B - P;     CHECK_MATRIX_EQ(Q, (a + 1.0) * b - transpose(b) * transpose(a));
    }
}

TEST_CASE("Column-major: the target may appear in the expression")
{
    C S = random_matrix<C>(5, 5);
    Ref s = Ref::from(S);
    C T(5, 5);

    T = S; T = transpose(T);       CHECK_MATRIX_EQ(T, transpose(s));
    T = S; T = T * T + T;          CHECK_MATRIX_EQ(T, s * s + s);
    T = S; T = T + transpose(T);   CHECK_MATRIX_EQ(T, s + transpose(s));
    T = S; T = -transpose(T * S);  CHECK_MATRIX_EQ(T, -1.0 * transpose(s * s));
}

TEST_CASE("Column-major: reductions and output")
{
    C A = {{1, -2}, {3, 4}};
    R B = {{1, 0}, {0, 1}};
    CHECK_EQ(sum(A), 6.0);
    CHECK_EQ(sum(A + B), 8.0);                 // mixed storage orders
    CHECK_EQ(trace(A), 5.0);
    CHECK_EQ(norm_1(A), 6.0);
    CHECK_EQ(norm_inf(A), 7.0);
    CHECK_EQ(dot(A, B), 5.0);

    std::ostringstream os;
    os << A;
    CHECK_EQ(os.str(), std::string("1  -2\n3   4"));
}

TEST_CASE("Column-major: temporaries of column-major expressions are column-major")
{
    using CE = expr::Elementwise<expr::Add, C, C>;
    static_assert(std::is_same_v<detail::result_matrix_t<CE>, C>);
    static_assert(std::is_same_v<detail::result_matrix_t<expr::Transpose<R>>, C>);
    static_assert(std::is_same_v<detail::result_matrix_t<expr::Elementwise<expr::Add, C, R>>, R>,
                  "mixed storage orders use the default (row-major)");

    C A = random_matrix<C>(20, 20), B = random_matrix<C>(20, 20);
    Ref a = Ref::from(A), b = Ref::from(B);
    C P(20, 20);
    P = (A + B) * (A - B);                    // operands evaluated into column-major temporaries
    CHECK_MATRIX_EQ(P, (a + b) * (a - b));
}
