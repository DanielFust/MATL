/**************************************************************************************************
 * File: test_errors.cpp
 *
 * Description: Debug checks enabled by matl.hpp (MATL_BOUNDS_CHECK, MATL_DIMENSION_CHECK).
 *              RowMajorMatrix and expression dimension errors abort(); FixedSizeRowMajorMatrix
 *              throws std::out_of_range. Abort checks run in a forked child process.
 *              (Dimensions known at compile time are checked at compile time, see
 *              test_types.cpp.)
 **************************************************************************************************/
#include "matl_test.hpp"

using namespace matl;

#ifdef MATL_BOUNDS_CHECK
TEST_CASE("Bounds check: RowMajorMatrix aborts on out-of-range indices")
{
    RowMajorMatrix<double> A(2, 3);
    A = 0.0;
    CHECK_ABORTS(A(0) = 1.0);
    CHECK_ABORTS(A(7) = 1.0);
    CHECK_ABORTS(A(0, 1) = 1.0);
    CHECK_ABORTS(A(3, 1) = 1.0);
    CHECK_ABORTS(A(1, 4) = 1.0);
    RowMajorMatrix<double> const& C = A;
    CHECK_ABORTS((void)C(2, 4));
}

TEST_CASE("Bounds check: FixedSizeRowMajorMatrix throws on out-of-range indices")
{
    FixedSizeRowMajorMatrix<double, 2, 3> F;
    CHECK_THROWS_AS(F(0), std::out_of_range);
    CHECK_THROWS_AS(F(7), std::out_of_range);
    CHECK_THROWS_AS(F(3, 1), std::out_of_range);
    CHECK_THROWS_AS(F(1, 4), std::out_of_range);
}

TEST_CASE("Bounds check: valid indices do not trigger the check")
{
    RowMajorMatrix<double> A(2, 3);
    A(1, 1) = 1.0;
    A(2, 3) = 1.0;
    A(6)    = 1.0;
    CHECK_EQ(A(2, 3), 1.0);
}
#endif

#ifdef MATL_DIMENSION_CHECK
TEST_CASE("Dimension check: assigning an expression of the wrong size aborts")
{
    RowMajorMatrix<double> A(2, 2), B(3, 3);
    A = 0.0;
    B = 1.0;
    CHECK_ABORTS(A = B + B);
    CHECK_ABORTS(A = B * B);
}

TEST_CASE("Dimension check: elementwise operands of different sizes abort")
{
    RowMajorMatrix<double> A(2, 2), B(3, 3), R(2, 2);
    A = 1.0;
    B = 2.0;
    CHECK_ABORTS(R = A + B);
    CHECK_ABORTS(R = A - B);
    CHECK_ABORTS(R = A + B * 2.0);
    CHECK_ABORTS(R = transpose(A - B));
    CHECK_ABORTS(R = A + transpose(A + B));         // mismatch nested inside the expression
    CHECK_ABORTS(Matrix<double> S = A + B);         // also when initializing a matrix

    FixedSizeRowMajorMatrix<double, 2, 2> F;
    F = 1.0;
    CHECK_ABORTS(R = F + B);                        // fixed-size with dynamic
}

TEST_CASE("Dimension check: products with different inner dimensions abort")
{
    RowMajorMatrix<double> A(2, 3), B(2, 3), R(2, 3), S(2, 2);
    A = 1.0;
    B = 1.0;
    CHECK_ABORTS(S = A * B);                        // 2x3 * 2x3
    CHECK_ABORTS(R = A * B + A);                    // nested product
    CHECK_ABORTS(S = (A * B) * transpose(A));       // nested operand of a product

    FixedSizeRowMajorMatrix<double, 2, 2> F;
    CHECK_ABORTS(F = A * B);                        // fixed-size target, dynamic operands
}

TEST_CASE("Dimension check: valid dimensions pass")
{
    RowMajorMatrix<double> A(2, 3), B(3, 2), R(2, 3), S(2, 2);
    A = 1.0;
    B = 2.0;
    R = A + transpose(B);
    S = A * B + 1.0;
    S = (A + 1.0) * (B * 2.0) - S;
    R = 2.0 - A;
    CHECK_EQ(S(1,1), 24.0 - 7.0);                   // (2 * 4) * 3 - (1*2*3 + 1)
}

TEST_CASE("Dimension check: a stored expression is checked when it is assigned")
{
    RowMajorMatrix<double> A(2, 2), B(2, 2), R(2, 2);
    A = 1.0;
    B = 2.0;
    auto e = A + B;
    R = e;                                          // dimensions fit
    allocate(B, 3, 3);                              // B changes size after e was built
    B = 2.0;
    CHECK_ABORTS(R = e);
}

TEST_CASE("Dimension check: assigning a matrix of the wrong size aborts")
{
    RowMajorMatrix<double> A(2, 2), B(3, 3);
    A = 0.0;
    B = 1.0;
    CHECK_ABORTS(A = B);
}
#endif

#ifdef MATL_BOUNDS_CHECK
TEST_CASE("Dimension check: fixed-size assignment of the wrong size throws")
{
    FixedSizeRowMajorMatrix<double, 2, 2> F;
    RowMajorMatrix<double> D(3, 3);
    D = 1.0;
    CHECK_THROWS_AS(F = D + D, std::out_of_range);
}
#endif
