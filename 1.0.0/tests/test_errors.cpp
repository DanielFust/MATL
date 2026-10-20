/**************************************************************************************************
 * File: test_errors.cpp
 *
 * Description: Debug checks enabled by matl.hpp (MATL_BOUNDS_CHECK, MATL_DIMENSION_CHECK).
 *              RowMajorMatrix reports errors with abort(); FixedSizeRowMajorMatrix throws
 *              std::out_of_range. Abort checks run in a forked child process.
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
