/**************************************************************************************************
 * File: test_errors.cpp
 *
 * Description: Error reporting (see _matl_errors_.hpp) and the debug checks enabled by matl.hpp
 *              (MATL_BOUNDS_CHECK, MATL_DIMENSION_CHECK). Errors throw matl::index_error or
 *              matl::dimension_error, or abort with MATL_ABORT_ON_ERROR (CHECK_MATL_ERROR checks
 *              either; abort checks run in a forked child process).
 *              (Dimensions known at compile time are checked at compile time, see
 *              test_types.cpp.)
 **************************************************************************************************/
#include <type_traits>
#include <vector>

#include "matl_test.hpp"

using namespace matl;

#ifdef MATL_BOUNDS_CHECK
TEST_CASE("Bounds check: RowMajorMatrix reports out-of-range indices")
{
    RowMajorMatrix<double> A(2, 3);
    A = 0.0;
    CHECK_MATL_ERROR(index_error, A(0) = 1.0);
    CHECK_MATL_ERROR(index_error, A(7) = 1.0);
    CHECK_MATL_ERROR(index_error, A(0, 1) = 1.0);
    CHECK_MATL_ERROR(index_error, A(3, 1) = 1.0);
    CHECK_MATL_ERROR(index_error, A(1, 4) = 1.0);
    RowMajorMatrix<double> const& C = A;
    CHECK_MATL_ERROR(index_error, (void)C(2, 4));
}

TEST_CASE("Bounds check: FixedSizeRowMajorMatrix reports out-of-range indices")
{
    FixedSizeRowMajorMatrix<double, 2, 3> F;
    CHECK_MATL_ERROR(index_error, (void)F(0));
    CHECK_MATL_ERROR(index_error, (void)F(7));
    CHECK_MATL_ERROR(index_error, (void)F(3, 1));
    CHECK_MATL_ERROR(index_error, (void)F(1, 4));
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
TEST_CASE("Dimension check: assigning an expression of the wrong size is reported")
{
    RowMajorMatrix<double> A(2, 2), B(3, 3);
    A = 0.0;
    B = 1.0;
    CHECK_MATL_ERROR(dimension_error, A = B + B);
    CHECK_MATL_ERROR(dimension_error, A = B * B);
}

TEST_CASE("Dimension check: elementwise operands of different sizes are reported")
{
    RowMajorMatrix<double> A(2, 2), B(3, 3), R(2, 2);
    A = 1.0;
    B = 2.0;
    CHECK_MATL_ERROR(dimension_error, R = A + B);
    CHECK_MATL_ERROR(dimension_error, R = A - B);
    CHECK_MATL_ERROR(dimension_error, R = A + B * 2.0);
    CHECK_MATL_ERROR(dimension_error, R = transpose(A - B));
    CHECK_MATL_ERROR(dimension_error, R = A + transpose(A + B));         // mismatch nested inside the expression
    CHECK_MATL_ERROR(dimension_error, Matrix<double> S = A + B);         // also when initializing a matrix

    FixedSizeRowMajorMatrix<double, 2, 2> F;
    F = 1.0;
    CHECK_MATL_ERROR(dimension_error, R = F + B);                        // fixed-size with dynamic
}

TEST_CASE("Dimension check: products with different inner dimensions are reported")
{
    RowMajorMatrix<double> A(2, 3), B(2, 3), R(2, 3), S(2, 2);
    A = 1.0;
    B = 1.0;
    CHECK_MATL_ERROR(dimension_error, S = A * B);                        // 2x3 * 2x3
    CHECK_MATL_ERROR(dimension_error, R = A * B + A);                    // nested product
    CHECK_MATL_ERROR(dimension_error, S = (A * B) * transpose(A));       // nested operand of a product

    FixedSizeRowMajorMatrix<double, 2, 2> F;
    CHECK_MATL_ERROR(dimension_error, F = A * B);                        // fixed-size target, dynamic operands
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
    CHECK_MATL_ERROR(dimension_error, R = e);
}

TEST_CASE("Dimension check: assigning a matrix of the wrong size is reported")
{
    RowMajorMatrix<double> A(2, 2), B(3, 3);
    A = 0.0;
    B = 1.0;
    CHECK_MATL_ERROR(dimension_error, A = B);
}
#endif

#ifdef MATL_DIMENSION_CHECK
TEST_CASE("Dimension check: fixed-size assignment of the wrong size is reported")
{
    FixedSizeRowMajorMatrix<double, 2, 2> F;
    RowMajorMatrix<double> D(3, 3);
    D = 1.0;
    CHECK_MATL_ERROR(dimension_error, F = D + D);
    CHECK_MATL_ERROR(dimension_error, F = D);
}

TEST_CASE("Dimension check: moves rebind matrices of any size, except over external memory")
{
    // A moved-to matrix takes the other matrix's data (no reallocation), so sizes may differ, as
    // when a std::vector of matrices of different sizes shifts its elements
    std::vector<RowMajorMatrix<double>> v;
    v.emplace_back(2, 2);
    v.emplace_back(3, 3);
    v[1] = 3.0;
    v.erase(v.begin());
    CHECK_EQ(v[0].rows(), 3ul);
    CHECK_EQ(v[0](3, 3), 3.0);

    // A matrix over external memory keeps its buffer, so the sizes must match. Move assignment is
    // noexcept, so the error terminates the program in either error mode.
    double buffer[4] = {0, 0, 0, 0};
    RowMajorMatrix<double> E(2, 2, buffer);
    CHECK_ABORTS(E = RowMajorMatrix<double>(3, 3));
}
#endif

/******************************* Exceptions and messages *************************************/

#if MATL_THROWS_ON_ERROR
namespace
{
    // The message of the MATL error the statement throws ("" if it throws none)
    template<typename F>
    std::string error_message(F&& f)
    {
        try { f(); }
        catch (std::exception const& e) { return e.what(); }
        return "";
    }
}

TEST_CASE("Errors: exception types derive from the standard exceptions")
{
    static_assert(std::is_base_of_v<std::invalid_argument, dimension_error>);
    static_assert(std::is_base_of_v<std::logic_error,      dimension_error>);
    static_assert(std::is_base_of_v<std::out_of_range,     index_error>);
    static_assert(std::is_base_of_v<std::logic_error,      index_error>);
    CHECK(true);
}

TEST_CASE("Errors: messages name the operation, the dimensions and the indices")
{
    RowMajorMatrix<double> A(2, 2), B(3, 3), R(2, 2), C(2, 3);
    A = 1.0; B = 1.0; C = 1.0;
    FixedSizeRowMajorMatrix<double, 2, 2> F;

    #ifdef MATL_DIMENSION_CHECK
    CHECK_EQ(error_message([&]{ R = A + B; }),
             std::string("incompatible matrix dimensions in elementwise operation: 2x2 and 3x3"));
    CHECK(error_message([&]{ R = C * A; }).find("matrix product") != std::string::npos);
    CHECK_EQ(error_message([&]{ R = B * 2.0; }), std::string("cannot assign a 3x3 value to a 2x2 matrix"));
    #endif
    #ifdef MATL_BOUNDS_CHECK
    CHECK_EQ(error_message([&]{ (void)C(3, 1); }),
             std::string("index (3,1) is out of range for a 2x3 matrix (indices are 1-based)"));
    CHECK_EQ(error_message([&]{ (void)C(7); }),
             std::string("index 7 is out of range for a matrix with 6 elements (indices are 1-based)"));
    CHECK_EQ(error_message([&]{ (void)F(0); }),
             std::string("index 0 is out of range for a matrix with 4 elements (indices are 1-based)"));
    #endif
    CHECK_EQ(error_message([&]{ allocate(F, 3, 3); }),
             std::string("cannot allocate a 2x2 fixed-size matrix as 3x3"));
}

TEST_CASE("Errors: a caught error leaves the program able to continue")
{
    RowMajorMatrix<double> A(2, 2), B(3, 3), R(2, 2);
    A = 1.0; B = 1.0;
    #ifdef MATL_DIMENSION_CHECK
    CHECK_THROWS_AS(R = A + B, dimension_error);
    #endif
    R = A + A;                                      // temporaries were released, evaluation works
    CHECK_EQ(R(2, 2), 2.0);
    CHECK_EQ(MemoryManager::temporary_arena().bytes_in_use(), 0ul);
}
#endif
