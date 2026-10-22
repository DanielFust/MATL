/**************************************************************************************************
 * File: test_construction.cpp
 *
 * Description: Constructing matrices from initializer lists, zeros and ones, and writing matrices
 *              and expressions to streams (operator<<, print).
 **************************************************************************************************/
#include <iomanip>
#include <sstream>
#include <type_traits>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    // Text written by operator<<
    template<typename X>
    std::string text(X const& x)
    {
        std::ostringstream os;
        os << x;
        return os.str();
    }
}

/********************************** Initializer lists ***************************************/

TEST_CASE("Construction: dynamic matrix from an initializer list")
{
    Matrix<double> A = {{1, 2, 3},
                        {4, 5, 6}};
    CHECK_EQ(A.rows(), 2ul);
    CHECK_EQ(A.cols(), 3ul);
    CHECK_EQ(A(1,1), 1.0); CHECK_EQ(A(1,3), 3.0);
    CHECK_EQ(A(2,1), 4.0); CHECK_EQ(A(2,3), 6.0);

    Matrix<int> I = {{7}};                          // 1x1
    CHECK_EQ(I.size(), 1ul);
    CHECK_EQ(I(1,1), 7);

    RowMajorMatrix<double> S{3, 4};                 // braces with two numbers: still the sized constructor
    CHECK_EQ(S.rows(), 3ul);
    CHECK_EQ(S.cols(), 4ul);
}

TEST_CASE("Construction: initializer-list assignment")
{
    Matrix<double> A;
    A = {{1, 2}, {3, 4}};                           // an empty matrix is allocated to fit
    CHECK_MATRIX_EQ(A, Ref::from(Matrix<double>{{1, 2}, {3, 4}}));

    const double* data = &A(1);
    A = {{5, 6}, {7, 8}};                           // same dimensions: written in place
    CHECK_EQ(&A(1), data);
    CHECK_EQ(A(2,1), 7.0);

    A = {};                                         // empties the matrix
    CHECK_EQ(A.size(), 0ul);
}

TEST_CASE("Construction: initializer lists in expressions")
{
    Matrix<double> A = {{1, 2}, {3, 4}}, R(2, 2);
    R = A * Matrix<double>{{0, 1}, {1, 0}};         // temporary matrix operand
    CHECK_EQ(R(1,1), 2.0); CHECK_EQ(R(1,2), 1.0);
    CHECK_EQ(R(2,1), 4.0); CHECK_EQ(R(2,2), 3.0);
}

TEST_CASE("Construction: fixed-size matrix from an initializer list")
{
    Matrix<double, 2, 3> F = {{1, 2, 3},
                              {4, 5, 6}};
    CHECK_EQ(F(2,3), 6.0);
    F = {{6, 5, 4}, {3, 2, 1}};
    CHECK_EQ(F(1,1), 6.0);
    CHECK_EQ(F(2,3), 1.0);
}

TEST_CASE("Construction: initializer lists of the wrong shape are reported")
{
    // Always checked, also without debug checks
    CHECK_MATL_ERROR(dimension_error, Matrix<double> A = {{1, 2}, {3}});            // ragged
    CHECK_MATL_ERROR(dimension_error, (Matrix<double, 2, 2>{{1, 2, 3}, {4, 5, 6}})); // wrong size
    CHECK_MATL_ERROR(dimension_error, (Matrix<double, 2, 2>{{1, 2}}));

    Matrix<double> B(2, 2);
    CHECK_MATL_ERROR(dimension_error, B = {{1, 2, 3}});                              // allocated 2x2

    Matrix<double, 2, 2> F;
    CHECK_MATL_ERROR(dimension_error, F = {{1, 2}, {3, 4}, {5, 6}});

    #if MATL_THROWS_ON_ERROR
    std::string message;
    try { Matrix<double> A = {{1, 2}, {3}}; }
    catch (dimension_error const& e) { message = e.what(); }
    CHECK_EQ(message, std::string("the rows of a matrix initializer list must have the same length: "
                                  "row 2 has 1 elements, row 1 has 2"));
    #endif
}

/************************************** zeros and ones ***************************************/

TEST_CASE("Construction: zeros and ones")
{
    CHECK_MATRIX_EQ(zeros<double>(2, 3), Ref(2, 3, 0.0));
    CHECK_MATRIX_EQ(zeros<double>(3),    Ref(3, 3, 0.0));
    CHECK_MATRIX_EQ(ones<double>(4, 2),  Ref(4, 2, 1.0));
    CHECK_MATRIX_EQ(ones<double>(2),     Ref(2, 2, 1.0));

    auto Z = zeros<float, 2, 3>();
    auto O = ones<int, 3>();
    static_assert(std::is_same_v<decltype(Z), FixedSizeRowMajorMatrix<float, 2, 3>>);
    static_assert(std::is_same_v<decltype(O), FixedSizeRowMajorMatrix<int, 3, 3>>);
    CHECK_MATRIX_EQ(Z, Ref(2, 3, 0.0));
    CHECK_MATRIX_EQ(O, Ref(3, 3, 1.0));

    Matrix<double> R(2, 2);
    R = 2.0 * ones<double>(2) + eye<double>(2);     // in an expression
    CHECK_MATRIX_EQ(R, Ref::from(Matrix<double>{{3, 2}, {2, 3}}));
}

/********************************** Writing to streams ***************************************/

TEST_CASE("Output: columns are right-aligned and separated by two spaces")
{
    Matrix<double> A = {{1, -2.5},
                        {10, 3}};
    CHECK_EQ(text(A), std::string(" 1  -2.5\n"
                                  "10     3"));     // no newline after the last row

    Matrix<int, 1, 3> F = {{1, 22, 333}};
    CHECK_EQ(text(F), std::string("1  22  333"));

    CHECK_EQ(text(Matrix<double>()), std::string(""));
}

TEST_CASE("Output: elements follow the stream's formatting")
{
    Matrix<double> A = {{1.0 / 3.0, 2.0}};
    std::ostringstream os;
    os << std::setprecision(3) << A;
    CHECK_EQ(os.str(), std::string("0.333  2"));

    std::ostringstream fixed;
    fixed << std::fixed << std::setprecision(2) << A;
    CHECK_EQ(fixed.str(), std::string("0.33  2.00"));
}

TEST_CASE("Output: expressions are evaluated before they are written")
{
    Matrix<double> A = {{1, 2}, {3, 4}};
    CHECK_EQ(text(transpose(A)), std::string("1  3\n2  4"));
    CHECK_EQ(text(A * A - 1.0),  std::string(" 6   9\n14  21"));

    std::ostringstream os;
    os << A << '\n' << 2.0 * A;                     // operator<< returns the stream
    CHECK_EQ(os.str(), std::string("1  2\n3  4\n2  4\n6  8"));
}

TEST_CASE("Output: print writes a line per row to standard output")
{
    Matrix<double> A = {{1, 2}, {3, 4}};
    std::ostringstream captured;
    std::streambuf* previous = std::cout.rdbuf(captured.rdbuf());
    print(A);
    print(A + A);
    std::cout.rdbuf(previous);
    #ifndef MATL_DEBUG
    CHECK_EQ(captured.str(), std::string("1  2\n3  4\n2  4\n6  8\n"));
    #else
    // The MATL_DEBUG trace also goes to std::cout
    CHECK(captured.str().find("1  2\n3  4\n") != std::string::npos);
    CHECK(captured.str().find("2  4\n6  8\n") != std::string::npos);
    #endif
}
