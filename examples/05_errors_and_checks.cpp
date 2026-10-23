/**************************************************************************************************
 * Example 5: errors and debug checks
 *
 * What MATL checks, when, and how errors are reported.
 *
 *   - Dimensions known at compile time (fixed-size matrices) are checked by the compiler.
 *   - Debug checks run in debug builds: MATL_DIMENSION_CHECK (operands and assignments fit) and
 *     MATL_BOUNDS_CHECK (indices are in range). They are on unless NDEBUG is defined (as in
 *     typical release builds) or MATL_NO_DEBUG is defined.
 *   - A few checks always run, e.g. the solvers' (singular matrices, dimensions).
 *
 * Errors throw matl::dimension_error, matl::index_error, matl::singular_matrix_error or
 * matl::not_positive_definite_error. Define MATL_ABORT_ON_ERROR to print the message and abort
 * instead, which stops a debugger at the faulty statement.
 *
 * Build and run from the repository root:
 *
 *     make examples && ./build/examples/05_errors_and_checks
 *
 * Compile it with -DNDEBUG to see which checks remain in a release build.
 **************************************************************************************************/
#include <iostream>
#include <matl.hpp>

using namespace matl;

int main()
{
    Matrix<double> A(2, 2), B(3, 3), R(2, 2);
    A = 1.0;
    B = 2.0;

    /**************************** Compile-time dimension checks *********************************/

    // Fixed-size matrices of the wrong dimensions do not compile. Uncomment to see the error:
    //
    //     Matrix<double, 2, 3> X;
    //     Matrix<double, 3, 2> Y;
    //     auto Z = X + Y;                  // error: no matching operator+

    /*********************************** Debug checks *******************************************/

    #ifdef MATL_DIMENSION_CHECK
    std::cout << "debug checks are on\n";
    try
    {
        R = A + B;                                     // 2x2 + 3x3
    }
    catch (dimension_error const& e)
    {
        std::cout << "caught dimension_error: " << e.what() << "\n";
    }
    #else
    std::cout << "debug checks are off (NDEBUG or MATL_NO_DEBUG): R = A + B would not be checked\n";
    #endif

    #ifdef MATL_BOUNDS_CHECK
    try
    {
        std::cout << A(3, 1) << "\n";                  // A has 2 rows
    }
    catch (index_error const& e)
    {
        std::cout << "caught index_error: " << e.what() << "\n";
    }

    try
    {
        auto rows = A(range(1, 3), all);               // a view past the last row
        std::cout << rows << "\n";
    }
    catch (index_error const& e)
    {
        std::cout << "caught index_error: " << e.what() << "\n";
    }
    #endif

    /******************************** Checks that always run ************************************/

    try
    {
        auto F = lu(Matrix<double>{{1, 2}, {2, 4}});   // singular
    }
    catch (singular_matrix_error const& e)
    {
        std::cout << "caught singular_matrix_error: " << e.what() << "\n";
    }

    try
    {
        Matrix<double> x = lu_solve(A + eye<double>(2), B);   // 2x2 system, 3x3 right-hand side
    }
    catch (dimension_error const& e)
    {
        std::cout << "caught dimension_error: " << e.what() << "\n";
    }

    // Every MATL exception derives from a standard one, so one handler can catch them all
    try
    {
        auto F = cholesky(-eye<double>(2));
    }
    catch (std::exception const& e)
    {
        std::cout << "caught std::exception: " << e.what() << "\n";
    }

    // After a caught error MATL is in a usable state (the target of a failed assignment may be
    // partly written)
    R = A * 2.0;
    std::cout << "\nstill working: R =\n" << R << "\n";

    return 0;
}
