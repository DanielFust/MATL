/**************************************************************************************************
 * Example 3: linear systems
 *
 * Solving A x = b with the LU factorization (any square matrix) and the Cholesky factorization
 * (symmetric positive definite matrices), checking how well conditioned a matrix is, triangular
 * solves, in-place factorization, and singular matrices.
 *
 * MATL names the algorithm instead of choosing one for you: lu, cholesky, solve_lower and
 * solve_upper. Factorize once, then solve for as many right-hand sides as needed.
 *
 * Build and run from the repository root:
 *
 *     make examples && ./build/examples/03_linear_systems
 **************************************************************************************************/
#include <cmath>
#include <iostream>
#include <matl.hpp>

using namespace matl;

int main()
{
    /************************************** LU *************************************************/

    Matrix<double> A = {{ 2,  1, -1},
                        {-3, -1,  2},
                        {-2,  1,  2}};
    Matrix<double> b = {{  8},
                        {-11},
                        { -3}};                        // a column vector is an n x 1 matrix

    // LU factorization with partial pivoting: P A = L U
    auto F = lu(A);
    Matrix<double> x = F.solve(b);
    std::cout << "x =\n" << x << "\n";                // (2, 3, -1)

    // Check the solution through its residual A x - b
    std::cout << "residual ||A x - b||_inf = " << norm_inf(A * x - b) << "\n\n";

    // The same factorization solves more systems: each column of B is one right-hand side
    Matrix<double> B = {{1, 0}, {0, 1}, {0, 0}};
    std::cout << "A X = B, X =\n" << F.solve(B) << "\n\n";

    // The factors, the determinant, and an estimate of how well conditioned A is
    std::cout << "L =\n" << F.L() << "\n\nU =\n" << F.U() << "\n\nP =\n" << F.P() << "\n\n";
    std::cout << "det(A) = " << F.determinant() << ",  rcond(A) = " << F.rcond() << "\n\n";

    // For a one-off solve, lu_solve factorizes and solves in one call
    std::cout << "lu_solve(A, b) =\n" << lu_solve(A, b) << "\n\n";

    /******************************** Condition numbers *****************************************/

    // rcond() estimates 1 / (||A||_1 ||A^-1||_1): close to 1 for a well-conditioned matrix, close
    // to 0 for a nearly singular one. A solution loses about log10(1/rcond) of the ~16 digits of
    // a double. MATL factorizes nearly singular matrices and leaves the decision to you:
    Matrix<double> N = {{1, 1},
                        {1, 1 + 1e-10}};
    auto Fn = lu(N);
    std::cout << "rcond of a nearly singular matrix: " << Fn.rcond() << "\n";
    if (Fn.rcond() < 1e-8)
        std::cout << "  -> about " << int(std::log10(1.0 / Fn.rcond())) << " digits of the solution may be lost\n\n";

    /*********************************** Cholesky ***********************************************/

    // A symmetric positive definite matrix (e.g. a covariance or stiffness matrix) can be
    // factorized as S = L L^T: no pivoting, and half the work of LU. Only the lower triangle of S
    // is read.
    Matrix<double> S = {{ 4, 12, -16},
                        {12, 37, -43},
                        {-16, -43, 98}};
    auto C = cholesky(S);
    std::cout << "Cholesky factor L =\n" << C.L() << "\n\n";
    std::cout << "S x = b, x =\n" << C.solve(b) << "\n\n";

    /******************************* Triangular systems ****************************************/

    // Forward and back substitution read only the named triangle of the matrix
    Matrix<double> T = {{2, 0, 0},
                        {1, 3, 0},
                        {4, 5, 6}};
    std::cout << "T y = b, y =\n" << solve_lower(T, b) << "\n\n";
    std::cout << "T^T z = b, z =\n" << solve_upper(transpose(T), b) << "\n\n";

    /************************************* In place *********************************************/

    // For large matrices, factorize in the matrix's own memory instead of a copy: afterwards K
    // holds the factors (and must outlive the factorization), and the solution overwrites the
    // right-hand side.
    Matrix<double> K = S;
    Matrix<double> y = b;
    auto Ck = cholesky_inplace(K);                    // K now holds L
    Ck.solve_inplace(y);                              // y now holds the solution of S y = b
    std::cout << "in place: y =\n" << y << "\n\n";

    /******************************** Singular matrices *****************************************/

    // A matrix that is exactly singular cannot be factorized; MATL reports it with an exception
    Matrix<double> singular = {{1, 2},
                               {2, 4}};
    try
    {
        auto Fs = lu(singular);
    }
    catch (singular_matrix_error const& e)
    {
        std::cout << "caught: " << e.what() << "\n";
    }

    // Likewise, Cholesky reports a matrix that is not positive definite
    try
    {
        auto Cs = cholesky(Matrix<double>{{1, 2}, {2, 1}});
    }
    catch (not_positive_definite_error const& e)
    {
        std::cout << "caught: " << e.what() << "\n";
    }

    return 0;
}
