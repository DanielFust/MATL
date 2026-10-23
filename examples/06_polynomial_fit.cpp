/**************************************************************************************************
 * Example 6: fitting a polynomial (an application)
 *
 * Fits a quadratic y = c1 + c2 x + c3 x^2 to measured points by least squares, combining several
 * features: views to build the design matrix column by column, elementwise powers, the QR
 * factorization, matrix products, a Cholesky solve, condition estimates and reductions.
 *
 * The least-squares coefficients c minimize ||V c - y||_2, where row i of the design matrix V is
 * (1, x_i, x_i^2). There are two standard ways to compute them:
 *
 *   - QR: factorize V = Q R, then solve R c = (Q^T y)(1:3). This is what qr(V).solve(y) does.
 *   - The normal equations (V^T V) c = V^T y. V^T V is symmetric positive definite (when the x_i
 *     are distinct), so Cholesky applies. This is cheaper, but V^T V has the square of V's
 *     condition number, so it loses twice as many digits.
 *
 * Build and run from the repository root:
 *
 *     make examples && ./build/examples/06_polynomial_fit
 **************************************************************************************************/
#include <cmath>
#include <iostream>
#include <matl.hpp>

using namespace matl;

int main()
{
    // Measurements: y = 1 + 2x - 3x^2, with a small, deterministic perturbation as "noise"
    const index_t n = 11;
    Matrix<double> x(n, 1), y(n, 1);
    for (index_t i = 1; i <= n; ++i)
    {
        x(i) = double(i - 1) / double(n - 1);                       // 0, 0.1, ..., 1
        y(i) = 1.0 + 2.0 * x(i) - 3.0 * x(i) * x(i) + 0.01 * std::sin(7.0 * double(i));
    }

    // Design matrix: column k holds x^(k-1). Each column is filled through a view.
    const index_t degree = 2;
    Matrix<double> V(n, degree + 1);
    for (index_t k = 1; k <= degree + 1; ++k)
        V(all, k) = pow(x, double(k - 1));

    /******************************** Least squares by QR ***************************************/

    auto F = qr(V);
    Matrix<double> c = F.solve(y);
    std::cout << "fitted coefficients (exact: 1, 2, -3):\n" << c << "\n\n";

    // How well does the fit match? The residual r = y - V c, and its root mean square
    Matrix<double> r = y - V * c;
    std::cout << "rms residual = " << std::sqrt(mean(times(r, r))) << "\n";
    std::cout << "largest residual = " << max(abs(r)) << "\n\n";

    // At the least-squares solution, the residual is orthogonal to the columns of V
    std::cout << "||V^T r|| = " << norm_fro(transpose(V) * r) << " (zero up to rounding)\n\n";

    /******************************** The normal equations **************************************/

    Matrix<double> VtV = transpose(V) * V;
    auto C = cholesky(VtV);
    Matrix<double> c_normal = C.solve(transpose(V) * y);
    std::cout << "normal equations give the same fit here: max difference "
              << max(abs(c_normal - c)) << "\n";

    // The condition estimates show why QR is preferred: rcond of V^T V is about the square of
    // that of V (R from the QR factorization has the condition number of V)
    std::cout << "rcond(V) ~ rcond(R) = " << F.rcond() << ",  rcond(V^T V) = " << C.rcond() << "\n\n";

    // Evaluate the fitted polynomial between the measurements
    Matrix<double> t = {{0.25}, {0.5}, {0.75}};
    Matrix<double> p = c(1) + c(2) * t + c(3) * times(t, t);
    std::cout << "fit at x = 0.25, 0.5, 0.75:\n" << p << "\n";

    return 0;
}
