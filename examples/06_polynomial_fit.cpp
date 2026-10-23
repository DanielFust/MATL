/**************************************************************************************************
 * Example 6: fitting a polynomial (an application)
 *
 * Fits a quadratic y = c1 + c2 x + c3 x^2 to measured points by least squares, combining several
 * features: views to build the design matrix column by column, elementwise powers, matrix
 * products, a Cholesky solve, the condition estimate, and reductions.
 *
 * The least-squares coefficients c minimize ||V c - y||_2, where row i of the design matrix V is
 * (1, x_i, x_i^2). They solve the normal equations
 *
 *     (V^T V) c = V^T y,
 *
 * whose matrix V^T V is symmetric positive definite (when the x_i are distinct), so Cholesky
 * applies. (The normal equations square the condition number of V; for badly conditioned fits a
 * QR factorization is the better tool.)
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

    // Normal equations, solved by Cholesky
    Matrix<double> VtV = transpose(V) * V;
    Matrix<double> Vty = transpose(V) * y;
    auto C = cholesky(VtV);
    Matrix<double> c = C.solve(Vty);

    std::cout << "fitted coefficients (exact: 1, 2, -3):\n" << c << "\n\n";
    std::cout << "rcond(V^T V) = " << C.rcond() << "\n";

    // How well does the fit match? The residual r = y - V c, and its root mean square
    Matrix<double> r = y - V * c;
    std::cout << "rms residual = " << std::sqrt(mean(times(r, r))) << "\n";
    std::cout << "largest residual = " << max(abs(r)) << "\n\n";

    // Evaluate the fitted polynomial between the measurements
    Matrix<double> t = {{0.25}, {0.5}, {0.75}};
    Matrix<double> p = c(1) + c(2) * t + c(3) * times(t, t);
    std::cout << "fit at x = 0.25, 0.5, 0.75:\n" << p << "\n";

    return 0;
}
