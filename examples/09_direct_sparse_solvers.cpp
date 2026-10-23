/**************************************************************************************************
 * Example 9: direct solvers for structured sparse matrices, and sparse arithmetic
 *
 * Tridiagonal and banded systems are solved directly, in time proportional to n times the
 * bandwidth squared, with no fill-in outside the band:
 *
 *   - the 1D heat equation u_t = u_xx, stepped implicitly (backward Euler): each time step solves
 *     a tridiagonal system (solve_tridiagonal, the Thomas algorithm);
 *   - the 2D Poisson equation of example 8, whose matrix has bandwidth k on a k x k grid: banded
 *     Cholesky solves it directly, compared here with CG.
 *
 * It also shows sparse arithmetic (sums, products, transposes of sparse matrices), used here to
 * build the 2D matrix from 1D ones.
 *
 * Build and run from the repository root:
 *
 *     make examples && ./build/examples/09_direct_sparse_solvers
 **************************************************************************************************/
#include <chrono>
#include <cmath>
#include <iostream>
#include <matl.hpp>

using namespace matl;

int main()
{
    /******************************** 1D heat equation ******************************************/

    // u_t = u_xx on [0, 1], u = 0 at both ends, u(x, 0) = sin(pi x). The exact solution decays
    // as exp(-pi^2 t). Backward Euler: (I + dt/h^2 T) u_new = u_old, with T = tridiag(-1, 2, -1).
    const index_t n = 99;                                  // interior points
    const double h = 1.0 / double(n + 1), dt = 1e-3, r = dt / (h * h);
    const double pi = 3.14159265358979323846;

    ColMajorMatrix<double> lower(n - 1, 1), diagonal(n, 1), upper(n - 1, 1), u(n, 1);
    lower = -r;
    upper = -r;
    diagonal = 1.0 + 2.0 * r;
    for (index_t i = 1; i <= n; ++i) u(i) = std::sin(pi * double(i) * h);

    const index_t steps = 100;                             // to t = 0.1
    for (index_t step = 1; step <= steps; ++step)
        u = solve_tridiagonal(lower, diagonal, upper, u);  // O(n) per step

    const double t = double(steps) * dt;
    std::cout << "heat equation at t = " << t << ": u(0.5) = " << u((n + 1) / 2)
              << ", exact " << std::exp(-pi * pi * t) << " (backward Euler is first order in dt)\n\n";

    /***************************** Sparse arithmetic: a 2D matrix *******************************/

    // The 2D Laplacian is the sum of 1D Laplacians in x and y: A = T (x) I + I (x) T (Kronecker
    // products). Here it is built from sparse products of "selection" matrices instead: the 1D
    // operator acting along rows, plus the same along columns via a permutation P (the transpose
    // of the grid): A = Ax + P^T Ax P.
    const index_t k = 60;                                  // a 60 x 60 grid: 3600 unknowns
    SparseBuilder<double> bx(k * k, k * k), bp(k * k, k * k);
    for (index_t i = 1; i <= k; ++i)
        for (index_t j = 1; j <= k; ++j)
        {
            const index_t p = (i - 1) * k + j;
            bx.add(p, p, 2.0);                             // the 1D stencil along j
            if (j > 1) bx.add(p, p - 1, -1.0);
            if (j < k) bx.add(p, p + 1, -1.0);
            bp.add(p, (j - 1) * k + i, 1.0);               // P: grid point (i,j) -> (j,i)
        }
    const CSRMatrix<double> Ax = bx.build(), P = bp.build();
    const CSRMatrix<double> Pt = transpose(P);
    const CSRMatrix<double> A = Ax + Pt * Ax * P;          // the 5-point Laplacian
    std::cout << "2D Laplacian from sparse arithmetic: " << A.rows() << " x " << A.cols() << ", "
              << A.nnz() << " stored entries, A(2,2) = " << A(2, 2) << "\n";

    /************************** 2D Poisson: banded Cholesky versus CG ***************************/

    const double hk = 1.0 / double(k + 1);
    const ColMajorMatrix<double> b = (hk * hk) * ones<double>(k * k, 1);

    const auto start = std::chrono::steady_clock::now();
    const BandedMatrix<double> AB = banded(A);              // bandwidth k: nonzeros within k of the diagonal
    auto C = cholesky(AB);                                  // O(n k^2), the fill-in stays in the band
    const ColMajorMatrix<double> x_direct = C.solve(b);
    const auto middle = std::chrono::steady_clock::now();
    ColMajorMatrix<double> x_cg = zeros<double>(k * k, 1);
    auto result = cg(A, b, x_cg, {.tolerance = 1e-12}, ILU0Preconditioner(A));
    const auto end = std::chrono::steady_clock::now();

    const auto ms = [](auto d) {return std::chrono::duration<double, std::milli>(d).count();};
    std::cout << "banded Cholesky (bandwidth " << AB.lower_bandwidth() << "): " << ms(middle - start) << " ms\n";
    std::cout << "CG with ILU(0): " << result.iterations << " iterations, " << ms(end - middle) << " ms\n";
    std::cout << "max |x_direct - x_cg| = " << max(abs(x_direct - x_cg))
              << ", u at the center = " << x_direct((k / 2 - 1) * k + k / 2) << "\n";
    std::cout << "rcond(A) = " << C.rcond() << "\n";

    return 0;
}
