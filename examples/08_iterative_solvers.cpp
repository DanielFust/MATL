/**************************************************************************************************
 * Example 8: iterative solvers
 *
 * Solving large sparse systems with conjugate gradients (CG), BiCGSTAB and GMRES, preconditioning
 * with Jacobi and ILU(0), checking convergence, and solving with an operator that is never stored
 * as a matrix.
 *
 * The problem: the 2D Poisson equation -(u_xx + u_yy) = f on the unit square, discretized with the
 * 5-point finite-difference stencil on a k x k grid of interior points. The matrix is symmetric
 * positive definite, so CG applies. Adding a convection term makes it nonsymmetric, which needs
 * BiCGSTAB or GMRES.
 *
 * Build and run from the repository root:
 *
 *     make examples && ./build/examples/08_iterative_solvers
 **************************************************************************************************/
#include <cmath>
#include <iostream>
#include <matl.hpp>

using namespace matl;

// The 5-point stencil on a k x k grid with spacing h = 1/(k+1), times h^2, plus a first-order
// upwind convection term c u_x (nonsymmetric for c != 0)
CSRMatrix<double> poisson_2d(index_t k, double c = 0.0)
{
    const index_t n = k * k;
    const double h = 1.0 / double(k + 1);
    SparseBuilder<double> B(n, n);
    B.reserve(5 * n);
    const auto node = [k](index_t i, index_t j) {return (i - 1) * k + j;};
    for (index_t i = 1; i <= k; ++i)
        for (index_t j = 1; j <= k; ++j)
        {
            B.add(node(i, j), node(i, j), 4.0 + c * h);
            if (i > 1) B.add(node(i, j), node(i - 1, j), -1.0);
            if (i < k) B.add(node(i, j), node(i + 1, j), -1.0);
            if (j > 1) B.add(node(i, j), node(i, j - 1), -1.0 - c * h);
            if (j < k) B.add(node(i, j), node(i, j + 1), -1.0);
        }
    return B.build();
}

int main()
{
    const index_t k = 100;                                 // 10,000 unknowns, 50,000 stored entries
    const index_t n = k * k;
    const double h = 1.0 / double(k + 1);
    const CSRMatrix<double> A = poisson_2d(k);
    std::cout << "A: " << n << " x " << n << ", " << A.nnz() << " stored entries\n\n";

    // Right-hand side h^2 f for f = 1, and an initial guess of zero
    const ColMajorMatrix<double> b = (h * h) * ones<double>(n, 1);

    /******************************* CG, with and without preconditioning **********************/

    // A solver returns whether it converged, the iterations it took, and the relative residual
    // ||b - A x|| / ||b||. Not converging is not an error, so check the result.
    ColMajorMatrix<double> x = zeros<double>(n, 1);
    auto plain = cg(A, b, x);
    std::cout << "CG:                 converged " << plain.converged << " in " << plain.iterations
              << " iterations, residual " << plain.residual << "\n";

    x = 0.0;
    auto jacobi = cg(A, b, x, {}, JacobiPreconditioner(A));
    std::cout << "CG with Jacobi:     converged " << jacobi.converged << " in " << jacobi.iterations
              << " iterations (the Laplacian's diagonal is constant, so Jacobi changes nothing here)\n";

    x = 0.0;
    auto ilu = cg(A, b, x, {}, ILU0Preconditioner(A));
    std::cout << "CG with ILU(0):     converged " << ilu.converged << " in " << ilu.iterations << " iterations\n";
    std::cout << "u at the center of the square: " << x((k / 2 - 1) * k + k / 2) << "\n\n";
    const ColMajorMatrix<double> solution = x;

    // Options: the tolerance, an iteration limit, and (for GMRES) the restart length
    x = 0.0;
    auto rough = cg(A, b, x, {.tolerance = 1e-4});
    std::cout << "CG to a tolerance of 1e-4: " << rough.iterations << " iterations\n";
    x = 0.0;
    auto limited = cg(A, b, x, {.max_iterations = 10});
    // (CG minimizes the error in the A-norm, so the residual need not decrease at every step)
    std::cout << "CG limited to 10 iterations: converged " << limited.converged << ", residual "
              << limited.residual << "\n\n";

    /***************************** A nonsymmetric system: BiCGSTAB, GMRES ***********************/

    const CSRMatrix<double> C = poisson_2d(k, 50.0);       // with convection
    ColMajorMatrix<double> y = zeros<double>(n, 1);
    auto bicg = bicgstab(C, b, y, {}, ILU0Preconditioner(C));
    // (each BiCGSTAB iteration costs two products with C, each GMRES iteration one)
    std::cout << "BiCGSTAB with ILU(0):        converged " << bicg.converged << " in " << bicg.iterations << " iterations\n";
    y = 0.0;
    auto gm = gmres(C, b, y, {.restart = 30}, ILU0Preconditioner(C));
    std::cout << "GMRES(30) with ILU(0):       converged " << gm.converged << " in " << gm.iterations << " iterations\n\n";

    /*************************** An operator that is not stored as a matrix *********************/

    // Any function op(in, out) that sets out = A in can stand for A. Here it applies the 5-point
    // stencil directly, as a matrix-free PDE code would.
    const auto stencil = [k](auto const& in, auto& out)
    {
        for (index_t i = 1; i <= k; ++i)
            for (index_t j = 1; j <= k; ++j)
            {
                const index_t p = (i - 1) * k + j;
                double v = 4.0 * in(p);
                if (i > 1) v -= in(p - k);
                if (i < k) v -= in(p + k);
                if (j > 1) v -= in(p - 1);
                if (j < k) v -= in(p + 1);
                out(p) = v;
            }
    };
    ColMajorMatrix<double> u = zeros<double>(n, 1);
    auto free = cg(stencil, b, u);
    std::cout << "matrix-free CG:     converged " << free.converged << " in " << free.iterations
              << " iterations; max |u - solution| = " << max(abs(u - solution)) << "\n";

    return 0;
}
