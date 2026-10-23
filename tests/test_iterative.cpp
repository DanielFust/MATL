/**************************************************************************************************
 * File: test_iterative.cpp
 *
 * Description: Iterative solvers (cg, bicgstab, gmres) and preconditioners (Jacobi, ILU(0)):
 *              convergence on the 2D Laplacian and a nonsymmetric convection-diffusion matrix,
 *              agreement with a dense LU solution, fewer iterations with preconditioning, dense
 *              matrices and operator functions as A, non-convergence, initial guesses, and errors.
 **************************************************************************************************/
#include <cmath>
#include <type_traits>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using M = RowMajorMatrix<double>;
    using V = ColMajorMatrix<double>;
    using S = CSRMatrix<double>;

    // The 2D Laplacian (5-point stencil) on a k x k grid, plus c times a first-order upwind
    // convection term in the x direction (nonsymmetric for c != 0)
    S convection_diffusion(index_t k, double c = 0.0)
    {
        const index_t n = k * k;
        SparseBuilder<double> B(n, n);
        const auto node = [k](index_t i, index_t j) {return (i - 1) * k + j;};
        for (index_t i = 1; i <= k; ++i)
            for (index_t j = 1; j <= k; ++j)
            {
                B.add(node(i, j), node(i, j), 4.0 + c);
                if (i > 1) B.add(node(i, j), node(i - 1, j), -1.0);
                if (i < k) B.add(node(i, j), node(i + 1, j), -1.0);
                if (j > 1) B.add(node(i, j), node(i, j - 1), -1.0 - c);
                if (j < k) B.add(node(i, j), node(i, j + 1), -1.0);
            }
        return B.build();
    }

    V random_vector(index_t n)
    {
        V v(n, 1);
        matl_test::fill_random(v);
        return v;
    }

    // ||b - A x|| / ||b||
    double relative_residual(S const& a, V const& b, V const& x)
    {
        return norm_fro(b - a * x) / norm_fro(b);
    }
}

TEST_CASE("Iterative: CG on the 2D Laplacian, with and without preconditioning")
{
    const S A = convection_diffusion(20);                  // n = 400, symmetric positive definite
    const V b = random_vector(A.rows());
    const V exact = lu(dense(A)).solve(b);

    V x = zeros<double>(A.rows(), 1);
    auto plain = cg(A, b, x);
    CHECK(plain.converged);
    CHECK(plain.residual <= 1e-10);
    CHECK(relative_residual(A, b, x) <= 1e-10);
    CHECK_MATRIX_NEAR(x, Ref::from(exact), 1e-8);
    static_assert(std::is_same_v<decltype(plain), IterativeResult<double>>);

    V xj = zeros<double>(A.rows(), 1), xi = zeros<double>(A.rows(), 1);
    auto jacobi = cg(A, b, xj, {}, JacobiPreconditioner(A));
    auto ilu = cg(A, b, xi, {}, ILU0Preconditioner(A));
    CHECK(jacobi.converged);
    CHECK(ilu.converged);
    CHECK(ilu.iterations < plain.iterations);              // ILU(0) roughly halves the iterations here
    CHECK_MATRIX_NEAR(xi, Ref::from(exact), 1e-8);

    // A good initial guess needs few iterations; the exact solution needs none
    V warm = exact;
    auto again = cg(A, b, warm);
    CHECK(again.converged);
    CHECK_EQ(again.iterations, 0ul);
}

TEST_CASE("Iterative: BiCGSTAB and GMRES on a nonsymmetric matrix")
{
    const S A = convection_diffusion(15, 2.0);             // n = 225, nonsymmetric
    const V b = random_vector(A.rows());
    const V exact = lu(dense(A)).solve(b);

    for (int preconditioned = 0; preconditioned <= 1; ++preconditioned)
    {
        V x1 = zeros<double>(A.rows(), 1), x2 = zeros<double>(A.rows(), 1);
        auto r1 = preconditioned ? bicgstab(A, b, x1, {}, ILU0Preconditioner(A)) : bicgstab(A, b, x1);
        auto r2 = preconditioned ? gmres(A, b, x2, {}, ILU0Preconditioner(A))    : gmres(A, b, x2);
        CHECK(r1.converged);
        CHECK(r2.converged);
        CHECK(relative_residual(A, b, x1) <= 1e-10);
        CHECK(relative_residual(A, b, x2) <= 1e-10);
        CHECK_MATRIX_NEAR(x1, Ref::from(exact), 1e-8);
        CHECK_MATRIX_NEAR(x2, Ref::from(exact), 1e-8);
    }

    // GMRES with a small restart still converges, in more iterations; full GMRES (restart = n)
    // converges in at most n iterations
    V xs = zeros<double>(A.rows(), 1), xf = zeros<double>(A.rows(), 1);
    auto small = gmres(A, b, xs, {.restart = 5});
    auto full  = gmres(A, b, xf, {.restart = A.rows()});
    CHECK(small.converged);
    CHECK(full.converged);
    CHECK(full.iterations <= small.iterations);
}

TEST_CASE("Iterative: dense matrices, operator functions, and other vector types")
{
    const S A = convection_diffusion(8);
    const M D = dense(A);
    const V b = random_vector(A.rows());
    const V exact = lu(D).solve(b);

    V x = zeros<double>(A.rows(), 1);
    CHECK(cg(D, b, x).converged);                          // a dense matrix
    CHECK_MATRIX_NEAR(x, Ref::from(exact), 1e-8);

    // A matrix-free operator: out = A in, applied stencil by stencil would work the same way
    index_t applications = 0;
    const auto op = [&](auto const& in, auto& out) {out = A * in; ++applications;};
    V y = zeros<double>(A.rows(), 1);
    auto result = gmres(op, b, y);
    CHECK(result.converged);
    CHECK(applications >= result.iterations);
    CHECK_MATRIX_NEAR(y, Ref::from(exact), 1e-8);

    // x as a row-major matrix or a view, b as an expression
    M xr = zeros<double>(A.rows(), 1);
    CHECK(bicgstab(A, 2.0 * b, xr).converged);
    CHECK_MATRIX_NEAR(xr, Ref::from(exact) * 2.0, 1e-8);
    M both = zeros<double>(A.rows(), 2);
    CHECK(cg(A, b, both(all, 2)).converged);
    CHECK_MATRIX_NEAR(both(all, 2), Ref::from(exact), 1e-8);
    CHECK_EQ(both(1, 1), 0.0);

    // float, and b = 0 (the solution is x = 0)
    CSRMatrix<float> Af = sparse(Matrix<float>(D));
    ColMajorMatrix<float> bf = ColMajorMatrix<float>(b), xf = zeros<float>(A.rows(), 1);
    auto rf = cg(Af, bf, xf, {.tolerance = 1e-5});
    static_assert(std::is_same_v<decltype(rf), IterativeResult<float>>);
    CHECK(rf.converged);
    V z = ones<double>(A.rows(), 1);
    auto r0 = cg(A, zeros<double>(A.rows(), 1), z);
    CHECK(r0.converged);
    CHECK_EQ(norm_fro(z), 0.0);
}

TEST_CASE("Iterative: non-convergence is reported, not thrown")
{
    const S A = convection_diffusion(20);
    const V b = random_vector(A.rows());
    V x = zeros<double>(A.rows(), 1);
    auto result = cg(A, b, x, {.tolerance = 1e-12, .max_iterations = 3});
    CHECK(!result.converged);
    CHECK_EQ(result.iterations, 3ul);
    CHECK(result.residual > 1e-12);
    CHECK(std::fabs(result.residual - relative_residual(A, b, x)) < 1e-14);   // the true residual
}

TEST_CASE("Iterative: errors")
{
    const S A = convection_diffusion(4);
    V b = random_vector(16), x = zeros<double>(16, 1);
    CHECK_MATL_ERROR(dimension_error, (void)cg(A, zeros<double>(15, 1), x));
    CHECK_MATL_ERROR(dimension_error, (void)gmres(A, b, zeros<double>(16, 2)));
    CHECK_MATL_ERROR(dimension_error, (void)bicgstab(S(16, 15), b, x));

    // CG needs a positive definite matrix
    CHECK_MATL_ERROR(not_positive_definite_error, (void)cg(-1.0 * A, b, x));

    // Preconditioners need the diagonal
    S no_diagonal(2, 2, {0, 1, 2}, {1, 0}, {1.0, 1.0});
    CHECK_MATL_ERROR(singular_matrix_error, (void)JacobiPreconditioner(no_diagonal));
    CHECK_MATL_ERROR(singular_matrix_error, (void)ILU0Preconditioner(no_diagonal));
    CHECK_MATL_ERROR(dimension_error, (void)JacobiPreconditioner(S(2, 3)));
}

TEST_CASE("Iterative: ILU(0) of a tridiagonal matrix is its exact LU factorization")
{
    // With no fill-in possible, ILU(0) is exact, so one preconditioned iteration solves
    SparseBuilder<double> B(50, 50);
    for (index_t i = 1; i <= 50; ++i)
    {
        B.add(i, i, 2.5);
        if (i > 1)  B.add(i, i - 1, -1.0);
        if (i < 50) B.add(i, i + 1, -1.2);
    }
    const S T = B.build();
    const ILU0Preconditioner P(T);
    const V b = random_vector(50);
    V z(50, 1);
    P.apply(b, z);
    CHECK(relative_residual(T, b, z) < 1e-13);
    V x = zeros<double>(50, 1);
    auto result = gmres(T, b, x, {}, P);
    CHECK(result.converged);
    CHECK(result.iterations <= 2ul);
}
