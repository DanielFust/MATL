/**************************************************************************************************
 * File: test_structured.cpp
 *
 * Description: Sparse arithmetic (sums, differences, products, the explicit transpose, the
 *              diagonal) and the direct solvers for structured sparse matrices: tridiagonal
 *              (Thomas) and banded (BandedMatrix, banded LU and Cholesky). Results are compared
 *              with the same operations on dense matrices.
 **************************************************************************************************/
#include <cmath>
#include <tuple>
#include <type_traits>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using M = RowMajorMatrix<double>;
    using S = CSRMatrix<double>;

    M random_matrix(index_t m, index_t n)
    {
        M A(m, n);
        matl_test::fill_random(A);
        return A;
    }

    S random_sparse(index_t m, index_t n, double fill = 0.25)
    {
        M D = random_matrix(m, n);
        std::uniform_real_distribution<double> u(0.0, 1.0);
        for (index_t k = 1; k <= D.size(); ++k)
            if (u(matl_test::rng()) > fill) D(k) = 0.0;
        return sparse(D);
    }

    // A random banded matrix with kl subdiagonals and ku superdiagonals (diagonally weighted, so
    // that it is not singular)
    BandedMatrix<double> random_banded(index_t n, index_t kl, index_t ku)
    {
        BandedMatrix<double> B(n, kl, ku);
        std::uniform_int_distribution<int> d(-5, 5);
        for (index_t i = 1; i <= n; ++i)
            for (index_t j = (i > kl ? i - kl : 1); j <= std::min(n, i + ku); ++j)
                B.entry(i, j) = double(d(matl_test::rng())) + (i == j ? 0.5 : 0.0);
        return B;
    }

    // ||D x - b|| / (||D|| ||x|| + ||b||): small for a backward-stable solve, however ill-conditioned D is
    double relative_residual(M const& D, M const& x, M const& b)
    {
        return norm_fro(D * x - b) / (norm_fro(D) * norm_fro(x) + norm_fro(b));
    }

    // The 2D Laplacian on a k x k grid
    S laplacian_2d(index_t k)
    {
        SparseBuilder<double> B(k * k, k * k);
        const auto node = [k](index_t i, index_t j) {return (i - 1) * k + j;};
        for (index_t i = 1; i <= k; ++i)
            for (index_t j = 1; j <= k; ++j)
            {
                B.add(node(i, j), node(i, j), 4.0);
                if (i > 1) B.add(node(i, j), node(i - 1, j), -1.0);
                if (i < k) B.add(node(i, j), node(i + 1, j), -1.0);
                if (j > 1) B.add(node(i, j), node(i, j - 1), -1.0);
                if (j < k) B.add(node(i, j), node(i, j + 1), -1.0);
            }
        return B.build();
    }
}

/******************************************* Sparse arithmetic ***********************************/

TEST_CASE("Structured: sparse sums, differences, products and transposes")
{
    for (index_t n : {1, 6, 40})
    {
        const S A = random_sparse(n + 2, n), B = random_sparse(n + 2, n), C = random_sparse(n, n + 3);
        const Ref a = Ref::from(dense(A)), b = Ref::from(dense(B)), c = Ref::from(dense(C));

        CHECK_MATRIX_EQ(dense(A + B), a + b);
        CHECK_MATRIX_EQ(dense(A - B), a - b);
        CHECK_MATRIX_EQ(dense(A * C), a * c);
        const S At = transpose(A);
        CHECK_MATRIX_EQ(dense(At), transpose(a));
        CHECK_MATRIX_EQ(dense(transpose(A)), transpose(a));
        CHECK_MATRIX_EQ(dense(At * A + At * A), 2.0 * (transpose(a) * a));
    }

    // The result patterns: the union for sums (entries that cancel are kept), all reachable (i,j)
    // for products
    S A(2, 2, {0, 1, 2}, {0, 1}, {1.0, 2.0}), B(2, 2, {0, 1, 2}, {0, 0}, {-1.0, 3.0});
    S D = A + B;
    CHECK_EQ(D.nnz(), 3ul);                             // (1,1) cancels to 0 but is stored
    CHECK(D.is_stored(1, 1));
    CHECK_EQ(D(1, 1), 0.0);
    CHECK_EQ((A * B).nnz(), 2ul);

    // Element types promote; the diagonal
    CSRMatrix<int> I = sparse(Matrix<int>{{1, 0}, {0, 2}});
    static_assert(std::is_same_v<decltype(I + A), S>);
    CHECK_MATRIX_EQ(diagonal(A + B), Ref::from(M{{0}, {2}}));
    CHECK_MATRIX_EQ(diagonal(random_sparse(3, 5) * 0.0), Ref(3, 1, 0.0));

    CHECK_MATL_ERROR(dimension_error, (void)(A + S(2, 3)));
    CHECK_MATL_ERROR(dimension_error, (void)(A * S(3, 2)));
}

/********************************************* Tridiagonal ***************************************/

TEST_CASE("Structured: tridiagonal systems")
{
    for (index_t n : {1, 2, 9, 200})
    {
        // Diagonally dominant: the Thomas algorithm needs no pivoting
        M lower(n > 0 ? n - 1 : 0, 1), diag(n, 1), upper(n > 0 ? n - 1 : 0, 1);
        matl_test::fill_random(lower);
        matl_test::fill_random(upper);
        for (index_t i = 1; i <= n; ++i) diag(i) = 12.0 + double(i % 3);
        M D = zeros<double>(n, n);
        for (index_t i = 1; i <= n; ++i)
        {
            D(i, i) = diag(i);
            if (i > 1) D(i, i - 1) = lower(i - 1);
            if (i < n) D(i, i + 1) = upper(i);
        }
        const M b = random_matrix(n, 3);
        const M x = solve_tridiagonal(lower, diag, upper, b);
        CHECK_MATRIX_NEAR(D * x, Ref::from(b), 1e-12);
    }

    // A zero pivot (the matrix needs pivoting), and dimensions
    M one = ones<double>(1, 1), zero = zeros<double>(2, 1);
    CHECK_MATL_ERROR(singular_matrix_error, (void)solve_tridiagonal(one, zero, one, ones<double>(2, 1)));
    CHECK_MATL_ERROR(dimension_error, (void)solve_tridiagonal(one, ones<double>(3, 1), one, ones<double>(3, 1)));
    CHECK_MATL_ERROR(dimension_error, (void)solve_tridiagonal(one, zero + 1.0, one, ones<double>(3, 1)));
}

/*********************************************** Banded ******************************************/

TEST_CASE("Structured: banded matrices")
{
    BandedMatrix<double> B(5, 1, 2);
    CHECK_EQ(B.lower_bandwidth(), 1ul);
    CHECK_EQ(B.upper_bandwidth(), 2ul);
    B.entry(3, 2) = 1.0;
    B.entry(1, 3) = 2.0;
    CHECK_EQ(B(3, 2), 1.0);
    CHECK_EQ(B(5, 1), 0.0);                              // outside the band: 0
    CHECK_MATL_ERROR(index_error, B.entry(5, 1) = 1.0);  // and not writable
    CHECK_EQ(B.band().size(), 5ul * 4ul);

    // Conversions, and bandwidths found from the nonzeros
    const S L = laplacian_2d(6);
    const BandedMatrix<double> LB = banded(L);
    CHECK_EQ(LB.lower_bandwidth(), 6ul);
    CHECK_EQ(LB.upper_bandwidth(), 6ul);
    CHECK_MATRIX_EQ(dense(LB), Ref::from(dense(L)));
    CHECK_MATRIX_EQ(dense(banded(dense(L))), Ref::from(dense(L)));
    CHECK_EQ(BandedMatrix<double>(3, 7, 7).upper_bandwidth(), 2ul);    // at most n - 1

    // Products
    const M X = random_matrix(36, 2);
    CHECK_MATRIX_EQ(LB * X, Ref::from(dense(L)) * Ref::from(X));
}

TEST_CASE("Structured: banded LU and Cholesky against dense factorizations")
{
    for (auto [n, kl, ku] : {std::tuple<index_t, index_t, index_t>{1, 0, 0}, {8, 1, 1}, {30, 3, 2}, {50, 0, 4}, {50, 5, 0}, {120, 7, 9}})
    {
        const BandedMatrix<double> B = random_banded(n, kl, ku);
        const M D = dense(B);
        const M b = random_matrix(n, 2);
        auto F = lu(B);
        auto G = lu(D);
        CHECK(relative_residual(D, F.solve(b), b) < 1e-14);
        CHECK_MATRIX_EQ(F.solve(b), Ref::from(G.solve(b)));          // the same arithmetic as dense LU
        CHECK(std::fabs(F.determinant() / G.determinant() - 1.0) < 1e-9);
        CHECK(std::fabs(F.rcond() / G.rcond() - 1.0) < 1e-6);
        M X = b;
        F.solve_inplace(X(all, 2));
        CHECK_MATRIX_EQ(X(all, 2), Ref::from(M(G.solve(b)(all, 2))));
    }

    // Cholesky of the 2D Laplacian: a direct solver for a 2D PDE matrix, O(n k^2)
    const S L = laplacian_2d(15);                        // n = 225, bandwidth 15
    const M b = random_matrix(225, 1);
    auto C = cholesky(banded(L));
    auto Cd = cholesky(dense(L));
    CHECK_MATRIX_NEAR(C.solve(b), Ref::from(Cd.solve(b)), 1e-11);
    CHECK(std::fabs(C.rcond() / Cd.rcond() - 1.0) < 1e-6);
    CHECK(std::fabs(std::log(C.determinant()) - std::log(Cd.determinant())) < 1e-9);
    CHECK_MATRIX_NEAR(L * cholesky_solve(banded(L), b), Ref::from(b), 1e-10);
    CHECK_MATRIX_NEAR(L * lu_solve(banded(L), b), Ref::from(b), 1e-10);

    // Integer matrices are factorized in double
    BandedMatrix<int> Bi(3, 1, 1);
    for (index_t i = 1; i <= 3; ++i) Bi.entry(i, i) = 2;
    static_assert(std::is_same_v<decltype(lu(Bi)), BandedLUFactorization<double>>);
    CHECK_EQ(lu(Bi).determinant(), 8.0);
}

TEST_CASE("Structured: banded errors")
{
    BandedMatrix<double> Z(3, 1, 1);                     // all zero: singular
    CHECK_MATL_ERROR(singular_matrix_error, (void)lu(Z));
    BandedMatrix<double> N(2, 1, 1);
    N.entry(1, 1) = 1.0; N.entry(2, 1) = 2.0; N.entry(2, 2) = 1.0;   // not positive definite (lower band: [1; 2 1])
    CHECK_MATL_ERROR(not_positive_definite_error, (void)cholesky(N));
    CHECK_MATL_ERROR(dimension_error, (void)lu(random_banded(4, 1, 1)).solve(zeros<double>(3, 1)));
    CHECK_MATL_ERROR(dimension_error, (void)banded(S(2, 3)));
    CHECK_MATL_ERROR(dimension_error, (void)(random_banded(4, 1, 1) * zeros<double>(3, 1)));
}
