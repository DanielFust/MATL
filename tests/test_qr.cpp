/**************************************************************************************************
 * File: test_qr.cpp
 *
 * Description: QR factorization and least squares: the factors (Q R = A, orthonormal Q,
 *              triangular R) of tall, square and wide matrices, small and blocked, in either
 *              storage order; least-squares solutions; products with Q; inputs of other kinds;
 *              in-place factorization; condition estimates; rank deficiency and errors.
 **************************************************************************************************/
#include <cmath>
#include <type_traits>
#include <utility>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using M = RowMajorMatrix<double>;
    using C = ColMajorMatrix<double>;

    template<typename MAT>
    MAT random_matrix(index_t m, index_t n)
    {
        MAT A(m, n);
        matl_test::fill_random(A);
        for (index_t i = 1; i <= std::min(m, n); ++i) A(i, i) += 0.5;     // full rank
        return A;
    }

    // Checks Q R = A, Q^T Q = I, R upper triangular, and the full Q orthogonal
    template<typename F, typename MAT>
    void check_factors(F const& factorization, MAT const& A, double tol)
    {
        const index_t m = A.rows(), n = A.cols(), k = std::min(m, n);
        M Q = factorization.Q(), R = factorization.R();
        CHECK_EQ(Q.rows(), m); CHECK_EQ(Q.cols(), k);
        CHECK_EQ(R.rows(), k); CHECK_EQ(R.cols(), n);
        CHECK_MATRIX_NEAR(Q * R, Ref::from(A), tol);
        CHECK_MATRIX_NEAR(transpose(Q) * Q, Ref::from(eye<double>(k)), tol);
        bool upper = true;
        for (index_t i = 1; i <= k; ++i)
            for (index_t j = 1; j < i; ++j) upper = upper && (R(i, j) == 0.0);
        CHECK(upper);
        M Qf = factorization.full_Q();
        CHECK_MATRIX_NEAR(transpose(Qf) * Qf, Ref::from(eye<double>(m)), tol);
        CHECK_MATRIX_NEAR(Qf(all, range(1, k)), Ref::from(Q), tol);
    }
}

TEST_CASE("QR: factors of tall, square and wide matrices, small and blocked, in either storage order")
{
    const std::pair<index_t, index_t> shapes[] = {{1, 1}, {5, 3}, {3, 5}, {4, 4}, {17, 17}, {40, 33},
                                                  {130, 70}, {70, 130}, {150, 150}, {300, 130}};
    for (auto [m, n] : shapes)
    {
        const double tol = 1e-12 * double(m + n);
        M A = random_matrix<M>(m, n);
        auto F = qr(A);
        static_assert(std::is_same_v<decltype(F), QRFactorization<double, StorageOrder::RowMajor>>);
        CHECK_EQ(F.rows(), m);
        CHECK_EQ(F.cols(), n);
        CHECK_EQ(F.tau().size(), std::min(m, n));
        check_factors(F, A, tol);

        C Ac = A;
        auto Fc = qr(Ac);
        static_assert(std::is_same_v<decltype(Fc), QRFactorization<double, StorageOrder::ColMajor>>);
        check_factors(Fc, Ac, tol);
        CHECK_MATRIX_NEAR(Fc.R(), Ref::from(F.R()), tol);
    }
}

TEST_CASE("QR: least squares")
{
    // A line through (0, 0), (1, 1), (2, 1): y = 1/6 + x/2
    M A = {{1, 0}, {1, 1}, {1, 2}}, y = {{0}, {1}, {1}};
    CHECK_MATRIX_NEAR(qr_solve(A, y), Ref::from(M{{1.0 / 6.0}, {0.5}}), 1e-15);

    // A consistent overdetermined system is solved exactly
    M T = random_matrix<M>(100, 80), x0 = random_matrix<M>(80, 3);
    M B = T * x0;
    CHECK_MATRIX_NEAR(qr(T).solve(B), Ref::from(x0), 1e-10);

    // An inconsistent one: the residual is orthogonal to the columns of A (normal equations)
    M b = random_matrix<M>(100, 2);
    auto F = qr(T);
    M x = F.solve(b);
    CHECK_MATRIX_NEAR(transpose(T) * (T * x - b), Ref(80, 2, 0.0), 1e-10);
    CHECK_MATRIX_NEAR(x, Ref::from(cholesky(transpose(T) * T).solve(transpose(T) * b)), 1e-9);

    // Square systems: the same solution as LU, in either storage order
    M S = random_matrix<M>(90, 90), c = random_matrix<M>(90, 2);
    CHECK_MATRIX_NEAR(qr_solve(S, c), Ref::from(lu_solve(S, c)), 1e-9);
    CHECK_MATRIX_NEAR(qr_solve(C(S), c), Ref::from(lu_solve(S, c)), 1e-9);
    static_assert(std::is_same_v<decltype(qr_solve(S, C(c))), C>, "the solution has the storage order of B");
}

TEST_CASE("QR: products with Q")
{
    M A = random_matrix<M>(9, 4), B = random_matrix<M>(9, 3);
    auto F = qr(A);
    M Q = F.full_Q();
    CHECK_MATRIX_NEAR(F.Qt_times(B), Ref::from(M(transpose(Q) * B)), 1e-13);
    CHECK_MATRIX_NEAR(F.Q_times(B), Ref::from(M(Q * B)), 1e-13);
    CHECK_MATRIX_NEAR(F.Q_times(F.Qt_times(B)), Ref::from(B), 1e-13);
    CHECK_MATRIX_NEAR(F.Qt_times(A), Ref::from(M(transpose(Q) * A)), 1e-13);
    CHECK_MATRIX_NEAR(F.Qt_times(A)(range(1, 4), all), Ref::from(F.R()), 1e-13);    // Q^T A = [R; 0]
}

TEST_CASE("QR: views, expressions, integer and float matrices, in place")
{
    M A = random_matrix<M>(8, 6);
    check_factors(qr(A(range(2, 7), range(1, 4))), M(A(range(2, 7), range(1, 4))), 1e-13);
    check_factors(qr(A + 1.0), M(A + 1.0), 1e-13);

    Matrix<int> I = {{1, 0}, {1, 1}, {1, 2}};
    auto Fi = qr(I);
    static_assert(std::is_same_v<decltype(Fi), QRFactorization<double, StorageOrder::RowMajor>>);
    CHECK_MATRIX_NEAR(Fi.solve(Matrix<int>{{0}, {1}, {1}}), Ref::from(M{{1.0 / 6.0}, {0.5}}), 1e-15);

    Matrix<float> Fl = {{1, 0}, {1, 1}, {1, 2}};
    auto Ff = qr(Fl);
    static_assert(std::is_same_v<decltype(Ff), QRFactorization<float, StorageOrder::RowMajor>>);
    CHECK_MATRIX_NEAR(Ff.solve(Matrix<float>{{0}, {1}, {1}}), Ref::from(M{{1.0 / 6.0}, {0.5}}), 1e-6);

    M P = random_matrix<M>(100, 70), P0 = P;
    auto Fp = qr_inplace(P);
    CHECK_EQ(Fp.factors().data(), P.data());
    CHECK_MATRIX_EQ(P, Ref::from(qr(P0).factors()));
    check_factors(Fp, P0, 1e-11);

    Matrix<double, 3, 2> G = {{3, 1}, {4, 2}, {0, 5}};
    auto Fg = qr_inplace(G);
    CHECK(std::fabs(std::fabs(G(1, 1)) - 5.0) < 1e-15);              // |R(1,1)| = norm of column 1
    check_factors(Fg, M{{3, 1}, {4, 2}, {0, 5}}, 1e-14);
}

TEST_CASE("QR: condition estimates")
{
    // R of an orthogonal matrix's QR is (+-1) on the diagonal: rcond 1
    CHECK(std::fabs(qr(eye<double>(6)).rcond() - 1.0) < 1e-15);
    M D = {{1, 0}, {0, 1e-6}, {0, 0}};
    CHECK(std::fabs(qr(D).rcond() - 1e-6) < 1e-18);

    for (index_t n : {10, 60, 100})
    {
        M A = random_matrix<M>(n + 20, n);
        auto F = qr(A);
        M R = F.R();
        const double exact = 1.0 / (norm_1(R) * norm_1(solve_upper(R, eye<double>(n))));
        CHECK(F.rcond() >= exact * (1 - 1e-6));
        CHECK(F.rcond() <= 10.0 * exact);
    }

    // Nearly rank deficient: factorized and solved, and rcond tells
    M N = {{1, 1}, {1, 1 + 1e-10}, {1, 1}};
    auto F = qr(N);
    CHECK(F.rcond() < 1e-9);
    CHECK(F.rcond() > 0.0);
}

TEST_CASE("QR: rank deficiency and errors")
{
    // A zero column: the factorization exists, but there is no unique least-squares solution
    M Z = {{1, 0}, {2, 0}, {3, 0}};
    auto F = qr(Z);
    check_factors(F, Z, 1e-15);
    M R = F.R();
    CHECK_EQ(R(2, 2), 0.0);
    CHECK_MATL_ERROR(singular_matrix_error, (void)F.solve(M{{1}, {1}, {1}}));
    #if MATL_THROWS_ON_ERROR
    std::string message;
    try { (void)F.solve(M{{1}, {1}, {1}}); }
    catch (singular_matrix_error const& e) { message = e.what(); }
    CHECK_EQ(message, std::string("qr solve: the matrix is rank deficient (R(2,2) is 0)"));
    #endif

    CHECK_MATL_ERROR(dimension_error, (void)qr(M{{1, 2, 3}, {4, 5, 6}}).solve(M{{1}, {1}}));    // underdetermined
    CHECK_MATL_ERROR(dimension_error, (void)qr(M{{1, 2}, {3, 4}}).solve(zeros<double>(3, 1)));
    CHECK_MATL_ERROR(dimension_error, (void)qr(M{{1, 2}, {3, 4}}).Qt_times(zeros<double>(3, 1)));
    CHECK_MATL_ERROR(dimension_error, (void)qr(M{{1, 2}, {3, 4}}).Q_times(zeros<double>(1, 3)));

    auto E = qr(M());                                                  // empty
    CHECK_EQ(E.R().size(), 0ul);
    CHECK_EQ(E.rcond(), 1.0);
}
