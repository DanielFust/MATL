/**************************************************************************************************
 * File: test_solvers.cpp
 *
 * Description: Linear systems: triangular solves, LU and Cholesky factorizations (factors,
 *              solutions, in-place variants, determinants, condition estimates), and their errors
 *              (singular and non-positive-definite matrices, dimensions). Solutions are checked
 *              by their residuals A X - B.
 **************************************************************************************************/
#include <cmath>
#include <limits>
#include <type_traits>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using M = RowMajorMatrix<double>;
    using C = ColMajorMatrix<double>;

    const double NaN = std::numeric_limits<double>::quiet_NaN();

    template<typename MAT>
    MAT random_matrix(index_t m, index_t n)
    {
        MAT A(m, n);
        matl_test::fill_random(A);
        return A;
    }

    // A random square matrix that needs pivoting (no dominant diagonal), and is not exactly
    // singular (the 0.5 on the diagonal cannot cancel integer elements)
    template<typename MAT>
    MAT general_matrix(index_t n)
    {
        MAT A = random_matrix<MAT>(n, n);
        for (index_t i = 1; i <= n; ++i) A(i, i) += 0.5;
        return A;
    }

    // A random symmetric positive definite matrix
    template<typename MAT>
    MAT spd_matrix(index_t n)
    {
        MAT B = random_matrix<MAT>(n, n);
        MAT A(n, n);
        A = B * transpose(B) + double(n) * eye<double>(n);
        return A;
    }

    // Reciprocal condition number from the exact inverse (for checking the estimates)
    template<typename F, typename MAT>
    double exact_rcond(F const& factorization, MAT const& A)
    {
        M inverse = factorization.solve(eye<double>(A.rows()));
        return 1.0 / (norm_1(A) * norm_1(inverse));
    }
}

/*************************************** Triangular solves ***************************************/

TEST_CASE("Solvers: triangular solves")
{
    const index_t n = 6;
    M L(n, n), U(n, n);
    L = 0.0; U = 0.0;
    for (index_t i = 1; i <= n; ++i)
        for (index_t j = 1; j <= n; ++j)
        {
            const double v = double((3 * i + 5 * j) % 7) - 3.0;
            if (i > j) L(i, j) = v;
            if (i < j) U(i, j) = v;
            if (i == j) {L(i, j) = 4.0 + double(i); U(i, j) = -3.0 - double(i);}
        }
    M B = random_matrix<M>(n, 3);
    Ref b = Ref::from(B);

    CHECK_MATRIX_NEAR(L * solve_lower(L, B), b, 1e-12);
    CHECK_MATRIX_NEAR(U * solve_upper(U, B), b, 1e-12);

    // Only the named triangle is read
    M LU = L + U;
    for (index_t i = 1; i <= n; ++i) LU(i, i) = L(i, i);
    M Lgarbage = L, Ugarbage = U;
    for (index_t i = 1; i <= n; ++i)
        for (index_t j = 1; j <= n; ++j)
        {
            if (i < j) Lgarbage(i, j) = NaN;
            if (i > j) Ugarbage(i, j) = NaN;
        }
    CHECK_MATRIX_NEAR(solve_lower(Lgarbage, B), Ref::from(solve_lower(L, B)), 0.0);
    CHECK_MATRIX_NEAR(solve_upper(Ugarbage, B), Ref::from(solve_upper(U, B)), 0.0);

    // Unit diagonal: the stored diagonal is ignored
    M L1 = L;
    for (index_t i = 1; i <= n; ++i) L1(i, i) = 1.0;
    Lgarbage = L;
    for (index_t i = 1; i <= n; ++i) Lgarbage(i, i) = NaN;
    CHECK_MATRIX_NEAR(L1 * solve_lower(Lgarbage, B, Diagonal::Unit), b, 1e-12);

    // In place, into a view; a triangle given by an expression; a column-major right-hand side
    M X = B;
    solve_upper_inplace(U, X(all, range(1, 2)));
    CHECK_MATRIX_NEAR(U * X(all, range(1, 2)), Ref::from(M(B(all, range(1, 2)))), 1e-12);
    CHECK_MATRIX_NEAR(X(all, 3), Ref::from(M(B(all, 3))), 0.0);
    CHECK_MATRIX_NEAR(L * solve_lower(L * 1.0, C(B)), b, 1e-12);
    static_assert(std::is_same_v<decltype(solve_lower(L, C(B))), C>);
    static_assert(std::is_same_v<decltype(solve_lower(Matrix<int>(2, 2), Matrix<int>(2, 1))), M>,
                  "integer systems are solved in double");
}

TEST_CASE("Solvers: triangular solve errors")
{
    M L = {{2, 0}, {1, 0}}, B = {{1}, {1}};
    CHECK_MATL_ERROR(singular_matrix_error, (void)solve_lower(L, B));
    CHECK_MATL_ERROR(dimension_error, (void)solve_lower(zeros<double>(2, 3), B));
    CHECK_MATL_ERROR(dimension_error, (void)solve_upper(M{{1, 0}, {0, 1}}, zeros<double>(3, 1)));
}

/********************************************** LU ***********************************************/

TEST_CASE("Solvers: LU of a matrix that needs pivoting")
{
    M A = {{0, 1}, {1, 1}};
    auto F = lu(A);
    CHECK_EQ(F.pivots()[0], 2ul);
    CHECK_EQ(F.pivots()[1], 2ul);
    CHECK_MATRIX_EQ(F.P(), Ref::from(M{{0, 1}, {1, 0}}));
    CHECK_MATRIX_EQ(F.L(), Ref::from(M{{1, 0}, {0, 1}}));
    CHECK_MATRIX_EQ(F.U(), Ref::from(M{{1, 1}, {0, 1}}));
    CHECK_MATRIX_EQ(F.solve(M{{2}, {3}}), Ref::from(M{{1}, {2}}));
    CHECK_EQ(F.determinant(), -1.0);
    CHECK_EQ(F.size(), 2ul);

    M S = {{2, 1, 1}, {4, 3, 3}, {8, 7, 9}};
    auto G = lu(S);
    CHECK_MATRIX_NEAR(G.L() * G.U(), Ref::from(M(G.P() * S)), 1e-14);
    CHECK(std::fabs(G.determinant() - 4.0) < 1e-12);
}

TEST_CASE("Solvers: LU solves random systems, small and blocked, in either storage order")
{
    for (index_t n : {1, 5, 15, 16, 17, 63, 64, 65, 150, 260})
    {
        M A = general_matrix<M>(n);
        M B = random_matrix<M>(n, 3);
        auto F = lu(A);
        static_assert(std::is_same_v<decltype(F), LUFactorization<double, StorageOrder::RowMajor>>);
        const double tol = 1e-10 * double(n);

        CHECK_MATRIX_NEAR(F.L() * F.U(), Ref::from(M(F.P() * A)), tol);
        CHECK_MATRIX_NEAR(A * F.solve(B), Ref::from(B), tol);
        for (index_t k = 1; k <= n; ++k)
            CHECK(F.pivots()[k - 1] >= k && F.pivots()[k - 1] <= n);

        C Ac = A;
        auto Fc = lu(Ac);
        static_assert(std::is_same_v<decltype(Fc), LUFactorization<double, StorageOrder::ColMajor>>);
        CHECK_MATRIX_NEAR(Fc.factors(), Ref::from(F.factors()), tol);
        CHECK_MATRIX_NEAR(Ac * Fc.solve(B), Ref::from(B), tol);
        CHECK_MATRIX_NEAR(A * lu_solve(A, B), Ref::from(B), tol);
    }
}

TEST_CASE("Solvers: LU of views, expressions, integer and float matrices")
{
    M A = general_matrix<M>(6);
    M b = random_matrix<M>(3, 1);
    CHECK_MATRIX_NEAR(A(range(2, 4), range(2, 4)) * lu(A(range(2, 4), range(2, 4))).solve(b), Ref::from(b), 1e-12);
    M S = A + transpose(A), R = random_matrix<M>(6, 2);
    CHECK_MATRIX_NEAR(S * lu(A + transpose(A)).solve(R), Ref::from(R), 1e-10);

    Matrix<int> I = {{2, 1}, {1, 3}};
    auto Fi = lu(I);
    static_assert(std::is_same_v<decltype(Fi), LUFactorization<double, StorageOrder::RowMajor>>);
    CHECK_EQ(Fi.determinant(), 5.0);
    CHECK_MATRIX_NEAR(Fi.solve(Matrix<int>{{3}, {4}}), Ref::from(M{{1}, {1}}), 1e-15);

    Matrix<float> Ff = {{4, 1}, {2, 3}};
    auto Ffl = lu(Ff);
    static_assert(std::is_same_v<decltype(Ffl), LUFactorization<float, StorageOrder::RowMajor>>);
    CHECK_MATRIX_NEAR(Ffl.solve(Matrix<float>{{5}, {5}}), Ref::from(M{{1}, {1}}), 1e-6);

    Matrix<double, 3, 3> F = {{2, 0, 1}, {1, 3, 0}, {0, 1, 4}};
    CHECK_MATRIX_NEAR(F * lu(F).solve(M{{3}, {4}, {5}}), Ref::from(M{{3}, {4}, {5}}), 1e-14);
}

TEST_CASE("Solvers: LU in place")
{
    M A = general_matrix<M>(70);
    M A0 = A;
    M B = random_matrix<M>(70, 2);
    auto F = lu_inplace(A);
    CHECK_EQ(F.factors().data(), A.data());                        // the factors are in A's memory
    CHECK_MATRIX_EQ(A, Ref::from(lu(A0).factors()));
    CHECK_MATRIX_NEAR(A0 * F.solve(B), Ref::from(B), 1e-9);

    M X = B;
    F.solve_inplace(X(all, 2));                                     // into a view
    CHECK_MATRIX_NEAR(A0 * X(all, 2), Ref::from(M(B(all, 2))), 1e-9);
    CHECK_MATRIX_EQ(X(all, 1), Ref::from(M(B(all, 1))));

    Matrix<double, 2, 2> G = {{0, 2}, {3, 1}};
    auto Fg = lu_inplace(G);
    CHECK_EQ(G(1, 1), 3.0);                                         // rows interchanged in G
    CHECK_MATRIX_NEAR(Fg.solve(M{{2}, {4}}), Ref::from(M{{1}, {1}}), 1e-15);
}

TEST_CASE("Solvers: LU of singular and nearly singular matrices")
{
    CHECK_MATL_ERROR(singular_matrix_error, (void)lu(M{{1, 2}, {2, 4}}));
    CHECK_MATL_ERROR(singular_matrix_error, (void)lu(zeros<double>(3, 3)));
    CHECK_MATL_ERROR(dimension_error, (void)lu(zeros<double>(2, 3)));
    CHECK_MATL_ERROR(dimension_error, (void)lu(M{{1, 0}, {0, 1}}).solve(zeros<double>(3, 1)));

    #if MATL_THROWS_ON_ERROR
    std::string message;
    try { (void)lu(M{{1, 2, 3}, {2, 4, 6}, {1, 0, 1}}); }
    catch (singular_matrix_error const& e) { message = e.what(); }
    CHECK_EQ(message, std::string("lu: the matrix is singular (column 3 has no nonzero pivot)"));
    #endif

    // Nearly singular: factorized, and rcond tells
    M N = {{1, 1}, {1, 1 + 1e-12}};
    auto F = lu(N);
    CHECK(F.rcond() < 1e-11);
    CHECK(F.rcond() > 0.0);
    CHECK_MATRIX_NEAR(N * F.solve(M{{2}, {2}}), Ref::from(M{{2}, {2}}), 1e-9);
}

/********************************** Determinants and conditions **********************************/

TEST_CASE("Solvers: condition estimates")
{
    CHECK_EQ(lu(eye<double>(5)).rcond(), 1.0);
    CHECK_EQ(cholesky(eye<double>(5)).rcond(), 1.0);
    M D = {{1, 0}, {0, 1e-6}};
    CHECK(std::fabs(lu(D).rcond() - 1e-6) < 1e-18);
    CHECK_EQ(lu(M()).rcond(), 1.0);

    // The estimate never overestimates ||A^-1||, so it is at least the exact rcond; it is usually
    // exact (it is within a factor of 3 almost always; the check allows 10)
    M H(8, 8);                                                     // Hilbert matrix (cond ~ 1e10)
    for (index_t i = 1; i <= 8; ++i)
        for (index_t j = 1; j <= 8; ++j) H(i, j) = 1.0 / double(i + j - 1);
    for (auto const& A : {H, general_matrix<M>(30), general_matrix<M>(100)})
    {
        auto F = lu(A);
        const double exact = exact_rcond(F, A);
        CHECK(F.rcond() >= exact * (1 - 1e-6));
        CHECK(F.rcond() <= 10.0 * exact);
    }
    M S = spd_matrix<M>(40);
    auto Cs = cholesky(S);
    CHECK(Cs.rcond() >= exact_rcond(Cs, S) * (1 - 1e-6));
    CHECK(Cs.rcond() <= 10.0 * exact_rcond(Cs, S));
    CHECK(lu(H).rcond() < 1e-9);
}

TEST_CASE("Solvers: determinants")
{
    CHECK_EQ(lu(M{{2, 1}, {1, 3}}).determinant(), 5.0);
    CHECK_EQ(lu(M{{0, 1, 0}, {1, 0, 0}, {0, 0, 1}}).determinant(), -1.0);   // a permutation
    CHECK_EQ(lu(M()).determinant(), 1.0);
    M S = {{4, 2, 0}, {2, 5, 1}, {0, 1, 3}};
    CHECK(std::fabs(cholesky(S).determinant() - 44.0) < 1e-12);
    CHECK(std::fabs(lu(S).determinant() - 44.0) < 1e-12);
}

/******************************************* Cholesky ********************************************/

TEST_CASE("Solvers: Cholesky factors and solves, small and blocked, in either storage order")
{
    for (index_t n : {1, 7, 16, 17, 64, 65, 150, 260})
    {
        M A = spd_matrix<M>(n);
        M B = random_matrix<M>(n, 2);
        const double tol = 1e-10 * double(n) * double(n);
        auto F = cholesky(A);
        static_assert(std::is_same_v<decltype(F), CholeskyFactorization<double, StorageOrder::RowMajor>>);
        M L = F.L();
        bool lower = true;
        for (index_t i = 1; i <= n; ++i)
            for (index_t j = i + 1; j <= n; ++j) lower = lower && (L(i, j) == 0.0);
        CHECK(lower);
        CHECK_MATRIX_NEAR(L * transpose(L), Ref::from(A), tol);
        CHECK_MATRIX_NEAR(A * F.solve(B), Ref::from(B), tol);

        C Ac = A;
        auto Fc = cholesky(Ac);
        CHECK_MATRIX_NEAR(Fc.factor(), Ref::from(L), tol);
        CHECK_MATRIX_NEAR(A * cholesky_solve(Ac, B), Ref::from(B), tol);

        // Only the lower triangle is read
        M G = A;
        for (index_t i = 1; i <= n; ++i)
            for (index_t j = i + 1; j <= n; ++j) G(i, j) = NaN;
        CHECK_MATRIX_EQ(cholesky(G).factor(), Ref::from(L));
    }
}

TEST_CASE("Solvers: Cholesky in place, and matrices that are not positive definite")
{
    M A = spd_matrix<M>(80), A0 = A;
    auto F = cholesky_inplace(A);
    CHECK_EQ(F.factor().data(), A.data());
    CHECK_MATRIX_EQ(A, Ref::from(cholesky(A0).factor()));           // A is now L
    CHECK_MATRIX_NEAR(A0 * F.solve(eye<double>(80)), Ref::from(eye<double>(80)), 1e-9);

    Matrix<double, 2, 2> G = {{4, 2}, {2, 5}};
    auto Fg = cholesky_inplace(G);
    CHECK_EQ(G(2, 1), 1.0);
    CHECK_EQ(G(1, 2), 0.0);
    CHECK_MATRIX_NEAR(Fg.solve(M{{6}, {7}}), Ref::from(M{{1}, {1}}), 1e-15);

    CHECK_MATL_ERROR(not_positive_definite_error, (void)cholesky(M{{1, 2}, {2, 1}}));
    CHECK_MATL_ERROR(not_positive_definite_error, (void)cholesky(M{{-1}}));
    CHECK_MATL_ERROR(not_positive_definite_error, (void)cholesky(M{{1, 1}, {1, 1}}));     // singular
    CHECK_MATL_ERROR(dimension_error, (void)cholesky(zeros<double>(3, 2)));
    #if MATL_THROWS_ON_ERROR
    std::string message;
    try { (void)cholesky(M{{4, 2}, {2, 1}}); }
    catch (not_positive_definite_error const& e) { message = e.what(); }
    CHECK_EQ(message, std::string("cholesky: the matrix is not positive definite (pivot 2 is not positive)"));
    #endif
}

/***************************************** Fixed-size ********************************************/

namespace
{
    // A fixed-size matrix with the elements of a dynamic one (the same dimensions)
    template<index_t N, StorageOrder ORDER = StorageOrder::RowMajor>
    FixedSizeMatrix<double, N, N, ORDER> fixed(M const& d) {return FixedSizeMatrix<double, N, N, ORDER>(d);}

    template<index_t N, StorageOrder ORDER>
    void check_fixed_lu()
    {
        const M D = general_matrix<M>(N), S = spd_matrix<M>(N), b = random_matrix<M>(N, 2);
        const auto F = fixed<N, ORDER>(D);
        const FixedSizeMatrix<double, N, 2, ORDER> bf = b;

        auto lf = lu(F);
        auto ld = lu(DynamicMatrix<double, ORDER>(D));
        static_assert(std::is_same_v<decltype(lf), FixedLUFactorization<double, N, ORDER>>);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(lf.pivots())>, std::array<index_t, N>>);
        static_assert(std::is_same_v<decltype(lf.solve(bf)), FixedSizeMatrix<double, N, 2, ORDER>>, "fixed in, fixed out");
        static_assert(std::is_same_v<decltype(lf.solve(b)), M>, "a dynamic right-hand side gives a dynamic solution");
        static_assert(std::is_same_v<decltype(lf.L()), FixedSizeMatrix<double, N, N, ORDER>>);
        static_assert(std::is_same_v<decltype(solve_lower(lf.L(), bf)), FixedSizeMatrix<double, N, 2, ORDER>>);
        CHECK_MATRIX_NEAR(solve_upper(lf.U(), solve_lower(lf.L(), lf.P() * bf, Diagonal::Unit)), Ref::from(ld.solve(b)), 1e-12);

        // The same algorithm as the dynamic factorization: the same factors, bit for bit. The
        // solves and estimates may differ in the last bit: the compiler can fuse multiply-adds
        // differently in the fixed-size code (seen with Apple Clang -O3 on ARM).
        const auto close = [](double x, double y) {return std::abs(x - y) <= 1e-13 * std::abs(y);};
        CHECK_MATRIX_EQ(lf.factors(), Ref::from(ld.factors()));
        CHECK_MATRIX_NEAR(lf.solve(bf), Ref::from(ld.solve(b)), 1e-13);
        CHECK_MATRIX_NEAR(lf.L() * lf.U(), Ref::from(M(lf.P() * D)), 1e-12);
        CHECK(close(lf.determinant(), ld.determinant()));
        CHECK(close(lf.rcond(), ld.rcond()));

        auto cf = cholesky(fixed<N, ORDER>(S));
        static_assert(std::is_same_v<decltype(cf), FixedCholeskyFactorization<double, N, ORDER>>);
        auto cd = cholesky(DynamicMatrix<double, ORDER>(S));
        CHECK_MATRIX_EQ(cf.factor(), Ref::from(cd.factor()));
        CHECK_MATRIX_NEAR(cf.solve(bf), Ref::from(cd.solve(b)), 1e-13);
        CHECK(close(cf.rcond(), cd.rcond()));
    }
}

TEST_CASE("Solvers: fixed-size LU and Cholesky factorizations")
{
    check_fixed_lu<4, StorageOrder::RowMajor>();
    check_fixed_lu<4, StorageOrder::ColMajor>();
    check_fixed_lu<16, StorageOrder::RowMajor>();
    check_fixed_lu<17, StorageOrder::RowMajor>();     // beyond the recursion base: views and products
    check_fixed_lu<24, StorageOrder::ColMajor>();

    // Fixed-size expressions are factorized as fixed-size matrices
    const Matrix<double, 3, 3> F = {{4, 1, 0}, {1, 3, 1}, {0, 1, 2}};
    static_assert(std::is_same_v<decltype(lu(transpose(F) + F)), FixedLUFactorization<double, 3, StorageOrder::RowMajor>>);
    static_assert(std::is_same_v<decltype(cholesky(F * 1.0)), FixedCholeskyFactorization<double, 3, StorageOrder::RowMajor>>);

    // A fixed-size matrix from a dynamic result
    const Matrix<double, 3, 1> x = lu(M(F)).solve(M{{5}, {5}, {3}});
    CHECK_MATRIX_NEAR(F * x, Ref::from(M{{5}, {5}, {3}}), 1e-14);
    CHECK_MATL_ERROR(singular_matrix_error, (void)lu(Matrix<double, 2, 2>{{1, 2}, {2, 4}}));
    CHECK_MATL_ERROR(not_positive_definite_error, (void)cholesky(Matrix<double, 2, 2>{{1, 2}, {2, 1}}));
}
