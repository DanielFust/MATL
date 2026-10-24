/**************************************************************************************************
 * File: test_eigen_svd.cpp
 *
 * Description: Determinants and inverses (closed forms for fixed sizes up to 3, LU otherwise),
 *              symmetric eigenvalues and eigenvectors (eig_sym, eigenvalues_sym), and the
 *              singular value decomposition with the pseudo-inverse (svd, singular_values, pinv,
 *              rank). Decompositions are checked by their defining properties: A V = V diag(lambda)
 *              and orthonormal V; U diag(sigma) V^T = A with orthonormal U and V; the Penrose
 *              conditions for the pseudo-inverse.
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

    // A random symmetric matrix
    template<typename MAT>
    MAT symmetric_matrix(index_t n)
    {
        const MAT A = random_matrix<MAT>(n, n);
        return MAT(A + transpose(A));
    }

    // X with column j multiplied by d(j) (X diag(d))
    template<typename X, typename D>
    M scale_columns(X const& x, D const& d)
    {
        M y(x);
        for (index_t j = 1; j <= y.cols(); ++j)
            for (index_t i = 1; i <= y.rows(); ++i) y(i,j) *= d(j,1);
        return y;
    }

    // ||Q^T Q - I||_F (0 for orthonormal columns)
    template<typename Q>
    double orthonormality_error(Q const& q)
    {return norm_fro(M(transpose(q) * q) - eye<double>(q.cols()));}

    template<typename D>
    bool ascending(D const& d)
    {
        for (index_t k = 1; k < d.rows(); ++k) if (d(k,1) > d(k + 1,1)) return false;
        return true;
    }

    template<typename D>
    bool descending_nonnegative(D const& d)
    {
        for (index_t k = 1; k <= d.rows(); ++k) if (d(k,1) < 0) return false;
        for (index_t k = 1; k < d.rows(); ++k) if (d(k,1) < d(k + 1,1)) return false;
        return true;
    }

    // Checks an eigen-decomposition of the symmetric matrix a
    template<typename A, typename E>
    void check_eigen(A const& a, E const& e)
    {
        const double scale = std::max(1.0, norm_fro(a));
        CHECK(norm_fro(M(a * e.vectors()) - scale_columns(e.vectors(), e.values())) <= 1e-13 * scale);
        CHECK(orthonormality_error(e.vectors()) <= 1e-13);
        CHECK(ascending(e.values()));
        CHECK(std::abs(sum(e.values()) - trace(a)) <= 1e-12 * scale);
    }

    // Checks a singular value decomposition of a
    template<typename A, typename S>
    void check_svd(A const& a, S const& s)
    {
        const index_t k = std::min(a.rows(), a.cols());
        const double scale = std::max(1.0, norm_fro(a));
        CHECK_EQ(s.U().rows(), a.rows());  CHECK_EQ(s.U().cols(), k);
        CHECK_EQ(s.V().rows(), a.cols());  CHECK_EQ(s.V().cols(), k);
        CHECK(orthonormality_error(s.U()) <= 1e-13);
        CHECK(orthonormality_error(s.V()) <= 1e-13);
        CHECK(descending_nonnegative(s.singular_values()));
        CHECK(norm_fro(M(s.U() * s.S() * transpose(s.V())) - a) <= 1e-13 * scale);
    }

    // The Penrose conditions: A X A = A, X A X = X, A X and X A symmetric
    template<typename A, typename X>
    void check_pseudo_inverse(A const& a, X const& x)
    {
        const double scale = std::max(1.0, norm_fro(a)) * std::max(1.0, norm_fro(x));
        const M ax = a * x, xa = x * a;
        CHECK(norm_fro(M(a * x * a) - a) <= 1e-12 * scale * std::max(1.0, norm_fro(a)));
        CHECK(norm_fro(M(x * a * x) - x) <= 1e-12 * scale * std::max(1.0, norm_fro(x)));
        CHECK(norm_fro(ax - transpose(ax)) <= 1e-12 * scale);
        CHECK(norm_fro(xa - transpose(xa)) <= 1e-12 * scale);
    }
}

/*************************************** det and inv *********************************************/

TEST_CASE("Inverse: det and inv of fixed-size matrices of order 1 to 3 (closed forms)")
{
    const Matrix<double, 1, 1> A1 = {{4}};
    const Matrix<double, 2, 2> A2 = {{4, 3}, {6, 3}};
    const Matrix<double, 3, 3> A3 = {{2, 1, 0}, {0, 3, 1}, {1, 0, 4}};
    CHECK_EQ(det(A1), 4.0);
    CHECK_EQ(det(A2), -6.0);
    CHECK_EQ(det(A3), 25.0);
    CHECK(std::abs(det(A3) - lu(A3).determinant()) <= 1e-13);

    static_assert(std::is_same_v<decltype(inv(A3)), Matrix<double, 3, 3>>);
    static_assert(std::is_same_v<decltype(inv(FixedSizeColMajorMatrix<double, 2, 2>{})), FixedSizeColMajorMatrix<double, 2, 2>>);
    static_assert(std::is_same_v<decltype(det(Matrix<int, 2, 2>{})), double>, "integer matrices give double results");
    static_assert(std::is_same_v<decltype(inv(Matrix<float, 3, 3>{})), Matrix<float, 3, 3>>);

    CHECK_MATRIX_NEAR(inv(A1), Ref::from(M{{0.25}}), 0.0);
    CHECK_MATRIX_NEAR(inv(A2), Ref::from(M{{-0.5, 0.5}, {1.0, -2.0 / 3.0}}), 1e-15);
    CHECK_MATRIX_NEAR(A3 * inv(A3), Ref::from(eye<double>(3)), 1e-15);
    CHECK_MATRIX_NEAR(inv(A3), Ref::from(lu(A3).solve(eye<double, 3, 3>())), 1e-15);
    CHECK_MATRIX_NEAR(inv(FixedSizeColMajorMatrix<double, 3, 3>(A3)), Ref::from(inv(A3)), 0.0);

    // Expressions and integer matrices
    CHECK_EQ(det(A3 + transpose(A3)), det(Matrix<double, 3, 3>(A3 + transpose(A3))));
    CHECK_EQ(det(Matrix<int, 2, 2>{{1, 2}, {3, 4}}), -2.0);
    CHECK_MATRIX_NEAR(inv(Matrix<int, 2, 2>{{1, 2}, {3, 4}}), Ref::from(M{{-2, 1}, {1.5, -0.5}}), 1e-15);

    // Singular matrices: det is 0, inv reports the error
    const Matrix<double, 3, 3> S3 = {{1, 2, 3}, {2, 4, 6}, {1, 0, 1}};
    CHECK_EQ(det(S3), 0.0);
    CHECK_MATL_ERROR(singular_matrix_error, (void)inv(S3));
    CHECK_MATL_ERROR(singular_matrix_error, (void)inv(Matrix<double, 2, 2>{{1, 2}, {2, 4}}));
    CHECK_MATL_ERROR(singular_matrix_error, (void)inv(Matrix<double, 1, 1>{{0}}));
}

TEST_CASE("Inverse: det and inv of larger fixed-size and dynamic matrices (LU)")
{
    for (index_t n : {1, 2, 3, 4, 7, 16, 17, 40})
    {
        const M A = random_matrix<M>(n, n) + 10.0 * eye<double>(n);
        CHECK(std::abs(det(A) - lu(A).determinant()) <= 1e-12 * std::abs(lu(A).determinant()));
        CHECK_MATRIX_NEAR(A * inv(A), Ref::from(eye<double>(n)), 1e-12);
        const C Ac = A;
        CHECK_EQ(det(Ac), lu(Ac).determinant());
        static_assert(std::is_same_v<decltype(inv(Ac)), C>);
        CHECK_MATRIX_NEAR(Ac * inv(Ac), Ref::from(eye<double>(n)), 1e-12);
    }

    // Fixed-size matrices beyond 3: LU, fixed-size results
    Matrix<double, 5, 5> F5;
    matl_test::fill_random(F5);
    F5 = F5 + 10.0 * eye<double, 5, 5>();
    static_assert(std::is_same_v<decltype(inv(F5)), Matrix<double, 5, 5>>);
    CHECK_EQ(det(F5), lu(F5).determinant());
    CHECK_MATRIX_NEAR(F5 * inv(F5), Ref::from(eye<double>(5)), 1e-13);

    // Views and expressions
    M B = random_matrix<M>(6, 6) + 10.0 * eye<double>(6);
    CHECK_EQ(det(B(range(2, 4), range(2, 4))), det(M(B(range(2, 4), range(2, 4)))));
    CHECK_MATRIX_NEAR(inv(2.0 * B) * B, Ref::from(M(0.5 * eye<double>(6))), 1e-13);

    // Singular matrices: det is exactly 0 (also when the zero pivot comes late, in the recursion)
    M S = random_matrix<M>(30, 30) + 10.0 * eye<double>(30);
    for (index_t i = 1; i <= 30; ++i) S(i, 20) = 0.0;
    CHECK_EQ(det(S), 0.0);
    CHECK_EQ(det(zeros<double>(4, 4)), 0.0);
    CHECK_MATL_ERROR(singular_matrix_error, (void)inv(S));
    CHECK_EQ(det(M(0, 0)), 1.0);

    CHECK_MATL_ERROR(dimension_error, (void)det(M(2, 3)));
    CHECK_MATL_ERROR(dimension_error, (void)inv(M(3, 2)));
}

/**************************************** eig_sym ************************************************/

TEST_CASE("Eigen: eig_sym of random symmetric matrices, small and larger, in either storage order")
{
    for (index_t n : {1, 2, 3, 5, 16, 40})
    {
        const M A = symmetric_matrix<M>(n);
        const auto E = eig_sym(A);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(E)>, SymmetricEigen<double, StorageOrder::RowMajor>>);
        CHECK_EQ(E.size(), n);
        check_eigen(A, E);
        CHECK_MATRIX_EQ(eigenvalues_sym(A), Ref::from(E.values()));     // the same rotations

        const C Ac = A;
        const auto Ec = eig_sym(Ac);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(Ec.vectors())>, C>);
        check_eigen(Ac, Ec);
        CHECK_MATRIX_EQ(Ec.values(), Ref::from(E.values()));             // the same arithmetic
    }
}

TEST_CASE("Eigen: known spectra, the lower triangle, and repeated eigenvalues")
{
    // Only the lower triangle is read
    const M A = {{2, NaN, NaN}, {1, 2, NaN}, {1, 1, 2}};
    const M Afull = {{2, 1, 1}, {1, 2, 1}, {1, 1, 2}};
    const auto E = eig_sym(A);
    CHECK_MATRIX_NEAR(E.values(), Ref::from(M{{1}, {1}, {4}}), 1e-14);  // a repeated eigenvalue
    check_eigen(Afull, E);

    CHECK_MATRIX_NEAR(eigenvalues_sym(M{{2, 1}, {1, 2}}), Ref::from(M{{1}, {3}}), 1e-15);
    CHECK_MATRIX_EQ(eigenvalues_sym(M{{3, 0, 0}, {0, -1, 0}, {0, 0, 2}}), Ref::from(M{{-1}, {2}, {3}}));

    // The identity and the zero matrix: no rotations
    const auto I = eig_sym(eye<double>(4));
    CHECK_MATRIX_EQ(I.values(), Ref::from(M{{1}, {1}, {1}, {1}}));
    CHECK_MATRIX_EQ(I.vectors(), Ref::from(eye<double>(4)));
    CHECK_MATRIX_EQ(eigenvalues_sym(zeros<double>(3, 3)), Ref::from(M{{0}, {0}, {0}}));

    // A small eigenvalue keeps its relative accuracy: det = 1e-20 and trace = 1 + 2e-20, so the
    // eigenvalues are 1e-20 and 1 (to 1e-20); an error of eps ||A|| would be 1e-16
    const M G = {{1, 1e-10}, {1e-10, 2e-20}};
    const auto Eg = eig_sym(G);
    CHECK(std::abs(Eg.values()(1) - 1e-20) <= 1e-14 * 1e-20);
    CHECK(std::abs(Eg.values()(2) - 1.0) <= 1e-15);

    // Empty matrices, and matrices that are not square
    CHECK_EQ(eig_sym(M(0, 0)).size(), index_t(0));
    CHECK_MATL_ERROR(dimension_error, (void)eig_sym(M(2, 3)));
    CHECK_MATL_ERROR(dimension_error, (void)eigenvalues_sym(M(3, 2)));
}

TEST_CASE("Eigen: fixed-size matrices, element types and principal stresses")
{
    // Principal stresses and directions of a stress tensor (MPa)
    const Matrix<double, 3, 3> sigma = {{50, 30, 20}, {30, -20, -10}, {20, -10, 10}};
    const auto P = eig_sym(sigma);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(P)>, FixedSymmetricEigen<double, 3, StorageOrder::RowMajor>>);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(P.values())>, Matrix<double, 3, 1>>);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(P.vectors())>, Matrix<double, 3, 3>>);
    check_eigen(sigma, P);
    CHECK(std::abs(P.values()(1) * P.values()(2) * P.values()(3) - det(sigma)) <= 1e-9);  // the third invariant

    // The same results as the dynamic decomposition, bit for bit
    const auto D = eig_sym(M(sigma));
    CHECK_MATRIX_EQ(P.values(), Ref::from(D.values()));
    CHECK_MATRIX_EQ(P.vectors(), Ref::from(D.vectors()));
    static_assert(std::is_same_v<decltype(eigenvalues_sym(sigma)), Matrix<double, 3, 1>>);
    CHECK_MATRIX_EQ(eigenvalues_sym(sigma), Ref::from(P.values()));

    // Expressions, integer and float matrices
    static_assert(std::is_same_v<decltype(eig_sym(sigma + transpose(sigma)).values()), Matrix<double, 3, 1> const&>);
    const auto Ei = eig_sym(Matrix<int>{{2, 1}, {1, 2}});
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(Ei.values())>, M>);
    CHECK_MATRIX_NEAR(Ei.values(), Ref::from(M{{1}, {3}}), 1e-15);
    const auto Ef = eig_sym(Matrix<float>{{2, 1}, {1, 2}});
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(Ef.values())>, Matrix<float>>);
    CHECK(std::abs(Ef.values()(1) - 1.0f) <= 1e-6f);
}

/****************************************** SVD **************************************************/

TEST_CASE("SVD: tall, wide and square matrices in either storage order")
{
    const index_t shapes[][2] = {{1, 1}, {1, 4}, {4, 1}, {5, 3}, {3, 5}, {12, 12}, {20, 7}, {7, 20}, {40, 17}};
    for (auto const& shape : shapes)
    {
        const M A = random_matrix<M>(shape[0], shape[1]);
        const auto S = svd(A);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(S)>, SVD<double, StorageOrder::RowMajor>>);
        CHECK_EQ(S.rows(), A.rows());
        CHECK_EQ(S.cols(), A.cols());
        check_svd(A, S);
        CHECK_MATRIX_EQ(singular_values(A), Ref::from(S.singular_values()));      // the same rotations

        // sigma^2 are the eigenvalues of A^T A (or A A^T)
        const M G = (A.rows() >= A.cols()) ? M(transpose(A) * A) : M(A * transpose(A));
        const M lambda = eigenvalues_sym(G);
        const index_t k = S.singular_values().rows();
        for (index_t j = 1; j <= k; ++j)
            CHECK(std::abs(S.singular_values()(j) * S.singular_values()(j) - lambda(k + 1 - j)) <= 1e-12 * norm_fro(G));

        const C Ac = A;
        const auto Sc = svd(Ac);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(Sc.U())>, C>);
        check_svd(Ac, Sc);
        CHECK_MATRIX_EQ(Sc.singular_values(), Ref::from(S.singular_values()));
    }

    // Known singular values
    CHECK_MATRIX_NEAR(singular_values(M{{3, 0}, {4, 5}}), Ref::from(M{{3 * std::sqrt(5.0)}, {std::sqrt(5.0)}}), 1e-14);
    const auto S = svd(M{{3, 0}, {4, 5}});
    CHECK(std::abs(S.cond() - 3.0) <= 1e-14);
    CHECK_EQ(S.rank(), index_t(2));
}

TEST_CASE("SVD: rank deficiency, the pseudo-inverse and minimum-norm solutions")
{
    // A rank-3 8 x 6 matrix
    const M A = random_matrix<M>(8, 3) * random_matrix<M>(3, 6);
    const auto S = svd(A);
    check_svd(A, S);
    CHECK_EQ(S.rank(), index_t(3));
    CHECK_EQ(rank(A), index_t(3));
    CHECK_EQ(rank(transpose(A)), index_t(3));
    check_pseudo_inverse(A, pinv(A));
    check_pseudo_inverse(M(transpose(A)), pinv(transpose(A)));
    CHECK_MATRIX_NEAR(pinv(transpose(A)), Ref::from(M(transpose(pinv(A)))), 1e-12);

    // Minimum-norm least squares: X = A+ B, orthogonal to the null space of A (in A's row space)
    const M B = random_matrix<M>(8, 2);
    const M X = S.solve(B);
    CHECK_MATRIX_NEAR(X, Ref::from(M(pinv(A) * B)), 1e-12);
    CHECK(norm_fro(M(transpose(A) * (A * X - B))) <= 1e-10 * norm_fro(A) * norm_fro(B));   // normal equations
    const M null_space = S.V()(all, range(4, last));                                        // V's columns for zero sigma
    CHECK(norm_fro(M(transpose(null_space) * X)) <= 1e-12 * norm_fro(X) + 1e-12);

    // Full rank: the inverse for a square matrix, least squares (as QR) for a tall one, and the
    // smallest-norm solution A^T (A A^T)^-1 b for a wide one
    const M Q = random_matrix<M>(5, 5) + 10.0 * eye<double>(5);
    CHECK_MATRIX_NEAR(pinv(Q), Ref::from(inv(Q)), 1e-14);
    const M T = random_matrix<M>(9, 4), b = random_matrix<M>(9, 1);
    CHECK_MATRIX_NEAR(svd(T).solve(b), Ref::from(qr_solve(T, b)), 1e-12);
    const M W = random_matrix<M>(3, 7), c = random_matrix<M>(3, 1);
    CHECK_MATRIX_NEAR(svd(W).solve(c), Ref::from(M(transpose(W) * lu_solve(W * transpose(W), c))), 1e-12);

    // Exactly zero singular values: U and V stay orthonormal
    M Z = random_matrix<M>(5, 3);
    for (index_t i = 1; i <= 5; ++i) Z(i, 2) = 0.0;
    check_svd(Z, svd(Z));
    CHECK_EQ(svd(Z).singular_values()(3), 0.0);
    CHECK_EQ(svd(Z).cond(), std::numeric_limits<double>::infinity());
    const auto S0 = svd(zeros<double>(4, 3));
    check_svd(zeros<double>(4, 3), S0);
    CHECK_EQ(S0.rank(), index_t(0));
    CHECK_MATRIX_EQ(S0.pinv(), Ref::from(zeros<double>(3, 4)));
    check_svd(zeros<double>(2, 5), svd(zeros<double>(2, 5)));

    // An explicit tolerance truncates small singular values
    const M D = {{1, 0, 0}, {0, 1e-3, 0}, {0, 0, 1e-10}};
    CHECK_EQ(rank(D), index_t(3));
    CHECK_EQ(rank(D, 1e-6), index_t(2));
    CHECK_MATRIX_NEAR(pinv(D, 1e-6), Ref::from(M{{1, 0, 0}, {0, 1e3, 0}, {0, 0, 0}}), 1e-9);
    CHECK_MATRIX_NEAR(svd(D).solve(M{{1}, {1}, {1}}, 1e-6), Ref::from(M{{1}, {1e3}, {0}}), 1e-9);

    // Empty matrices and dimension errors
    CHECK_EQ(svd(M(0, 3)).singular_values().rows(), index_t(0));
    CHECK_EQ(rank(M(3, 0)), index_t(0));
    CHECK_MATL_ERROR(dimension_error, (void)S.solve(M(7, 1)));
}

TEST_CASE("SVD: fixed-size matrices and element types")
{
    Matrix<double, 4, 3> F;
    matl_test::fill_random(F);
    const auto S = svd(F);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(S)>, FixedSVD<double, 4, 3, StorageOrder::RowMajor>>);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(S.U())>, Matrix<double, 4, 3>>);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(S.V())>, Matrix<double, 3, 3>>);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(S.singular_values())>, Matrix<double, 3, 1>>);
    static_assert(std::is_same_v<decltype(S.pinv()), Matrix<double, 3, 4>>);
    static_assert(std::is_same_v<decltype(S.solve(Matrix<double, 4, 2>{})), Matrix<double, 3, 2>>, "fixed in, fixed out");
    static_assert(std::is_same_v<decltype(S.solve(M(4, 2))), M>);
    check_svd(F, S);

    // The same results as the dynamic decomposition, bit for bit
    const auto D = svd(M(F));
    CHECK_MATRIX_EQ(S.singular_values(), Ref::from(D.singular_values()));
    CHECK_MATRIX_EQ(S.U(), Ref::from(D.U()));
    CHECK_MATRIX_EQ(S.V(), Ref::from(D.V()));
    CHECK_MATRIX_EQ(S.pinv(), Ref::from(D.pinv()));

    // Wide fixed-size matrices
    const Matrix<double, 2, 5> W = {{1, 2, 3, 4, 5}, {2, 0, 1, 0, 2}};
    const auto Sw = svd(W);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(Sw.U())>, Matrix<double, 2, 2>>);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(Sw.V())>, Matrix<double, 5, 2>>);
    check_svd(W, Sw);
    static_assert(std::is_same_v<decltype(pinv(W)), Matrix<double, 5, 2>>);
    check_pseudo_inverse(W, pinv(W));
    static_assert(std::is_same_v<decltype(singular_values(W)), Matrix<double, 2, 1>>);
    CHECK_MATRIX_EQ(singular_values(W), Ref::from(Sw.singular_values()));

    // Integer and float matrices
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(svd(Matrix<int>{{1}}).U())>, M>);
    const auto Sf = svd(Matrix<float>{{3, 0}, {4, 5}});
    CHECK(std::abs(Sf.singular_values()(2) - std::sqrt(5.0f)) <= 1e-5f);
}
