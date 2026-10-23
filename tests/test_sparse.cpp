/**************************************************************************************************
 * File: test_sparse.cpp
 *
 * Description: Sparse matrices (CSRMatrix): assembly with SparseBuilder, construction from CSR
 *              arrays and its checks, element access and changing values, conversions to and from
 *              dense matrices, products with dense vectors and matrices (also transposed, nested,
 *              and with targets the product reads), scaling, output, and errors. Results are
 *              compared with the same operations on dense(S).
 **************************************************************************************************/
#include <cstdint>
#include <sstream>
#include <type_traits>
#include <vector>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using M = RowMajorMatrix<double>;
    using C = ColMajorMatrix<double>;
    using S = CSRMatrix<double>;

    template<typename MAT>
    MAT random_matrix(index_t m, index_t n)
    {
        MAT A(m, n);
        matl_test::fill_random(A);
        return A;
    }

    // A random sparse matrix with about 'fill' of its elements stored
    template<class I = index_t>
    CSRMatrix<double, I> random_sparse(index_t m, index_t n, double fill = 0.2)
    {
        M D = random_matrix<M>(m, n);
        std::uniform_real_distribution<double> u(0.0, 1.0);
        for (index_t k = 1; k <= D.size(); ++k)
            if (u(matl_test::rng()) > fill) D(k) = 0.0;
        return sparse<I>(D);
    }

    // The 2D Laplacian (5-point stencil) on a k x k grid: the standard test matrix of PDE solvers
    S laplacian_2d(index_t k)
    {
        const index_t n = k * k;
        SparseBuilder<double> B(n, n);
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

    // The elements of a span, as a vector (to compare). Built element by element: constructing a
    // vector from span iterators trips a false -Wfree-nonheap-object in GCC 15 at -O3.
    template<typename T>
    std::vector<std::remove_const_t<T>> as_vector(std::span<T> s)
    {
        std::vector<std::remove_const_t<T>> v;
        v.reserve(s.size());
        for (auto x : s) v.push_back(x);
        return v;
    }

    template<typename F>
    std::size_t temporary_bytes(F f)
    {
        MemoryManager::temporary_arena().reset_high_water_mark();
        f();
        return MemoryManager::temporary_arena().high_water_mark();
    }
}

/********************************************* Building ******************************************/

TEST_CASE("Sparse: building from entries in any order, with duplicates summed")
{
    SparseBuilder<double> B(3, 4);
    B.add(3, 4, 4.0);
    B.add(1, 3, 1.5);
    B.add(1, 1, 1.0);
    B.add(3, 2, 3.0);
    B.add(1, 3, 0.5);                                   // (1,3) again: 1.5 + 0.5
    B.add(2, 2, 0.0);                                   // an explicit zero is stored
    CHECK_EQ(B.entries(), 6ul);
    S A = B.build();

    CHECK_EQ(A.rows(), 3ul);
    CHECK_EQ(A.cols(), 4ul);
    CHECK_EQ(A.nnz(), 5ul);
    CHECK(as_vector(A.row_offsets()) == std::vector<index_t>({0, 2, 3, 5}));
    CHECK(as_vector(A.column_indices()) == std::vector<index_t>({0, 2, 1, 1, 3}));
    CHECK(as_vector(A.values()) == std::vector<double>({1.0, 2.0, 0.0, 3.0, 4.0}));
    CHECK_MATRIX_EQ(dense(A), Ref::from(M{{1, 0, 2, 0}, {0, 0, 0, 0}, {0, 3, 0, 4}}));

    // The builder is unchanged, so it can build again, or be cleared and reused
    CHECK_EQ(B.build().nnz(), 5ul);
    B.clear();
    CHECK_EQ(B.entries(), 0ul);
    CHECK_EQ(B.build().nnz(), 0ul);

    // Empty matrices, empty rows and columns
    CHECK_EQ(SparseBuilder<double>(0, 0).build().nnz(), 0ul);
    CHECK_EQ(SparseBuilder<double>(5, 2).build().row_offsets().size(), 6ul);

    // A 32-bit index type
    SparseBuilder<float, std::int32_t> B32(2, 2);
    B32.add(2, 1, 5.0f);
    CSRMatrix<float, std::int32_t> A32 = B32.build();
    static_assert(std::is_same_v<decltype(A32)::index_type, std::int32_t>);
    CHECK_EQ(A32(2, 1), 5.0f);

    // Indices outside the matrix are always reported
    CHECK_MATL_ERROR(index_error, B.add(4, 1, 1.0));
    CHECK_MATL_ERROR(index_error, B.add(1, 0, 1.0));
}

TEST_CASE("Sparse: construction from CSR arrays is checked")
{
    S A(2, 3, {0, 1, 3}, {2, 0, 1}, {7.0, 8.0, 9.0});
    CHECK_MATRIX_EQ(dense(A), Ref::from(M{{0, 0, 7}, {8, 9, 0}}));
    CHECK_EQ(S(4, 5).nnz(), 0ul);                       // no stored entries: all zero
    CHECK_EQ(S(4, 5)(4, 5), 0.0);

    CHECK_MATL_ERROR(dimension_error, (S(2, 3, {0, 1}, {2}, {7.0})));                    // too few offsets
    CHECK_MATL_ERROR(dimension_error, (S(2, 3, {0, 1, 2}, {2, 0}, {7.0})));              // sizes differ
    CHECK_MATL_ERROR(index_error, (S(2, 3, {0, 1, 3}, {2, 3, 1}, {7.0, 8.0, 9.0})));     // column 3 (0-based)
    CHECK_MATL_ERROR(index_error, (S(2, 3, {0, 2, 2}, {1, 0}, {7.0, 8.0})));             // not increasing
    CHECK_MATL_ERROR(index_error, (S(2, 3, {0, 2, 1}, {0, 1}, {7.0, 8.0})));             // offsets decrease
    CHECK_MATL_ERROR(index_error, (CSRMatrix<double, int>(1, 3, {0, 1}, {-1}, {7.0})));   // negative
}

/********************************************* Elements ******************************************/

TEST_CASE("Sparse: element access and changing values")
{
    S A(2, 3, {0, 1, 3}, {2, 0, 1}, {7.0, 8.0, 9.0});
    CHECK_EQ(A(1, 3), 7.0);
    CHECK_EQ(A(1, 1), 0.0);                             // not stored
    CHECK(A.is_stored(2, 2));
    CHECK(!A.is_stored(2, 3));

    A.entry(2, 1) += 1.0;                               // e.g. re-assembly into a fixed pattern
    CHECK_EQ(A(2, 1), 9.0);
    for (double& v : A.values()) v *= 2.0;
    CHECK_EQ(A(2, 2), 18.0);

    CHECK_MATL_ERROR(index_error, A.entry(1, 1) = 1.0);  // the pattern cannot change
    CHECK_MATL_ERROR(index_error, (void)A.is_stored(3, 1));
    #ifdef MATL_BOUNDS_CHECK
    CHECK_MATL_ERROR(index_error, (void)A(1, 4));
    #endif
}

/******************************************** Conversions ****************************************/

TEST_CASE("Sparse: conversions to and from dense matrices")
{
    M D = {{0, 1.5, 0}, {0, 0, 0}, {-2, 0, 3}};
    S A = sparse(D);
    CHECK_EQ(A.nnz(), 3ul);                             // zeros are not stored
    CHECK_MATRIX_EQ(dense(A), Ref::from(D));
    static_assert(std::is_same_v<decltype(dense(A)), M>);

    CHECK_MATRIX_EQ(dense(sparse(D(range(1, 2), all))), Ref::from(M(D(range(1, 2), all))));  // a view
    CHECK_MATRIX_EQ(dense(sparse(transpose(D) * 2.0)), Ref::from(M(transpose(D) * 2.0)));     // an expression
    CHECK_MATRIX_EQ(dense(sparse(C(D))), Ref::from(D));

    auto I = sparse(Matrix<int>{{0, 4}, {5, 0}});
    static_assert(std::is_same_v<decltype(I), CSRMatrix<int>>);
    CHECK_EQ(I(2, 1), 5);
    auto F = sparse<std::int32_t>(Matrix<float>{{1, 0}});
    static_assert(std::is_same_v<decltype(F), CSRMatrix<float, std::int32_t>>);
}

/********************************************* Products ******************************************/

TEST_CASE("Sparse: products with dense vectors and matrices")
{
    for (index_t n : {1, 7, 60})
    {
        S A = random_sparse(n + 3, n);
        const Ref a = Ref::from(dense(A));
        M x = random_matrix<M>(n, 1), X = random_matrix<M>(n, 4);
        M z = random_matrix<M>(n + 3, 1), Z = random_matrix<M>(n + 3, 3);

        CHECK_MATRIX_EQ(M(A * x), a * Ref::from(x));
        CHECK_MATRIX_EQ(M(A * X), a * Ref::from(X));
        CHECK_MATRIX_EQ(M(transpose(A) * z), transpose(a) * Ref::from(z));
        CHECK_MATRIX_EQ(M(transpose(A) * Z), transpose(a) * Ref::from(Z));

        // Column-major operands and targets, views, and expressions as operands
        C Xc = X, Y(n + 3, 4);
        Y = A * Xc;                                         CHECK_MATRIX_EQ(Y, a * Ref::from(X));
        CHECK_MATRIX_EQ(M(A * X(all, range(2, 3))), a * Ref::from(M(X(all, range(2, 3)))));
        CHECK_MATRIX_EQ(M(A * (x + 1.0)), a * (Ref::from(x) + 1.0));
        M W(n + 5, 4);
        W = 0.0;
        W(range(2, n + 4), all) = A * X;                    // into a view
        CHECK_MATRIX_EQ(W(range(2, n + 4), all), a * Ref::from(X));
        CHECK_EQ(W(1, 1), 0.0);

        // Nested in larger expressions, and 32-bit indices
        CHECK_MATRIX_EQ(M(z - A * x), Ref::from(z) - a * Ref::from(x));
        CHECK_MATRIX_EQ(M(2.0 * (A * x) + transpose(transpose(A * x))), 3.0 * (a * Ref::from(x)));
        CSRMatrix<double, std::int32_t> A32 = sparse<std::int32_t>(dense(A));
        CHECK_MATRIX_EQ(M(A32 * X), a * Ref::from(X));
        CHECK_MATRIX_EQ(M(transpose(A32) * z), transpose(a) * Ref::from(z));
    }

    // The 2D Laplacian times a constant vector: 0 inside the grid, positive at the boundary
    S L = laplacian_2d(5);
    M y = L * ones<double>(25, 1);
    CHECK_EQ(y(13), 0.0);                                   // the center of the grid
    CHECK_EQ(y(1), 2.0);                                    // a corner (two neighbors)
    CHECK_EQ(y(3), 1.0);                                    // an edge (three neighbors)
    CHECK_EQ(L.nnz(), 25ul + 4ul * 20ul);
}

TEST_CASE("Sparse: products whose target is read by the product")
{
    S A = random_sparse(30, 30, 0.3);
    const Ref a = Ref::from(dense(A));
    M x = random_matrix<M>(30, 2), y(30, 2);
    const Ref x0 = Ref::from(x);

    CHECK_EQ(temporary_bytes([&] {y = A * x;}), 0ul);      // x is another matrix: no temporary
    CHECK_MATRIX_EQ(y, a * x0);
    CHECK(temporary_bytes([&] {x = A * x;}) > 0ul);         // x is read and written: a temporary
    CHECK_MATRIX_EQ(x, a * x0);
    x = A * x;
    x = transpose(A) * x;
    CHECK_MATRIX_EQ(x, transpose(a) * (a * (a * x0)));
    x(all, 1) = A * x(all, 2);                              // a view of the operand
    CHECK_MATRIX_NEAR(x(all, 1), a * Ref::from(M(x(all, 2))), 0.0);

    // A stored expression refers to the sparse matrix, and reads its current values
    M v = ones<double>(30, 1);
    auto e = A * v;
    for (double& value : A.values()) value = 1.0;
    M row_counts = e;
    CHECK_EQ(row_counts(1), double(A.row_offsets()[1] - A.row_offsets()[0]));
}

/***************************************** Scaling and output ************************************/

TEST_CASE("Sparse: scaling and output")
{
    S A(2, 3, {0, 1, 3}, {2, 0, 1}, {7.0, 8.0, 9.0});
    CHECK_MATRIX_EQ(dense(A * 2.0), Ref::from(M{{0, 0, 14}, {16, 18, 0}}));
    CHECK_MATRIX_EQ(dense(0.5 * A), Ref::from(M{{0, 0, 3.5}, {4, 4.5, 0}}));
    CHECK_MATRIX_EQ(dense(A / 2), Ref::from(M{{0, 0, 3.5}, {4, 4.5, 0}}));
    CHECK_MATRIX_EQ(dense(-A), Ref::from(M{{0, 0, -7}, {-8, -9, 0}}));
    static_assert(std::is_same_v<decltype(CSRMatrix<int>() * 2.5), CSRMatrix<double>>);
    A *= 3;
    A /= 3.0;
    CHECK_EQ(A(2, 2), 9.0);

    std::ostringstream os;
    os << A;
    CHECK_EQ(os.str(), std::string("2x3 sparse matrix, 3 stored entries\n(1,3)  7\n(2,1)  8\n(2,2)  9"));
}

TEST_CASE("Sparse: dimension errors")
{
    S A = random_sparse(4, 3);
    M y(4, 1);
    // Always checked, also without debug checks
    CHECK_MATL_ERROR(dimension_error, y = A * zeros<double>(4, 1));
    CHECK_MATL_ERROR(dimension_error, y = transpose(A) * zeros<double>(3, 1));
    CHECK_MATL_ERROR(dimension_error, (void)sum(A * zeros<double>(2, 1)));
    #ifdef MATL_DIMENSION_CHECK
    CHECK_MATL_ERROR(dimension_error, y = A * zeros<double>(3, 2));     // the target's dimensions
    #endif
}
