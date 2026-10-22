/**************************************************************************************************
 * File: test_views.cpp
 *
 * Description: Views (MatrixView): selecting rows and columns with range, all and last, reading
 *              and writing through views, views in expressions and products, views of views,
 *              column-major and fixed-size matrices, assignments whose target overlaps what they
 *              read, views of temporaries, and errors.
 **************************************************************************************************/
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

    // Reference: the given rows and columns of r (1-based, in the given order)
    Ref pick(Ref const& r, std::vector<index_t> const& rows, std::vector<index_t> const& cols)
    {
        Ref s(rows.size(), cols.size());
        for (index_t i = 1; i <= rows.size(); ++i)
            for (index_t j = 1; j <= cols.size(); ++j)
                s(i, j) = r(rows[i - 1], cols[j - 1]);
        return s;
    }

    // Reference: r with the given rows and columns replaced by the elements of v
    Ref put(Ref r, std::vector<index_t> const& rows, std::vector<index_t> const& cols, Ref const& v)
    {
        for (index_t i = 1; i <= rows.size(); ++i)
            for (index_t j = 1; j <= cols.size(); ++j)
                r(rows[i - 1], cols[j - 1]) = v(i, j);
        return r;
    }

    // A matrix with elements 1, 2, 3, ... (row by row)
    template<typename MAT>
    MAT sequence(index_t m, index_t n)
    {
        MAT A(m, n);
        for (index_t i = 1; i <= m; ++i)
            for (index_t j = 1; j <= n; ++j)
                A(i, j) = double((i - 1) * n + j);
        return A;
    }

    template<typename MAT>
    MAT random_matrix(index_t m, index_t n)
    {
        MAT A(m, n);
        matl_test::fill_random(A);
        return A;
    }

    // Bytes of evaluation temporaries used by f
    template<typename F>
    std::size_t temporary_bytes(F f)
    {
        MemoryManager::temporary_arena().reset_high_water_mark();
        f();
        return MemoryManager::temporary_arena().high_water_mark();
    }
}

/********************************************* Types *********************************************/

TEST_CASE("Views: types")
{
    M A(3, 4);
    M const& K = A;
    static_assert(std::is_same_v<decltype(A(all, 1)), MatrixView<double, StorageOrder::RowMajor>>);
    static_assert(std::is_same_v<decltype(K(all, 1)), MatrixView<const double, StorageOrder::RowMajor>>,
                  "views of const matrices are read-only");
    static_assert(std::is_same_v<decltype(M(3, 4)(all, 1)), MatrixView<const double, StorageOrder::RowMajor>>,
                  "views of temporaries are read-only");
    static_assert(std::is_same_v<decltype(C(3, 4)(all, 1)), MatrixView<const double, StorageOrder::ColMajor>>);
    static_assert(std::is_same_v<decltype(A(1, 2)), double&>, "two integers select an element");
    static_assert(std::is_same_v<decltype(A(last, 2)), double&>, "so do integers and last");

    using V = MatrixView<double, StorageOrder::RowMajor>;
    static_assert(ExpressionType<V> && !MatrixType<V> && WritableMatrix<V>);
    static_assert(!WritableMatrix<MatrixView<const double, StorageOrder::RowMajor>>);
    static_assert(std::is_same_v<detail::dtype_t<MatrixView<const float, StorageOrder::ColMajor>>, float>);
    static_assert(std::is_same_v<detail::result_matrix_t<V>, M>);
    static_assert(std::is_same_v<detail::result_matrix_t<MatrixView<double, StorageOrder::ColMajor>>, C>);
    CHECK(true);
}

/************************************* Selecting elements ****************************************/

TEST_CASE("Views: rows, columns and blocks")
{
    M A = sequence<M>(4, 5);
    Ref a = Ref::from(A);

    CHECK_MATRIX_EQ(A(range(1, 3), all), pick(a, {1, 2, 3}, {1, 2, 3, 4, 5}));
    CHECK_MATRIX_EQ(A(all, 2),           pick(a, {1, 2, 3, 4}, {2}));
    CHECK_MATRIX_EQ(A(3, all),           pick(a, {3}, {1, 2, 3, 4, 5}));
    CHECK_MATRIX_EQ(A(range(2, 3), range(2, 4)), pick(a, {2, 3}, {2, 3, 4}));
    CHECK_MATRIX_EQ(A(all, all),         a);

    auto v = A(range(2, 3), range(2, 4));
    CHECK_EQ(v.rows(), 2ul);
    CHECK_EQ(v.cols(), 3ul);
    CHECK_EQ(v.size(), 6ul);
}

TEST_CASE("Views: steps, backwards ranges and last")
{
    M A = sequence<M>(5, 6);
    Ref a = Ref::from(A);

    CHECK_MATRIX_EQ(A(range(1, 2, last), all),        pick(a, {1, 3, 5}, {1, 2, 3, 4, 5, 6}));
    CHECK_MATRIX_EQ(A(all, range(2, 2, last)),        pick(a, {1, 2, 3, 4, 5}, {2, 4, 6}));
    CHECK_MATRIX_EQ(A(range(1, 2, 6), all),           pick(a, {1, 3, 5}, {1, 2, 3, 4, 5, 6}));  // up to 6
    CHECK_MATRIX_EQ(A(range(last, -1, 1), all),       pick(a, {5, 4, 3, 2, 1}, {1, 2, 3, 4, 5, 6}));
    CHECK_MATRIX_EQ(A(range(4, -2, 1), range(last, -3, 1)), pick(a, {4, 2}, {6, 3}));
    CHECK_MATRIX_EQ(A(last, all),                     pick(a, {5}, {1, 2, 3, 4, 5, 6}));
    CHECK_MATRIX_EQ(A(range(2, last - 1), last - 2),  pick(a, {2, 3, 4}, {4}));
    CHECK_MATRIX_EQ(A(range(last - 1, last), range(1, 2)), pick(a, {4, 5}, {1, 2}));
}

TEST_CASE("Views: empty ranges")
{
    M A = sequence<M>(3, 4);
    auto e = A(range(3, 2), all);                   // the step points away from the last index
    CHECK_EQ(e.rows(), 0ul);
    CHECK_EQ(e.cols(), 4ul);
    CHECK_EQ(A(range(last, -1, last), all).rows(), 1ul);

    // An empty range may start just past the end (MATLAB A(k+1:end, :) with k = end)
    auto tail = A(range(4, last), all);
    CHECK_EQ(tail.rows(), 0ul);
    CHECK_EQ(sum(tail), 0.0);
    tail = 5.0;                                     // writes nothing
    tail = tail * 2.0;
    CHECK_MATRIX_EQ(A, Ref::from(sequence<M>(3, 4)));

    M copy = A(range(1, 0), all);
    CHECK_EQ(copy.size(), 0ul);
}

TEST_CASE("Views: element access")
{
    M A = sequence<M>(3, 4);
    auto v = A(range(2, 3), range(2, 4));
    CHECK_EQ(v(1, 1), 6.0);
    CHECK_EQ(v(2, 3), 12.0);
    CHECK_EQ(v(last, last), 12.0);
    CHECK_EQ(v(1, last - 1), 7.0);
    CHECK_EQ(v(4), 10.0);                           // storage order of A (row by row)

    v(1, 2) = -1.0;                                 // writes A(2,3)
    CHECK_EQ(A(2, 3), -1.0);

    auto column = A(all, 4);
    CHECK_EQ(column(3), 12.0);                      // vectors: element k
    auto row = A(2, all);
    CHECK_EQ(row(3), -1.0);

    // last in element access
    CHECK_EQ(A(last, last), 12.0);
    CHECK_EQ(A(last, 1), 9.0);
    CHECK_EQ(A(1, last - 1), 3.0);
    A(last, last) = 0.0;
    CHECK_EQ(A(3, 4), 0.0);
    M const& K = A;
    CHECK_EQ(K(last - 2, last), 4.0);

    Matrix<double, 2, 3> F = {{1, 2, 3}, {4, 5, 6}};
    CHECK_EQ(F(last, last - 1), 5.0);
    CHECK_EQ(F.data()[1], 2.0);
    CHECK_EQ(A.data()[0], 1.0);
}

/************************************ Writing through views **************************************/

TEST_CASE("Views: assignment writes the viewed elements")
{
    M A = sequence<M>(4, 5);
    Ref a = Ref::from(A);
    const std::vector<index_t> rows = {2, 3}, cols = {1, 3, 5};

    A(all, 1) = 0.0;
    a = put(a, {1, 2, 3, 4}, {1}, Ref(4, 1, 0.0));
    CHECK_MATRIX_EQ(A, a);

    M B = random_matrix<M>(2, 3);
    A(range(2, 3), range(1, 2, 5)) = B;                             // a matrix
    a = put(a, rows, cols, Ref::from(B));
    CHECK_MATRIX_EQ(A, a);

    A(range(2, 3), range(1, 2, 5)) = 2.0 * B + 1.0;                 // an expression
    a = put(a, rows, cols, 2.0 * Ref::from(B) + 1.0);
    CHECK_MATRIX_EQ(A, a);

    M D = random_matrix<M>(5, 3);
    A(range(2, 3), range(1, 2, 5)) = D(range(1, 2), all);           // a view of another matrix
    a = put(a, rows, cols, pick(Ref::from(D), {1, 2}, {1, 2, 3}));
    CHECK_MATRIX_EQ(A, a);

    A(1, all) = transpose(D(all, 2));                               // a transposed column
    a = put(a, {1}, {1, 2, 3, 4, 5}, transpose(pick(Ref::from(D), {1, 2, 3, 4, 5}, {2})));
    CHECK_MATRIX_EQ(A, a);
}

TEST_CASE("Views: arithmetic assignment")
{
    M A = sequence<M>(3, 4);
    Ref a = Ref::from(A);
    auto v = A(range(1, 2), range(2, 3));
    Ref w = pick(a, {1, 2}, {2, 3});
    M B = random_matrix<M>(2, 2);
    Ref b = Ref::from(B);

    v += 1.0;    w = w + 1.0;    CHECK_MATRIX_EQ(v, w);
    v -= 3;      w = w - 3.0;    CHECK_MATRIX_EQ(v, w);
    v *= 2.0;    w = w * 2.0;    CHECK_MATRIX_EQ(v, w);
    v /= 4.0;    w = w / 4.0;    CHECK_MATRIX_EQ(v, w);
    v += B;      w = w + b;      CHECK_MATRIX_EQ(v, w);
    v -= 2.0 * B;  w = w - 2.0 * b;  CHECK_MATRIX_EQ(v, w);
    v *= B;      w = w * b;      CHECK_MATRIX_EQ(v, w);            // matrix product
    CHECK_MATRIX_EQ(A, put(a, {1, 2}, {2, 3}, w));

    Matrix<int> I = {{7, 8}, {9, 10}};
    I(all, 1) /= 2;                                                 // integer division
    CHECK_EQ(I(1, 1), 3);
    CHECK_EQ(I(2, 1), 4);
    CHECK_EQ(I(2, 2), 10);
}

TEST_CASE("Views: copying a view copies the view, not the elements")
{
    M A = sequence<M>(3, 3);
    auto v = A(all, 1);
    auto w = v;                                     // another view of column 1
    w = 0.0;
    CHECK_EQ(A(2, 1), 0.0);
    CHECK_EQ(v(3), 0.0);

    M c = A(all, 2);                                // a new matrix: a copy of column 2
    c = 9.0;
    CHECK_EQ(A(1, 2), 2.0);

    v = A(all, 3);                                  // assignment writes column 1 of A
    CHECK_EQ(A(2, 1), 6.0);
    CHECK_EQ(v(2), 6.0);
}

/*********************************** Views in expressions ****************************************/

TEST_CASE("Views: operands of elementwise operations, functions and reductions")
{
    M A = random_matrix<M>(5, 6), B = random_matrix<M>(3, 3), R(3, 3);
    Ref a = Ref::from(A), b = Ref::from(B);
    auto v = A(range(1, 2, 5), range(2, 4));
    Ref w = pick(a, {1, 3, 5}, {2, 3, 4});

    R = v + B;                      CHECK_MATRIX_EQ(R, w + b);
    R = 2.0 * v - B * 3.0;          CHECK_MATRIX_EQ(R, 2.0 * w - b * 3.0);
    R = transpose(v);               CHECK_MATRIX_EQ(R, transpose(w));
    R = v + transpose(v);           CHECK_MATRIX_EQ(R, w + transpose(w));
    R = abs(v) + times(v, B);       CHECK_MATRIX_EQ(R, Ref::from(M(abs(v))) + Ref::from(M(times(v, B))));

    double total = 0, largest = w.data[0];
    for (double x : w.data) {total += x; largest = std::max(largest, x);}
    CHECK_EQ(sum(v), total);
    CHECK_EQ(max(v), largest);
    CHECK_EQ(sum(v + B), total + sum(B));
    CHECK_EQ(trace(v), w(1, 1) + w(2, 2) + w(3, 3));
    CHECK_EQ(dot(A(all, 1), A(all, 2)), dot(M(A(all, 1)), M(A(all, 2))));
}

TEST_CASE("Views: products, small and blocked")
{
    for (index_t n : {3, 70})
    {
        M A = random_matrix<M>(n + 2, n + 3), B = random_matrix<M>(n + 1, n);
        Ref a = Ref::from(A), b = Ref::from(B);
        std::vector<index_t> rows, cols;
        for (index_t k = 2; k <= n + 1; ++k) rows.push_back(k);
        for (index_t k = 1; k <= n + 1; ++k) cols.push_back(k);

        auto v = A(range(2, n + 1), range(1, n + 1));               // n x (n+1) view
        Ref w = pick(a, rows, cols);

        M P(n, n);
        P = v * B;                       CHECK_MATRIX_EQ(P, w * b);
        P = transpose(B) * transpose(v); CHECK_MATRIX_EQ(P, transpose(b) * transpose(w));
        P = (v + 1.0) * B;               CHECK_MATRIX_EQ(P, (w + 1.0) * b);

        // A product written into a view
        M T = random_matrix<M>(n + 1, n + 2);
        Ref t = Ref::from(T);
        std::vector<index_t> trows, tcols;
        for (index_t k = 1; k <= n; ++k) trows.push_back(k + 1);
        for (index_t k = 1; k <= n; ++k) tcols.push_back(k + 2);
        T(range(2, n + 1), range(3, n + 2)) = v * B;
        CHECK_MATRIX_EQ(T, put(t, trows, tcols, w * b));
    }
}

TEST_CASE("Views: views of views")
{
    M A = sequence<M>(6, 6);
    Ref a = Ref::from(A);
    auto v = A(range(2, 6), range(1, 2, 6));                       // rows 2-6, columns 1, 3, 5
    auto w = v(range(1, 2, last), range(2, 3));                    // rows 2, 4, 6; columns 3, 5
    CHECK_MATRIX_EQ(w, pick(a, {2, 4, 6}, {3, 5}));
    CHECK_MATRIX_EQ(v(last, all), pick(a, {6}, {1, 3, 5}));

    w = 0.0;
    CHECK_EQ(A(4, 5), 0.0);
    CHECK_EQ(A(4, 4), 22.0);
}

/****************************** Storage orders and fixed sizes ***********************************/

TEST_CASE("Views: column-major matrices")
{
    C A = sequence<C>(4, 5);
    Ref a = Ref::from(A);
    CHECK_MATRIX_EQ(A(range(2, 3), range(2, 2, last)), pick(a, {2, 3}, {2, 4}));
    CHECK_EQ(A(all, 2)(3), 12.0);
    CHECK_EQ(A(range(1, 2), range(1, 2))(2), 6.0);  // storage order of A (column by column)

    M R = random_matrix<M>(2, 2);
    A(range(1, 2), range(4, 5)) = R;                               // row-major into column-major
    CHECK_MATRIX_EQ(A, put(a, {1, 2}, {4, 5}, Ref::from(R)));

    M P(2, 2);
    P = A(range(1, 2), all) * transpose(A(range(3, 4), all));      // mixed storage orders
    Ref p_a = put(a, {1, 2}, {4, 5}, Ref::from(R));
    CHECK_MATRIX_EQ(P, pick(p_a, {1, 2}, {1, 2, 3, 4, 5}) * transpose(pick(p_a, {3, 4}, {1, 2, 3, 4, 5})));

    // External memory in Fortran order: a view of a block of a Fortran array
    double fortran[6] = {1, 4, 2, 5, 3, 6};                        // 2x3 column by column
    ColMajorMatrix<double> F(2, 3, fortran);
    F(all, range(2, 3)) = 0.0;
    CHECK_EQ(fortran[1], 4.0);
    CHECK_EQ(fortran[2], 0.0);
    CHECK_EQ(fortran[5], 0.0);
}

TEST_CASE("Views: fixed-size matrices")
{
    Matrix<double, 3, 3> F = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Ref f = Ref::from(F);
    CHECK_MATRIX_EQ(F(range(2, 3), all), pick(f, {2, 3}, {1, 2, 3}));

    F(all, 1) = F(all, 3);
    CHECK_EQ(F(3, 1), 9.0);

    Matrix<double, 3, 1> c = F(all, 2);                             // dimensions checked at run time
    CHECK_EQ(c(2), 5.0);
    FixedSizeColMajorMatrix<double, 2, 2> G = {{1, 2}, {3, 4}};
    CHECK_EQ(sum(G(last, all)), 7.0);
    CHECK_EQ(G(all, last)(1), 2.0);
}

/******************************* Assignments that overlap ****************************************/

TEST_CASE("Views: assignments whose target overlaps what they read")
{
    M A = sequence<M>(5, 4);
    Ref a = Ref::from(A);

    // Shift rows down and up (each row is read after the row above was written)
    M S = A;
    S(range(2, last), all) = S(range(1, last - 1), all);
    CHECK_MATRIX_EQ(S, put(a, {2, 3, 4, 5}, {1, 2, 3, 4}, pick(a, {1, 2, 3, 4}, {1, 2, 3, 4})));
    S = A;
    S(range(1, last - 1), all) = S(range(2, last), all) + 0.0;
    CHECK_MATRIX_EQ(S, put(a, {1, 2, 3, 4}, {1, 2, 3, 4}, pick(a, {2, 3, 4, 5}, {1, 2, 3, 4})));

    // Reverse the rows in place (target and operand are the same elements in other positions)
    S = A;
    S = S(range(last, -1, 1), all);
    CHECK_MATRIX_EQ(S, pick(a, {5, 4, 3, 2, 1}, {1, 2, 3, 4}));
    S = A;
    S(range(last, -1, 1), all) = S;
    CHECK_MATRIX_EQ(S, pick(a, {5, 4, 3, 2, 1}, {1, 2, 3, 4}));
    S = A;
    S(range(last, -1, 1), all) = S + 1.0;
    CHECK_MATRIX_EQ(S, pick(a, {5, 4, 3, 2, 1}, {1, 2, 3, 4}) + 1.0);

    // Transpose a square block in place, and transpose one block into an overlapping one
    S = A;
    S(range(1, 3), range(1, 3)) = transpose(S(range(1, 3), range(1, 3)));
    CHECK_MATRIX_EQ(S, put(a, {1, 2, 3}, {1, 2, 3}, transpose(pick(a, {1, 2, 3}, {1, 2, 3}))));
    S = A;
    S(range(1, 3), range(1, 3)) = transpose(S(range(2, 4), range(2, 4)));
    CHECK_MATRIX_EQ(S, put(a, {1, 2, 3}, {1, 2, 3}, transpose(pick(a, {2, 3, 4}, {2, 3, 4}))));
    S = A;
    S(range(1, 3), range(1, 3)) = S(range(2, 4), range(2, 4)) + transpose(S(range(1, 3), range(1, 3)));
    CHECK_MATRIX_EQ(S, put(a, {1, 2, 3}, {1, 2, 3},
                           pick(a, {2, 3, 4}, {2, 3, 4}) + transpose(pick(a, {1, 2, 3}, {1, 2, 3}))));

    // A product that reads the elements it writes
    M B = random_matrix<M>(4, 4);
    S = A;
    S(range(1, 2), all) = S(range(1, 2), all) * B;
    CHECK_MATRIX_EQ(S, put(a, {1, 2}, {1, 2, 3, 4}, pick(a, {1, 2}, {1, 2, 3, 4}) * Ref::from(B)));
    S = A;
    S(all, range(1, 2)) = S * B(all, range(3, 4));
    CHECK_MATRIX_EQ(S, put(a, {1, 2, 3, 4, 5}, {1, 2}, a * pick(Ref::from(B), {1, 2, 3, 4}, {3, 4})));

    // Columns of a row-major matrix interleave in memory: copying one into another is correct
    // (it uses a temporary, because the address ranges overlap)
    S = A;
    S(all, 1) = S(all, 2);
    CHECK_MATRIX_EQ(S, put(a, {1, 2, 3, 4, 5}, {1}, pick(a, {1, 2, 3, 4, 5}, {2})));
}

TEST_CASE("Views: no temporary when the target is read only at the positions it writes")
{
    M A = random_matrix<M>(6, 6), B = random_matrix<M>(3, 6);
    Ref a = Ref::from(A), b = Ref::from(B);
    auto v = A(range(1, 2, last), all);
    Ref w = pick(a, {1, 3, 5}, {1, 2, 3, 4, 5, 6});

    CHECK_EQ(temporary_bytes([&] {v = 2.0 * v + B;}), 0ul);
    CHECK_MATRIX_EQ(v, 2.0 * w + b);
    CHECK_EQ(temporary_bytes([&] {v = abs(v) - v;}), 0ul);
    CHECK_EQ(temporary_bytes([&] {A(range(1, 3), all) = B;}), 0ul);         // B is another matrix

    // Overlap at other positions uses a temporary
    CHECK(temporary_bytes([&] {A(range(2, 4), all) = A(range(1, 3), all);}) > 0ul);
}

/********************************************* swap ********************************************/

TEST_CASE("Views: swap exchanges the elements of two views")
{
    M A = sequence<M>(4, 3);
    Ref a = Ref::from(A);

    swap(A(1, all), A(3, all));                                     // rows (temporary views)
    CHECK_MATRIX_EQ(A, pick(a, {3, 2, 1, 4}, {1, 2, 3}));

    A = sequence<M>(4, 3);
    auto v = A(all, 1), w = A(all, 3);
    using std::swap;
    swap(v, w);                                                     // finds matl::swap, not std::swap
    CHECK_MATRIX_EQ(A, pick(a, {1, 2, 3, 4}, {3, 2, 1}));
    CHECK_MATRIX_EQ(v, pick(a, {1, 2, 3, 4}, {3}));                 // the views still view columns 1 and 3

    // Different matrices with different storage orders
    C B = random_matrix<C>(4, 5);
    Ref b = Ref::from(B);
    A = sequence<M>(4, 3);
    matl::swap(A(all, 2), B(all, 5));
    CHECK_MATRIX_EQ(A, put(a, {1, 2, 3, 4}, {2}, pick(b, {1, 2, 3, 4}, {5})));
    CHECK_MATRIX_EQ(B, put(b, {1, 2, 3, 4}, {5}, pick(a, {1, 2, 3, 4}, {2})));
}

TEST_CASE("Views: swap of overlapping views, of a view with itself, and of other dimensions")
{
    M A = sequence<M>(3, 2);
    Ref a = Ref::from(A);
    swap(A(range(1, 2), all), A(range(2, 3), all));                // t = rows 1-2; rows 1-2 = rows 2-3; rows 2-3 = t
    CHECK_MATRIX_EQ(A, pick(a, {2, 1, 2}, {1, 2}));

    A = sequence<M>(3, 2);
    swap(A(all, 1), A(all, 1));                                     // nothing to exchange
    CHECK_MATRIX_EQ(A, a);

    #ifdef MATL_DIMENSION_CHECK
    CHECK_MATL_ERROR(dimension_error, swap(A(all, 1), A(1, all)));
    #endif
    CHECK(true);
}

/********************************* Temporaries and output ****************************************/

TEST_CASE("Views: views of temporaries keep them alive")
{
    auto v = ones<double>(3, 4)(all, range(2, 3));                  // the temporary is owned by v
    CHECK_EQ(sum(v), 6.0);
    auto w = v(range(1, 2), 1);                                     // and shared by views of v
    CHECK_EQ(sum(w), 2.0);

    auto f = Matrix<double, 2, 2>{{1, 2}, {3, 4}}(last, all);
    CHECK_EQ(f(1), 3.0);
    CHECK_EQ(f(2), 4.0);

    auto e = sqrt(M{{4, 9}, {16, 25}}(all, 2)) + 1.0;              // a view inside a stored expression
    M R = e;
    CHECK_EQ(R(1, 1), 4.0);
    CHECK_EQ(R(2, 1), 6.0);
}

TEST_CASE("Views: output")
{
    M A = {{1, 2, 3}, {4, 5, 6}};
    std::ostringstream os;
    os << A(all, range(last, -1, 2));
    CHECK_EQ(os.str(), std::string("3  2\n6  5"));
}

/********************************************** Errors *******************************************/

TEST_CASE("Views: errors")
{
    M A = sequence<M>(3, 4);

    // A step of 0 is always reported (also without debug checks)
    CHECK_MATL_ERROR(index_error, (void)A(range(1, 0, 3), all));

    #ifdef MATL_BOUNDS_CHECK
    CHECK_MATL_ERROR(index_error, (void)A(range(0, 2), all));
    CHECK_MATL_ERROR(index_error, (void)A(range(2, 4), all));       // 4 rows: past the end
    CHECK_MATL_ERROR(index_error, (void)A(all, 5));
    CHECK_MATL_ERROR(index_error, (void)A(last - 3, all));
    CHECK_MATL_ERROR(index_error, (void)A(range(1, 2), all)(3, 1));
    CHECK_MATL_ERROR(index_error, (void)A(all, 1)(4));
    CHECK_MATL_ERROR(index_error, (void)A(last, last - 4));

    #if MATL_THROWS_ON_ERROR
    std::string message;
    try { (void)A(all, range(2, 6)); }
    catch (index_error const& e) { message = e.what(); }
    CHECK_EQ(message, std::string("column index 6 is out of range for a matrix with 4 columns (indices are 1-based)"));
    #endif
    #endif

    #ifdef MATL_DIMENSION_CHECK
    CHECK_MATL_ERROR(dimension_error, A(range(1, 2), all) = M(3, 4));
    CHECK_MATL_ERROR(dimension_error, A(all, 1) = A(1, all));
    CHECK_MATL_ERROR(dimension_error, A(all, 1) += A(1, all));
    #endif
    CHECK(true);
}
