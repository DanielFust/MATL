/**************************************************************************************************
 * File: test_reductions.cpp
 *
 * Description: Reductions (sum, prod, mean, min, max, trace, dot, norm_fro, norm_1, norm_inf) of
 *              matrices and expressions, compared with the reference implementation.
 **************************************************************************************************/
#include <algorithm>
#include <cmath>
#include <type_traits>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using M = RowMajorMatrix<double>;

    double ref_sum(Ref const& r) {double s = 0; for (double x : r.data) s += x; return s;}
    double ref_prod(Ref const& r) {double p = 1; for (double x : r.data) p *= x; return p;}
    double ref_min(Ref const& r) {return *std::min_element(r.data.begin(), r.data.end());}
    double ref_max(Ref const& r) {return *std::max_element(r.data.begin(), r.data.end());}
    double ref_trace(Ref const& r) {double t = 0; for (index_t i = 1; i <= r.rows; ++i) t += r(i, i); return t;}

    M random_matrix(index_t m, index_t n)
    {
        M A(m, n);
        matl_test::fill_random(A);
        return A;
    }
}

TEST_CASE("Reductions: sum, prod, mean, min and max")
{
    M A = random_matrix(4, 5);
    Ref a = Ref::from(A);
    CHECK_EQ(sum(A), ref_sum(a));
    CHECK_EQ(mean(A), ref_sum(a) / 20.0);
    CHECK_EQ(min(A), ref_min(a));
    CHECK_EQ(max(A), ref_max(a));

    Matrix<double> P = {{1, -2, 0.5}, {3, 4, -1}};
    CHECK_EQ(prod(P), ref_prod(Ref::from(P)));
    CHECK_EQ(prod(P), 12.0);
}

TEST_CASE("Reductions: expressions are reduced without a temporary")
{
    M A = random_matrix(4, 3), B = random_matrix(3, 4), C = random_matrix(4, 3);
    Ref a = Ref::from(A), b = Ref::from(B), c = Ref::from(C);

    MemoryManager::temporary_arena().reset_high_water_mark();
    CHECK_EQ(sum(A + C * 2.0), ref_sum(a + c * 2.0));
    CHECK_EQ(max(abs(A - C) * 3.0), 3.0 * ref_max(Ref::from(M(abs(A - C)))));
    CHECK_EQ(sum(A + transpose(B)), ref_sum(a + transpose(b)));     // mixed storage orders
    CHECK_EQ(min(transpose(B) - A), ref_min(transpose(b) - a));
    CHECK_EQ(MemoryManager::temporary_arena().high_water_mark(), 0ul);

    // Products are computed into a temporary, as in any evaluation
    CHECK_EQ(sum(A * B), ref_sum(a * b));
    CHECK_EQ(trace(B * A), ref_trace(b * a));
    CHECK_EQ(MemoryManager::temporary_arena().bytes_in_use(), 0ul);
}

TEST_CASE("Reductions: trace and dot")
{
    Matrix<double> A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    CHECK_EQ(trace(A), 15.0);
    CHECK_EQ(trace(transpose(A) + A), 30.0);
    CHECK_EQ(trace(Matrix<double>()), 0.0);

    Matrix<double> x = {{1}, {2}, {3}}, y = {{4}, {-5}, {6}};
    CHECK_EQ(dot(x, y), 12.0);                     // vectors
    CHECK_EQ(dot(A, A), 285.0);                    // Frobenius inner product
    CHECK_EQ(dot(transpose(x), transpose(y)), 12.0);

    Matrix<double, 2, 2> F = {{1, 2}, {3, 4}}, G = {{1, 0}, {0, 1}};
    CHECK_EQ(dot(F, G), 5.0);
}

TEST_CASE("Reductions: norms")
{
    Matrix<double> A = {{1, -2}, {-3, 4}, {5, -6}};
    CHECK_EQ(norm_fro(A), std::sqrt(91.0));
    CHECK_EQ(norm_1(A), 12.0);                     // columns: 9 and 12
    CHECK_EQ(norm_inf(A), 11.0);                   // rows: 3, 7 and 11
    CHECK_EQ(norm_1(transpose(A)), 11.0);          // the transpose swaps the two
    CHECK_EQ(norm_inf(transpose(A)), 12.0);

    Matrix<double> x = {{3}, {-4}};
    CHECK_EQ(norm_fro(x), 5.0);                    // 2-norm of a vector
    CHECK_EQ(norm_1(x), 7.0);                      // 1-norm of a vector
    CHECK_EQ(norm_inf(x), 4.0);                    // largest absolute value
}

TEST_CASE("Reductions: result types")
{
    Matrix<int> I = {{1, 2}, {3, 5}};
    static_assert(std::is_same_v<decltype(sum(I)), int>);
    static_assert(std::is_same_v<decltype(mean(I)), double>);
    static_assert(std::is_same_v<decltype(norm_fro(I)), double>);
    static_assert(std::is_same_v<decltype(norm_1(I)), int>);
    CHECK_EQ(sum(I), 11);
    CHECK_EQ(mean(I), 2.75);
    CHECK_EQ(norm_inf(-I), 8);

    Matrix<unsigned, 1, 2> U = {{3u, 7u}};
    CHECK_EQ(norm_1(U), 7u);
    CHECK_EQ(norm_inf(U), 10u);

    Matrix<float> F = {{1.5f, 2.5f}};
    static_assert(std::is_same_v<decltype(mean(F)), float>);
    CHECK_EQ(mean(F), 2.0f);
}

TEST_CASE("Reductions: empty matrices")
{
    Matrix<double> E;
    CHECK_EQ(sum(E), 0.0);
    CHECK_EQ(prod(E), 1.0);
    CHECK_EQ(norm_fro(E), 0.0);
    CHECK_EQ(norm_1(E), 0.0);
    CHECK_MATL_ERROR(dimension_error, (void)mean(E));
    CHECK_MATL_ERROR(dimension_error, (void)min(E));
    CHECK_MATL_ERROR(dimension_error, (void)max(E));
}

TEST_CASE("Reductions: dimensions are checked")
{
    #ifdef MATL_DIMENSION_CHECK
    M A(2, 3), B(3, 2), S(2, 2);
    A = 1.0; B = 1.0; S = 1.0;
    CHECK_MATL_ERROR(dimension_error, (void)trace(A));        // not square
    CHECK_MATL_ERROR(dimension_error, (void)dot(A, B));       // different dimensions
    CHECK_MATL_ERROR(dimension_error, (void)sum(A + B));      // inside the expression
    CHECK_MATL_ERROR(dimension_error, (void)norm_1(A * A));
    #endif
    CHECK(true);
}

TEST_CASE("Reductions: temporaries and stored expressions")
{
    CHECK_EQ(sum(ones<double>(3, 4)), 12.0);
    CHECK_EQ(trace(eye<double>(5) * 2.0), 10.0);

    M A = random_matrix(3, 3);
    auto e = A * 2.0;
    A = 1.0;                                       // e reads A's current values
    CHECK_EQ(sum(e), 18.0);
}
