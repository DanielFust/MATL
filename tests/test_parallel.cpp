/**************************************************************************************************
 * File: test_parallel.cpp
 *
 * Description: Thread parallelism (OpenMP): parallel_ranges covers every item once, exceptions
 *              inside parallel loops are rethrown, results do not depend on the number of threads
 *              (bit for bit: reductions, elementwise evaluation, sparse products, and iterative
 *              solves), and MATL runs serially inside a caller's parallel region. Built without
 *              OpenMP, the same tests check the serial code paths (make check-openmp runs them on
 *              4 threads).
 **************************************************************************************************/
#include <cmath>
#include <stdexcept>
#include <vector>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using V = ColMajorMatrix<double>;

    // A vector large enough to run in parallel, with values in (-1, 1)
    V large_vector(index_t n)
    {
        V v(n, 1);
        std::uniform_real_distribution<double> u(-1.0, 1.0);
        for (index_t i = 1; i <= n; ++i) v(i) = u(matl_test::rng());
        return v;
    }

    CSRMatrix<double> laplacian_2d(index_t k)
    {
        SparseBuilder<double> B(k * k, k * k);
        B.reserve(5 * k * k);
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

    // Runs f with the given number of OpenMP threads (no effect without OpenMP)
    template<typename F>
    void with_threads([[maybe_unused]] int threads, F f)
    {
        #if MATL_OPENMP
        const int previous = omp_get_max_threads();
        omp_set_num_threads(threads);
        f();
        omp_set_num_threads(previous);
        #else
        f();
        #endif
    }

    // Results of a computation that exercises the parallel paths
    struct Results
    {
        double sum, dot, norm, mean, min, max, prod;
        V      elementwise, product;
        index_t iterations;
        double  residual;
        V       solution;
    };

    Results compute(V const& a, V const& b, CSRMatrix<double> const& L, V const& rhs)
    {
        Results r;
        r.sum  = sum(a);
        r.dot  = dot(a, b);
        r.norm = norm_fro(a - 2.0 * b);
        r.mean = mean(abs(a));
        r.min  = min(a + b);
        r.max  = max(a);
        r.prod = prod(1.0 + 1e-6 * a);
        r.elementwise = 2.0 * a - sqrt(abs(b)) + times(a, b);
        r.product = L * rhs;
        r.solution = zeros<double>(rhs.rows(), 1);
        const auto result = cg(L, rhs, r.solution, {.tolerance = 1e-10}, JacobiPreconditioner(L));
        r.iterations = result.iterations;
        r.residual = result.residual;
        return r;
    }
}

TEST_CASE("Parallel: parallel_ranges covers every item once")
{
    for (index_t n : {index_t(0), index_t(1), index_t(7), 3 * detail::PARALLEL_THRESHOLD + 5})
    {
        std::vector<int> visits(n, 0);
        detail::parallel_ranges(n, detail::PARALLEL_THRESHOLD * 4, [&](index_t first, index_t last)
        {
            for (index_t k = first; k < last; ++k) ++visits[k];
        });
        bool once = true;
        for (int v : visits) once = once && (v == 1);
        CHECK(once);
    }

    // With a custom split (e.g. by the stored entries of a sparse matrix)
    index_t n = 100000;                                         // (not const: the lambda below captures it)
    std::vector<int> visits(n, 0);
    detail::parallel_ranges(n, n, [&](index_t first, index_t last)
    {
        for (index_t k = first; k < last; ++k) ++visits[k];
    }, [n](index_t t, index_t threads) {return n * t * t / (threads * threads);});
    bool once = true;
    for (int v : visits) once = once && (v == 1);
    CHECK(once);
}

TEST_CASE("Parallel: an exception inside a parallel loop is rethrown after it")
{
    #if MATL_PARALLEL_EXCEPTIONS
    const index_t n = 4 * detail::PARALLEL_THRESHOLD;
    CHECK_THROWS_AS(detail::parallel_ranges(n, n, [](index_t first, index_t last)
    {
        if (first <= 12345 && 12345 < last) throw std::runtime_error("item 12345");
    }), std::runtime_error);

    // MATL is usable afterwards
    V a = large_vector(n);
    CHECK(std::isfinite(sum(a)));
    #endif
    CHECK(true);
}

TEST_CASE("Parallel: results do not depend on the number of threads")
{
    const index_t n = 3 * detail::PARALLEL_THRESHOLD;
    const V a = large_vector(n), b = large_vector(n);
    const CSRMatrix<double> L = laplacian_2d(200);                 // n = 40000, 199200 stored entries
    const V rhs = large_vector(L.rows());

    Results one, four;
    with_threads(1, [&] {one = compute(a, b, L, rhs);});
    with_threads(4, [&] {four = compute(a, b, L, rhs);});

    // Bit for bit (the reductions add fixed blocks in a fixed order)
    CHECK_EQ(one.sum, four.sum);
    CHECK_EQ(one.dot, four.dot);
    CHECK_EQ(one.norm, four.norm);
    CHECK_EQ(one.mean, four.mean);
    CHECK_EQ(one.min, four.min);
    CHECK_EQ(one.max, four.max);
    CHECK_EQ(one.prod, four.prod);
    CHECK_MATRIX_NEAR(four.elementwise, Ref::from(one.elementwise), 0.0);
    CHECK_MATRIX_NEAR(four.product, Ref::from(one.product), 0.0);
    CHECK_EQ(one.iterations, four.iterations);
    CHECK_EQ(one.residual, four.residual);
    CHECK_MATRIX_NEAR(four.solution, Ref::from(one.solution), 0.0);

    // And the values are right
    double s = 0, d = 0;
    for (index_t i = 1; i <= n; ++i) {s += a(i); d += a(i) * b(i);}
    CHECK(std::fabs(one.sum - s) < 1e-9);
    CHECK(std::fabs(one.dot - d) < 1e-9);
    CHECK(one.residual <= 1e-10);
}

TEST_CASE("Parallel: MATL runs serially inside a caller's parallel region")
{
    const index_t n = 2 * detail::PARALLEL_THRESHOLD;
    const V a = large_vector(n);
    const double expected = sum(2.0 * a);

    std::vector<double> results(4, 0.0);
    std::size_t team = results.size();
    #if MATL_OPENMP
    #pragma omp parallel num_threads(4)
    {
        #pragma omp single
        team = static_cast<std::size_t>(omp_get_num_threads());
        const int t = omp_get_thread_num();
        V twice = 2.0 * a;                                         // each thread its own evaluation
        results[static_cast<std::size_t>(t)] = sum(twice);
    }
    #else
    for (double& r : results) r = sum(2.0 * a);
    #endif
    for (std::size_t t = 0; t < team; ++t) CHECK_EQ(results[t], expected);
}
