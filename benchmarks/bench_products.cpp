/**************************************************************************************************
 * File: bench_products.cpp
 *
 * Description: Measures matrix products and transposes: the time per assignment and, for
 *              products, the floating-point rate (a product of an MxK and a KxN matrix is
 *              2*M*N*K floating-point operations).
 *
 * Build and run:   make bench
 *
 * Note: compiled with -O2 -DNDEBUG (a release build: no bounds or dimension checks). Timings are
 *       the best of several runs and depend on the machine.
 **************************************************************************************************/
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>

#include "matl.hpp"

using namespace matl;

namespace
{
    // Best time per call of f, in nanoseconds
    template<typename F>
    double best_ns(F&& f, double flops_per_call)
    {
        // Enough calls per run for about 50 ms of work (assuming at least ~1 GFLOP/s)
        const long calls = std::max(1L, static_cast<long>(5.0e7 / std::max(flops_per_call, 1.0)));
        double best = 1e300;
        for (int run = 0; run < 5; ++run)
        {
            const auto start = std::chrono::steady_clock::now();
            for (long c = 0; c < calls; ++c) f();
            const auto stop = std::chrono::steady_clock::now();
            best = std::min(best, std::chrono::duration<double, std::nano>(stop - start).count() / calls);
        }
        return best;
    }

    template<typename MAT>
    void fill(MAT& m, std::mt19937& rng)
    {
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        for (index_t i = 1; i <= m.size(); ++i) m(i) = dist(rng);
    }

    void row(char const* label, index_t n, double ns, double flops)
    {
        if (flops > 0) std::printf("%-30s %5lu | %14.0f | %8.2f\n", label, n, ns, flops / ns);
        else           std::printf("%-30s %5lu | %14.0f |\n", label, n, ns);
    }

    // C = A*B and D = A*B + C for n x n dynamic matrices
    void dynamic_products(index_t n, std::mt19937& rng)
    {
        RowMajorMatrix<double> A(n, n), B(n, n), C(n, n), D(n, n);
        fill(A, rng); fill(B, rng); fill(C, rng);
        const double flops = 2.0 * n * n * n;
        row("C = A*B",     n, best_ns([&]{ C = A * B; },     flops), flops);
        row("D = A*B + C", n, best_ns([&]{ D = A * B + C; }, flops), flops);
    }

    // C = A*B for N x N fixed-size matrices
    template<index_t N>
    void fixed_product(std::mt19937& rng)
    {
        FixedSizeRowMajorMatrix<double, N, N> A, B, C;
        fill(A, rng); fill(B, rng);
        const double flops = 2.0 * N * N * N;
        row("C = A*B (fixed-size)", N, best_ns([&]{ C = A * B; }, flops), flops);
    }

    // Transposes: elements are read in a different order than they are written
    void transposes(index_t n, std::mt19937& rng)
    {
        RowMajorMatrix<double> A(n, n), B(n, n), C(n, n), R(n, 2 * n), T(2 * n, n);
        fill(A, rng); fill(B, rng); fill(R, rng);
        const double work = 1.0 * n * n;
        row("B = transpose(A)",     n, best_ns([&]{ B = transpose(A); },     work), 0);
        row("C = A + transpose(B)", n, best_ns([&]{ C = A + transpose(B); }, work), 0);
        row("T = transpose(R) (2n x n)", n, best_ns([&]{ T = transpose(R); }, work), 0);
    }
}

int main()
{
    std::mt19937 rng(42);
    std::printf("MATL products and transposes (best of 5 runs)\n\n");
    std::printf("%-30s %5s | %14s | %8s\n", "operation", "n", "ns / assign", "GFLOP/s");
    std::printf("-------------------------------------------------------------------\n");
    for (index_t n : {4, 16, 64, 128, 256, 512, 1024}) dynamic_products(n, rng);
    fixed_product<4>(rng);
    fixed_product<8>(rng);
    fixed_product<16>(rng);
    for (index_t n : {64, 512, 2048}) transposes(n, rng);
    return 0;
}
