/**************************************************************************************************
 * File: bench_temporaries.cpp
 *
 * Description: Measures the cost of the temporary matrices created during expression evaluation,
 *              with temporaries allocated from the recycled TemporaryArena (the default) and from
 *              the heap (MemoryManager::CACHE_TEMPORARY_MATRICES = false).
 *
 *              For each expression and matrix size it reports the time per assignment and the
 *              number of heap allocations per assignment.
 *
 * Build and run (from 1.0.0/):   make bench
 *
 * Note: compiled with MATL_NO_DEBUG (no bounds/dimension checks or trace output) and -O2.
 *       Timings depend on the machine; compare the two columns of a row, not across machines.
 **************************************************************************************************/
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <new>
#include <random>
#include <string>
#include <vector>

#include "matl.hpp"

using namespace matl;

/********************************** Heap allocation counter **********************************/

namespace
{
    long allocations = 0;
    void* counted_allocate(std::size_t bytes)
    {
        ++allocations;
        if (void* p = std::malloc(bytes == 0 ? 1 : bytes)) return p;
        throw std::bad_alloc();
    }
}

void* operator new  (std::size_t bytes) {return counted_allocate(bytes);}
void* operator new[](std::size_t bytes) {return counted_allocate(bytes);}
void  operator delete  (void* p) noexcept              {std::free(p);}
void  operator delete[](void* p) noexcept              {std::free(p);}
void  operator delete  (void* p, std::size_t) noexcept {std::free(p);}
void  operator delete[](void* p, std::size_t) noexcept {std::free(p);}

/*************************************** Benchmark ******************************************/

namespace
{
    struct Result
    {
        double ns_per_assignment;
        double allocations_per_assignment;
    };

    struct Operands
    {
        RowMajorMatrix<double> A, B, C, E, D;
        explicit Operands(index_t n) : A(n, n), B(n, n), C(n, n), E(n, n), D(n, n)
        {
            std::mt19937 rng(42);
            std::uniform_real_distribution<double> dist(-1.0, 1.0);
            for (RowMajorMatrix<double>* M : {&A, &B, &C, &E, &D})
                for (index_t i = 1; i <= M->size(); ++i) (*M)(i) = dist(rng);
        }
    };

    typedef std::function<void(Operands&)> Statement;

    // Time 'statement' with the arena on or off. The iteration count targets ~0.2 s per case.
    Result run(Statement const& statement, index_t n, bool use_arena)
    {
        MemoryManager::CACHE_TEMPORARY_MATRICES = use_arena;
        Operands ops(n);
        statement(ops);                                   // warm up

        const double flops_per_iteration = 2.0 * double(n) * double(n) * double(n);
        long iterations = long(2.0e8 / flops_per_iteration);
        if (iterations < 5)      iterations = 5;
        if (iterations > 200000) iterations = 200000;

        allocations = 0;
        const auto start = std::chrono::steady_clock::now();
        for (long k = 0; k < iterations; ++k) statement(ops);
        const auto stop  = std::chrono::steady_clock::now();
        const long counted = allocations;

        volatile double sink = ops.D(1,1);                // keep the result observable
        (void)sink;

        const double ns = std::chrono::duration<double, std::nano>(stop - start).count();
        return Result{ns / double(iterations), double(counted) / double(iterations)};
    }
}

int main()
{
    struct Case {const char* text; Statement statement;};
    const std::vector<Case> cases = {
        {"D = A*B + C",                   [](Operands& o){o.D = o.A*o.B + o.C;}},
        {"D = A*B + transpose(C*A) - E",  [](Operands& o){o.D = o.A*o.B + transpose(o.C*o.A) - o.E;}},
        {"D = (A + B)*(C - E)",           [](Operands& o){o.D = (o.A + o.B)*(o.C - o.E);}},
        {"D = A*D",                       [](Operands& o){o.D = o.A*o.D;}},
        {"D = A*B  (no temporaries)",     [](Operands& o){o.D = o.A*o.B;}},
    };
    const std::vector<index_t> sizes = {4, 16, 64, 128};

    std::printf("MATL temporaries benchmark: arena (CACHE_TEMPORARY_MATRICES = true) vs heap\n\n");
    std::printf("%-32s %5s | %12s %12s | %8s %8s | %7s\n",
                "expression", "n", "arena ns", "heap ns", "arena", "heap", "speedup");
    std::printf("%-32s %5s | %12s %12s | %8s %8s | %7s\n",
                "", "", "/assign", "/assign", "allocs", "allocs", "");
    std::printf("%s\n", std::string(100, '-').c_str());

    for (const Case& c : cases)
    {
        for (index_t n : sizes)
        {
            const Result arena = run(c.statement, n, true);
            const Result heap  = run(c.statement, n, false);
            std::printf("%-32s %5lu | %12.0f %12.0f | %8.2f %8.2f | %6.2fx\n",
                        c.text, n, arena.ns_per_assignment, heap.ns_per_assignment,
                        arena.allocations_per_assignment, heap.allocations_per_assignment,
                        heap.ns_per_assignment / arena.ns_per_assignment);
        }
    }
    MemoryManager::CACHE_TEMPORARY_MATRICES = true;
    return 0;
}
