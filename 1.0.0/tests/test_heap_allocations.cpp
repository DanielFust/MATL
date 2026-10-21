/**************************************************************************************************
 * File: test_heap_allocations.cpp
 *
 * Description: Counts heap allocations during expression evaluation. With the temporary arena
 *              (the default), assignments that need temporaries make no heap allocations once
 *              the arena has grown to fit them.
 *
 *              This file replaces the global operator new/delete for the test program. The
 *              replacements only count allocations while a CountAllocations object is active,
 *              and otherwise behave like the standard ones.
 **************************************************************************************************/
#include <cstdlib>
#include <new>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    thread_local bool counting    = false;
    thread_local long allocations = 0;

    // Not inlined: when GCC 12/13 inline a malloc-based operator new into std::vector's growth
    // path at -O2, they report a false-positive -Wstringop-overflow error
    [[gnu::noinline]] void* counted_allocate(std::size_t bytes)
    {
        if (counting) ++allocations;
        if (void* p = std::malloc(bytes == 0 ? 1 : bytes)) return p;
        throw std::bad_alloc();
    }

    // Counts heap allocations made by this thread while the object is alive
    struct CountAllocations
    {
        CountAllocations()  {allocations = 0; counting = true;}
        ~CountAllocations() {counting = false;}
        long count() const  {return allocations;}
    };

    struct CacheTemporaries
    {
        bool previous;
        explicit CacheTemporaries(bool value) : previous(MemoryManager::CACHE_TEMPORARY_MATRICES)
        {MemoryManager::CACHE_TEMPORARY_MATRICES = value;}
        ~CacheTemporaries() {MemoryManager::CACHE_TEMPORARY_MATRICES = previous;}
    };
}

void* operator new  (std::size_t bytes) {return counted_allocate(bytes);}
void* operator new[](std::size_t bytes) {return counted_allocate(bytes);}
void  operator delete  (void* p) noexcept              {std::free(p);}
void  operator delete[](void* p) noexcept              {std::free(p);}
void  operator delete  (void* p, std::size_t) noexcept {std::free(p);}
void  operator delete[](void* p, std::size_t) noexcept {std::free(p);}

TEST_CASE("Heap allocations: the counter works")
{
    // Call the allocation function directly: a 'new' expression may legally be optimized away
    CountAllocations counter;
    void* p = ::operator new(sizeof(int));
    ::operator delete(p);
    std::byte* q = new std::byte[64];
    q[0] = std::byte(1);
    CHECK(q[0] == std::byte(1));
    delete[] q;
    CHECK(counter.count() >= 1l);
}

TEST_CASE("Heap allocations: none during evaluation once the arena has grown")
{
    const index_t n = 16;
    RowMajorMatrix<double> A(n, n), B(n, n), C(n, n), D(n, n);
    matl_test::fill_random(A); matl_test::fill_random(B); matl_test::fill_random(C);

    // Warm up: the arena grows to fit the temporaries (and stream buffers are set up)
    D = A * B + transpose(C * A) - (A + B) * C;
    D = A * D;

    CountAllocations counter;
    for (int k = 0; k < 10; ++k)
    {
        D = A * B + transpose(C * A) - (A + B) * C;   // three products, two expression operands
        D = A * D;                                    // aliasing temporary
        D = (A * B) * (C * A) + D;
    }
    CHECK_EQ(counter.count(), 0l);
}

TEST_CASE("Heap allocations: reserve() avoids allocation on the first evaluation")
{
    const index_t n = 64;
    RowMajorMatrix<double> A(n, n), B(n, n), D(n, n);
    matl_test::fill_random(A); matl_test::fill_random(B);
    MemoryManager::temporary_arena().reserve(4 * n * n * sizeof(double));

    CountAllocations counter;
    D = A * B + B * A;
    CHECK_EQ(counter.count(), 0l);
}

TEST_CASE("Heap allocations: without the arena, each temporary is a heap allocation")
{
    const index_t n = 16;
    RowMajorMatrix<double> A(n, n), B(n, n), C(n, n), D(n, n);
    matl_test::fill_random(A); matl_test::fill_random(B); matl_test::fill_random(C);
    CacheTemporaries off(false);
    D = A * B + C * A;                                // warm up

    CountAllocations counter;
    D = A * B + C * A;                                // two temporaries
    CHECK_EQ(counter.count(), 2l);
}

TEST_CASE("Heap allocations: building and copying expressions allocates nothing")
{
    RowMajorMatrix<double> A(4, 4), B(4, 4);
    A = 1.0;
    B = 2.0;
    CountAllocations counter;
    auto e = A * B + transpose(A * B) - (A + B) * 2.0;
    auto f = e;
    (void)f;
    CHECK_EQ(counter.count(), 0l);
}
