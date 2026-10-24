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

    // As counted_allocate, but returns null instead of throwing
    [[gnu::noinline]] void* counted_allocate_nothrow(std::size_t bytes) noexcept
    {
        if (counting) ++allocations;
        return std::malloc(bytes == 0 ? 1 : bytes);
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

// All the unaligned forms are replaced, the nothrow ones too (std::stable_sort's buffer comes from
// operator new(n, std::nothrow)), so that no block is freed by a different allocator than made it
void* operator new  (std::size_t bytes) {return counted_allocate(bytes);}
void* operator new[](std::size_t bytes) {return counted_allocate(bytes);}
void* operator new  (std::size_t bytes, std::nothrow_t const&) noexcept {return counted_allocate_nothrow(bytes);}
void* operator new[](std::size_t bytes, std::nothrow_t const&) noexcept {return counted_allocate_nothrow(bytes);}
void  operator delete  (void* p) noexcept                        {std::free(p);}
void  operator delete[](void* p) noexcept                        {std::free(p);}
void  operator delete  (void* p, std::size_t) noexcept           {std::free(p);}
void  operator delete[](void* p, std::size_t) noexcept           {std::free(p);}
void  operator delete  (void* p, std::nothrow_t const&) noexcept {std::free(p);}
void  operator delete[](void* p, std::nothrow_t const&) noexcept {std::free(p);}

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
    // Two product results (n*n each), plus each product's packing buffers (at most n*n each for
    // the left and right operands, see _matl_product_kernel_.hpp)
    MemoryManager::temporary_arena().reserve(6 * n * n * sizeof(double));

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

    // Two product results, and each product's two packing buffers (see _matl_product_kernel_.hpp)
    CountAllocations counter;
    D = A * B + C * A;
    CHECK_EQ(counter.count(), 6l);
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

TEST_CASE("Heap allocations: large reductions make none once the arena has grown")
{
    // Reductions over many blocks keep the blocks' results in the temporary arena (they run on
    // several threads with OpenMP; this counts the calling thread's allocations)
    ColMajorMatrix<double> a(50000, 1), b(50000, 1);
    a = 1.5;
    b = 2.0;
    const double warm = sum(a) + dot(a, b) + norm_fro(a - b) + max(a);
    long allocations_made = 0;
    double result = 0;
    {
        CountAllocations counter;
        result = sum(a) + dot(a, b) + norm_fro(a - b) + max(a);
        allocations_made = counter.count();
    }
    CHECK_EQ(allocations_made, 0l);
    CHECK_EQ(result, warm);
}

TEST_CASE("Heap allocations: fixed-size factorizations and solves make none")
{
    // Element-level linear algebra: factorize and solve small fixed-size systems (also beyond the
    // recursion base, where views and products run on the temporary arena)
    Matrix<double, 6, 6> A6;
    Matrix<double, 24, 24> A24;
    for (index_t i = 1; i <= 6; ++i)
        for (index_t j = 1; j <= 6; ++j) A6(i, j) = (i == j) ? 8.0 : 1.0 / double(i + j);
    for (index_t i = 1; i <= 24; ++i)
        for (index_t j = 1; j <= 24; ++j) A24(i, j) = (i == j) ? 30.0 : 1.0 / double(i + j);
    const Matrix<double, 6, 1> b6 = ones<double, 6, 1>();
    const Matrix<double, 24, 1> b24 = ones<double, 24, 1>();

    const auto work = [&]
    {
        const Matrix<double, 6, 1>  x6  = lu(A6).solve(b6);
        const Matrix<double, 24, 1> x24 = lu(A24).solve(b24);
        const Matrix<double, 24, 1> y24 = cholesky(A24).solve(b24);
        return x6(1) + x24(1) + y24(1);
    };
    const double warm = work();                       // the arena grows to fit, once
    long allocations_made = 0;
    double result = 0;
    {
        CountAllocations counter;
        result = work();
        allocations_made = counter.count();
    }
    CHECK_EQ(allocations_made, 0l);
    CHECK_EQ(result, warm);
}

TEST_CASE("Heap allocations: fixed-size det, inv, eig_sym and svd make none")
{
    // Element-level work: Jacobians (det, inv), principal stresses (eig_sym), and pseudo-inverses
    const Matrix<double, 3, 3> J = {{2, 1, 0}, {0, 3, 1}, {1, 0, 4}};
    const Matrix<double, 3, 3> sigma = {{50, 30, 20}, {30, -20, -10}, {20, -10, 10}};
    Matrix<double, 6, 6> K;
    Matrix<double, 6, 4> B;
    for (index_t i = 1; i <= 6; ++i)
    {
        for (index_t j = 1; j <= 6; ++j) K(i, j) = (i == j) ? 8.0 : 1.0 / double(i + j);
        for (index_t j = 1; j <= 4; ++j) B(i, j) = double(i * j % 5) - 2.0;
    }
    const Matrix<double, 6, 1> b = ones<double, 6, 1>();

    const auto work = [&]
    {
        const Matrix<double, 3, 3> Ji = inv(J);
        const Matrix<double, 6, 6> Ki = inv(K);
        const auto P = eig_sym(sigma);
        const auto E = eig_sym(K);
        const auto S = svd(B);
        const Matrix<double, 4, 6> Bp = S.pinv();
        const Matrix<double, 4, 1> x = S.solve(b);
        const Matrix<double, 4, 6> Bt = transpose(B);
        const Matrix<double, 4, 1> wide = singular_values(Bt);              // a wide matrix
        return det(J) + det(K) + Ji(1, 1) + Ki(1, 1) + P.values()(1) + E.vectors()(1, 1) + Bp(1, 1) + x(1) +
               wide(1) + double(S.rank());
    };
    const double warm = work();
    long allocations_made = 0;
    double result = 0;
    {
        CountAllocations counter;
        result = work();
        allocations_made = counter.count();
    }
    CHECK_EQ(allocations_made, 0l);
    CHECK_EQ(result, warm);
}
