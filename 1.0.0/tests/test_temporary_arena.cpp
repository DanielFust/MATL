/**************************************************************************************************
 * File: test_temporary_arena.cpp
 *
 * Description: TemporaryArena (stack-ordered memory for evaluation temporaries) and how
 *              expression evaluation uses it: which assignments need temporaries, and that the
 *              memory is released and recycled.
 **************************************************************************************************/
#include <cstdint>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    bool aligned(std::byte* p)
    {return reinterpret_cast<std::uintptr_t>(p) % TemporaryArena::ALIGNMENT == 0;}

    // Bytes the arena uses for an n-element block of doubles (blocks are rounded to ALIGNMENT)
    index_t block_bytes(index_t elements)
    {
        const index_t a = TemporaryArena::ALIGNMENT;
        return ((elements * sizeof(double) + a - 1) / a) * a;
    }

    // Measures the arena's high water mark across one statement
    struct ArenaPeak
    {
        ArenaPeak()  {MemoryManager::temporary_arena().reset_high_water_mark();}
        index_t peak() const {return MemoryManager::temporary_arena().high_water_mark();}
    };

    // Temporarily changes MemoryManager::CACHE_TEMPORARY_MATRICES
    struct CacheTemporaries
    {
        bool previous;
        explicit CacheTemporaries(bool value) : previous(MemoryManager::CACHE_TEMPORARY_MATRICES)
        {MemoryManager::CACHE_TEMPORARY_MATRICES = value;}
        ~CacheTemporaries() {MemoryManager::CACHE_TEMPORARY_MATRICES = previous;}
    };
}

/************************************* TemporaryArena ***************************************/

TEST_CASE("TemporaryArena: blocks are aligned, distinct and writable")
{
    TemporaryArena arena;
    CHECK_EQ(arena.capacity(), 0ul);              // no memory until first use
    std::byte* a = arena.allocate(24);
    std::byte* b = arena.allocate(0);
    std::byte* c = arena.allocate(1000);
    CHECK(aligned(a) && aligned(b) && aligned(c));
    CHECK(a != b && b != c && a != c);
    for (int k = 0; k < 1000; ++k) c[k] = std::byte(k & 0xff);
    for (int k = 0; k < 24; ++k)   a[k] = std::byte(1);
    CHECK(c[999] == std::byte(999 & 0xff));       // blocks do not overlap
    CHECK_EQ(arena.blocks_in_use(), 3ul);
    arena.release(c);
    arena.release(b);
    arena.release(a);
    CHECK_EQ(arena.bytes_in_use(), 0ul);
    CHECK_EQ(arena.blocks_in_use(), 0ul);
}

TEST_CASE("TemporaryArena: memory is reused in stack order")
{
    TemporaryArena arena;
    std::byte* a = arena.allocate(256);
    arena.release(a);
    std::byte* b = arena.allocate(256);
    CHECK(a == b);                                // the same memory is handed out again
    arena.release(b);
}

TEST_CASE("TemporaryArena: grows when needed and consolidates when empty")
{
    TemporaryArena arena;
    const index_t big = 3 * TemporaryArena::MIN_CHUNK_SIZE;
    std::byte* a = arena.allocate(TemporaryArena::MIN_CHUNK_SIZE / 2);
    std::byte* b = arena.allocate(big);           // does not fit: a new chunk is added
    std::byte* c = arena.allocate(big);
    b[big - 1] = std::byte(7);
    c[0]       = std::byte(9);
    CHECK(b[big - 1] == std::byte(7));
    const index_t grown = arena.capacity();
    CHECK(grown >= TemporaryArena::MIN_CHUNK_SIZE / 2 + 2 * big);
    arena.release(c);
    arena.release(b);
    arena.release(a);

    // Merged into a single chunk of the same capacity: the same allocations now fit without
    // growing, in one contiguous chunk
    CHECK_EQ(arena.capacity(), grown);
    std::byte* a2 = arena.allocate(TemporaryArena::MIN_CHUNK_SIZE / 2);
    std::byte* b2 = arena.allocate(big);
    std::byte* c2 = arena.allocate(big);
    CHECK_EQ(arena.capacity(), grown);
    CHECK(b2 > a2 && c2 > b2);
    arena.release(c2);
    arena.release(b2);
    arena.release(a2);
}

TEST_CASE("TemporaryArena: out-of-order release is reclaimed when the blocks above are released")
{
    TemporaryArena arena;
    std::byte* a = arena.allocate(64);
    std::byte* b = arena.allocate(64);
    std::byte* c = arena.allocate(64);

    arena.release(b);                             // out of order: marked free
    CHECK_EQ(arena.bytes_in_use(), 128ul);
    CHECK_EQ(arena.blocks_in_use(), 3ul);         // not reclaimed while c is in use

    std::byte* d = arena.allocate(64);            // allocated above c, not in b's memory
    CHECK(d > c);
    arena.release(d);
    arena.release(c);                             // reclaims c and b
    CHECK_EQ(arena.blocks_in_use(), 1ul);
    std::byte* e = arena.allocate(64);
    CHECK(e == b);                                // b's memory is available again
    arena.release(e);
    arena.release(a);
    CHECK_EQ(arena.bytes_in_use(), 0ul);
}

TEST_CASE("TemporaryArena: invalid release is an error")
{
    TemporaryArena arena;
    std::byte* a = arena.allocate(16);
    std::byte not_from_arena[16];
    CHECK_THROWS_AS(arena.release(not_from_arena), std::logic_error);
    arena.release(a);
    CHECK_THROWS_AS(arena.release(a), std::logic_error);   // double release
}

TEST_CASE("TemporaryArena: reserve, statistics and release_memory")
{
    TemporaryArena arena;
    arena.reserve(1 << 20);
    CHECK(arena.capacity() >= index_t(1 << 20));
    const index_t reserved = arena.capacity();
    std::byte* a = arena.allocate(1 << 19);
    std::byte* b = arena.allocate(1 << 19);
    CHECK_EQ(arena.capacity(), reserved);         // fits in the reservation
    CHECK_EQ(arena.high_water_mark(), index_t(1 << 20));
    CHECK_THROWS_AS(arena.release_memory(), std::logic_error);  // still in use
    arena.release(b);
    arena.release(a);
    arena.reset_high_water_mark();
    CHECK_EQ(arena.high_water_mark(), 0ul);
    arena.release_memory();
    CHECK_EQ(arena.capacity(), 0ul);
}

/*************************** Temporaries used by evaluation ***********************************/

TEST_CASE("Arena: assignments that need no temporaries")
{
    const index_t n = 6;
    RowMajorMatrix<double> A(n, n), B(n, n), C(n, n), D(n, n);
    matl_test::fill_random(A); matl_test::fill_random(B); matl_test::fill_random(C);
    { ArenaPeak p; D = A * B;                     CHECK_EQ(p.peak(), 0ul); }  // product into target
    { ArenaPeak p; D = A + B * 2.0 - C;           CHECK_EQ(p.peak(), 0ul); }  // elementwise
    { ArenaPeak p; D = transpose(A) + C;          CHECK_EQ(p.peak(), 0ul); }
    { ArenaPeak p; D = A; D = transpose(D);       CHECK_EQ(p.peak(), 0ul); }  // in-place square transpose
    { ArenaPeak p; D = A * C; D = D + A * 2.0;    CHECK_EQ(p.peak(), 0ul); }  // target read elementwise
}

TEST_CASE("Arena: assignments that need temporaries use exactly what they need")
{
    const index_t n = 6;
    const index_t one = block_bytes(n * n);
    RowMajorMatrix<double> A(n, n), B(n, n), C(n, n), D(n, n);
    matl_test::fill_random(A); matl_test::fill_random(B); matl_test::fill_random(C);
    { ArenaPeak p; D = A * B + C;                 CHECK_EQ(p.peak(), one); }      // nested product
    { ArenaPeak p; D = A * B + transpose(C * A);  CHECK_EQ(p.peak(), 2 * one); }  // two nested products
    { ArenaPeak p; D = (A + B) * C;               CHECK_EQ(p.peak(), one); }      // expression operand
    { ArenaPeak p; D = (A + B) * (C - A);         CHECK_EQ(p.peak(), 2 * one); }
    { ArenaPeak p; D = A * D;                     CHECK_EQ(p.peak(), one); }      // target is an operand
    { ArenaPeak p; D = D + transpose(D);          CHECK_EQ(p.peak(), one); }      // target read transposed
    CHECK_EQ(MemoryManager::temporary_arena().bytes_in_use(), 0ul);
}

TEST_CASE("Arena: operands of a product are released before the product's result")
{
    // (A+B)*(C-A) nested in a sum: the product's result (1 block) is allocated first, then its two
    // expression operands (2 blocks), which are released as soon as the product is computed.
    const index_t n = 5;
    const index_t one = block_bytes(n * n);
    RowMajorMatrix<double> A(n, n), B(n, n), C(n, n), D(n, n);
    matl_test::fill_random(A); matl_test::fill_random(B); matl_test::fill_random(C);
    ArenaPeak p;
    D = (A + B) * (C - A) + C * A;                // second product allocated after the first's operands are gone
    CHECK_EQ(p.peak(), 3 * one);
    Ref a = Ref::from(A), b = Ref::from(B), c = Ref::from(C);
    CHECK_MATRIX_EQ(D, (a + b) * (c - a) + c * a);
}

TEST_CASE("Arena: memory is recycled across assignments")
{
    const index_t n = 8;
    RowMajorMatrix<double> A(n, n), B(n, n), C(n, n), D(n, n);
    matl_test::fill_random(A); matl_test::fill_random(B); matl_test::fill_random(C);
    TemporaryArena& arena = MemoryManager::temporary_arena();

    D = A * B + transpose(C * A) - (A + B) * C;   // warm up
    const index_t capacity = arena.capacity();
    for (int k = 0; k < 20; ++k)
    {
        D = A * B + transpose(C * A) - (A + B) * C;
        CHECK_EQ(arena.bytes_in_use(), 0ul);
    }
    CHECK_EQ(arena.capacity(), capacity);         // no growth after the first evaluation
    Ref a = Ref::from(A), b = Ref::from(B), c = Ref::from(C);
    CHECK_MATRIX_EQ(D, a * b + transpose(c * a) - (a + b) * c);
}

TEST_CASE("Arena: fixed-size temporaries need no arena memory")
{
    FixedSizeRowMajorMatrix<double, 3, 3> A, B, C, D;
    matl_test::fill_random(A); matl_test::fill_random(B); matl_test::fill_random(C);
    ArenaPeak p;
    D = A * B + C * A;
    CHECK_EQ(p.peak(), 0ul);
    CHECK_MATRIX_EQ(D, Ref::from(A) * Ref::from(B) + Ref::from(C) * Ref::from(A));
}

TEST_CASE("Arena: CACHE_TEMPORARY_MATRICES = false allocates temporaries on the heap")
{
    CacheTemporaries off(false);
    const index_t n = 6;
    RowMajorMatrix<double> A(n, n), B(n, n), C(n, n), D(n, n);
    matl_test::fill_random(A); matl_test::fill_random(B); matl_test::fill_random(C);
    ArenaPeak p;
    D = A * B + transpose(C * A);
    CHECK_EQ(p.peak(), 0ul);                      // arena unused
    Ref a = Ref::from(A), b = Ref::from(B), c = Ref::from(C);
    CHECK_MATRIX_EQ(D, a * b + transpose(c * a));
    D = A * D;                                    // aliasing temporary on the heap
    CHECK_MATRIX_EQ(D, a * (a * b + transpose(c * a)));
}

#ifdef MATL_BOUNDS_CHECK
TEST_CASE("Arena: an exception during evaluation releases the temporaries")
{
    // A (3x3) times B (2x2) has mismatched inner dimensions. The nested product's temporary is
    // allocated, then reading B(3,j) during the product throws (fixed-size bounds check).
    RowMajorMatrix<double> A(3, 3), C(3, 2), R(3, 2);
    FixedSizeRowMajorMatrix<double, 2, 2> B;
    A = 1.0;
    B = 2.0;
    C = 0.0;
    ArenaPeak p;
    CHECK_THROWS_AS(R = A * B + C, std::out_of_range);
    CHECK(p.peak() > 0);                                     // a temporary was in use...
    CHECK_EQ(MemoryManager::temporary_arena().bytes_in_use(), 0ul);   // ...and was released
}
#endif
