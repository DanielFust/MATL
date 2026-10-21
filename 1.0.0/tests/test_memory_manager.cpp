/**************************************************************************************************
 * File: test_memory_manager.cpp
 *
 * Description: MemoryManager / ScratchContainer behavior: allocating matrices in scratch memory,
 *              reuse of freed memory, externally supplied scratch buffers, and teardown.
 **************************************************************************************************/
#include <memory>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    // Every container's neighbor links must be mutually consistent
    bool neighbor_links_consistent()
    {
        for (ScratchContainer& c : MemoryManager::get_memory_manager().containers())
        {
            if (c.next_neighbor() != nullptr && c.next_neighbor()->prev_neighbor() != &c) return false;
            if (c.prev_neighbor() != nullptr && c.prev_neighbor()->next_neighbor() != &c) return false;
        }
        return true;
    }
}

TEST_CASE("MemoryManager: is enabled by default and scratch allocation is opt-in")
{
    CHECK(MemoryManager::get_memory_manager().enabled());
    CHECK(!MemoryManager::ALLOCATE_MATRICES_IN_SCRATCH);
}

TEST_CASE("MemoryManager: matrices can be allocated in scratch memory")
{
    matl_test::ScratchMode scratch;
    RowMajorMatrix<double> A(4, 4);
    matl_test::fill_sequence(A);
    CHECK_MATRIX_EQ(A, Ref::from(A));
    CHECK(!A.external_mem());
    CHECK(neighbor_links_consistent());
}

TEST_CASE("MemoryManager: many scratch matrices with interleaved frees stay intact")
{
    matl_test::ScratchMode scratch;
    std::vector<std::unique_ptr<RowMajorMatrix<double>>> mats;
    for (int k = 0; k < 40; ++k)
    {
        mats.emplace_back(new RowMajorMatrix<double>(16 + k, 16));
        *mats.back() = double(k);
    }
    for (int k = 0; k < 40; ++k) CHECK_MATRIX_EQ(*mats[k], Ref(16 + k, 16, double(k)));

    // Free every other matrix, then allocate new ones into the freed memory
    for (int k = 0; k < 40; k += 2) mats[k].reset();
    CHECK(neighbor_links_consistent());
    for (int k = 0; k < 40; k += 2)
    {
        mats[k].reset(new RowMajorMatrix<double>(8, 8));
        *mats[k] = -double(k);
    }
    for (int k = 0; k < 40; ++k)
    {
        if (k % 2) CHECK_MATRIX_EQ(*mats[k], Ref(16 + k, 16, double(k)));
        else       CHECK_MATRIX_EQ(*mats[k], Ref(8, 8, -double(k)));
    }
    CHECK(neighbor_links_consistent());
}

TEST_CASE("MemoryManager: expressions evaluate correctly with scratch-backed matrices")
{
    matl_test::ScratchMode scratch;
    RowMajorMatrix<double> A(6, 6), B(6, 6), C(6, 6);
    matl_test::fill_random(A);
    matl_test::fill_random(B);
    Ref a = Ref::from(A), b = Ref::from(B);

    C = A * B + A;           CHECK_MATRIX_EQ(C, a * b + a);
    C = A * C;               CHECK_MATRIX_EQ(C, a * (a * b + a));
    C = transpose(A * B);    CHECK_MATRIX_EQ(C, transpose(a * b));
    C = A; C = transpose(C); CHECK_MATRIX_EQ(C, transpose(a));
}

TEST_CASE("MemoryManager: move constructing a scratch-backed matrix transfers the pairing")
{
    matl_test::ScratchMode scratch;
    RowMajorMatrix<double> A(4, 5);
    matl_test::fill_sequence(A);
    Ref expected = Ref::from(A);
    RowMajorMatrix<double> B(std::move(A));
    CHECK_EQ(A.size(), 0ul);
    CHECK_MATRIX_EQ(B, expected);

    // B must still be the borrower: freeing it and reallocating must not disturb other matrices
    RowMajorMatrix<double> C(4, 5);
    C = 1.0;
    CHECK_MATRIX_EQ(B, expected);
}

TEST_CASE("MemoryManager: delete_scratch_memory empties borrowers and can be repeated")
{
    matl_test::ScratchMode scratch;
    MemoryManager& mm = MemoryManager::get_memory_manager();
    {
        RowMajorMatrix<double> A(8, 8);
        A = 1.0;
        mm.delete_scratch_memory();
        CHECK_EQ(A.size(), 0ul);              // borrowed memory was withdrawn
    }
    mm.delete_scratch_memory();
    RowMajorMatrix<double> B(8, 8);           // allocation still works afterwards
    B = 4.0;
    CHECK_MATRIX_EQ(B, Ref(8, 8, 4.0));
}

TEST_CASE("MemoryManager: externally supplied scratch memory is used and never freed")
{
    MemoryManager& mm = MemoryManager::get_memory_manager();
    mm.delete_scratch_memory();
    std::vector<double> buffer(1024, -1.0);
    mm.add_scratch_memory(buffer.data(), buffer.size());

    matl_test::ScratchMode scratch;
    RowMajorMatrix<double> A(4, 4), B(4, 4), C(4, 4);
    A = 1.0; B = 2.0; C = 3.0;

    // All three matrices live inside the supplied buffer
    const double* begin = buffer.data();
    const double* end   = buffer.data() + buffer.size();
    CHECK(&A(1) >= begin && &A(16) < end);
    CHECK(&B(1) >= begin && &B(16) < end);
    CHECK(&C(1) >= begin && &C(16) < end);
    CHECK(neighbor_links_consistent());

    // Removing the buffer withdraws it from every borrower and does not delete it
    mm.remove_scratch_memory(buffer.data());
    CHECK_EQ(A.size(), 0ul);
    CHECK_EQ(B.size(), 0ul);
    CHECK_EQ(C.size(), 0ul);
    buffer[0] = 0.0;                          // still valid memory (ASan would flag a free)
}

TEST_CASE("MemoryManager: restore_defaults keeps containers that are in use")
{
    matl_test::ScratchMode scratch;
    std::vector<std::unique_ptr<RowMajorMatrix<double>>> mats;
    for (int k = 0; k < 20; ++k) { mats.emplace_back(new RowMajorMatrix<double>(64, 64)); *mats.back() = double(k); }
    MemoryManager::get_memory_manager().restore_defaults();
    for (int k = 0; k < 20; ++k) CHECK_MATRIX_EQ(*mats[k], Ref(64, 64, double(k)));
}

TEST_CASE("MemoryManager: fixed-size allocate rejects mismatched dimensions")
{
    FixedSizeRowMajorMatrix<double, 2, 2> F;
    CHECK_THROWS_AS(MemoryManager::allocate(F, 3, 3), std::runtime_error);
    MemoryManager::allocate(F, 2, 2);         // matching dimensions are accepted
}
