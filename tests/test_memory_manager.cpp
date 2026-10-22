/**************************************************************************************************
 * File: test_memory_manager.cpp
 *
 * Description: Matrix memory ownership (heap memory owned by the matrix, or externally owned
 *              memory that the matrix uses but never frees), allocate/deallocate, and the
 *              per-thread MemoryManager settings.
 **************************************************************************************************/
#include <thread>
#include <type_traits>
#include <vector>

#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

/********************************** MemoryManager settings **********************************/

TEST_CASE("MemoryManager: is enabled by default and caches temporaries")
{
    CHECK(MemoryManager::get_memory_manager().enabled());
    CHECK(MemoryManager::CACHE_TEMPORARY_MATRICES);
    CHECK(MemoryManager::use_temporary_arena());
}

TEST_CASE("MemoryManager: disabling it allocates temporaries on the heap; restore_defaults resets")
{
    MemoryManager& mm = MemoryManager::get_memory_manager();
    mm.disable();
    CHECK(!MemoryManager::use_temporary_arena());
    mm.restore_defaults();
    CHECK(mm.enabled());

    MemoryManager::CACHE_TEMPORARY_MATRICES = false;
    CHECK(!MemoryManager::use_temporary_arena());
    mm.restore_defaults();
    CHECK(MemoryManager::CACHE_TEMPORARY_MATRICES);
    CHECK(MemoryManager::use_temporary_arena());
}

TEST_CASE("MemoryManager: each thread has its own memory manager and arena")
{
    RowMajorMatrix<double> A(4, 4), B(4, 4), R(4, 4);
    matl_test::fill_sequence(A);
    matl_test::fill_random(B);
    Ref expected = Ref::from(A) * Ref::from(B) + Ref::from(A);

    TemporaryArena* main_arena  = &MemoryManager::temporary_arena();
    TemporaryArena* other_arena = nullptr;
    bool            other_ok    = false;
    std::thread worker([&]
    {
        other_arena = &MemoryManager::temporary_arena();
        RowMajorMatrix<double> W(4, 4);
        W = A * B + A;                        // uses the worker thread's arena
        other_ok = (MemoryManager::temporary_arena().high_water_mark() > 0);
        for (index_t i = 1; i <= W.size(); ++i) other_ok = other_ok && (W(i) == expected.data[i-1]);
    });                                       // the worker's arena is freed when it exits
    worker.join();

    CHECK(other_arena != main_arena);
    CHECK(other_ok);
    R = A * B + A;
    CHECK_MATRIX_EQ(R, expected);
}

/********************************* Matrices own heap memory *********************************/

TEST_CASE("Memory: a matrix is a plain value type (no virtual functions)")
{
    static_assert(!std::is_polymorphic_v<RowMajorMatrix<double>>, "RowMajorMatrix has no vtable");
    static_assert(std::is_nothrow_move_constructible_v<RowMajorMatrix<double>>);
    static_assert(std::is_nothrow_move_assignable_v<RowMajorMatrix<double>>);
    static_assert(std::is_copy_assignable_v<RowMajorMatrix<double>>);
    CHECK(true);
}

TEST_CASE("Memory: copies are deep and moves transfer the data")
{
    RowMajorMatrix<double> A(3, 4);
    matl_test::fill_sequence(A);
    Ref a = Ref::from(A);

    RowMajorMatrix<double> B(A);              // copy constructor
    B(1,1) = -1.0;
    CHECK_MATRIX_EQ(A, a);

    const double* data = &A(1);
    RowMajorMatrix<double> C(std::move(A));   // move constructor takes the same memory
    CHECK_EQ(&C(1), data);
    CHECK_EQ(A.size(), 0ul);
    CHECK_MATRIX_EQ(C, a);

    RowMajorMatrix<double> D;
    D = std::move(C);                         // move assignment takes the same memory
    CHECK_EQ(&D(1), data);
    CHECK_EQ(C.size(), 0ul);
    CHECK(!D.external_mem());
    CHECK_MATRIX_EQ(D, a);
}

TEST_CASE("Memory: copy assignment from a const matrix")
{
    RowMajorMatrix<double> A(2, 3);
    matl_test::fill_sequence(A);
    RowMajorMatrix<double> const& c = A;

    RowMajorMatrix<double> B;                 // empty: allocated to fit
    B = c;
    CHECK_MATRIX_EQ(B, Ref::from(A));
    CHECK(&B(1) != &A(1));

    RowMajorMatrix<double> C(2, 3);           // allocated: copied into place
    const double* data = &C(1);
    C = c;
    CHECK_EQ(&C(1), data);
    CHECK_MATRIX_EQ(C, Ref::from(A));
}

TEST_CASE("Memory: allocate and deallocate")
{
    RowMajorMatrix<double> A;
    allocate(A, 3, 2);
    CHECK_EQ(A.rows(), 3ul);
    CHECK_EQ(A.cols(), 2ul);
    CHECK(!A.external_mem());
    A = 7.0;
    CHECK_MATRIX_EQ(A, Ref(3, 2, 7.0));

    allocate(A, 4, 4);                        // reallocating releases the old memory
    A = 1.0;
    CHECK_MATRIX_EQ(A, Ref(4, 4, 1.0));

    deallocate(A);
    CHECK_EQ(A.size(), 0ul);
    deallocate(A);                            // deallocating an empty matrix is harmless
    CHECK_EQ(A.size(), 0ul);
}

TEST_CASE("Memory: fixed-size allocate rejects mismatched dimensions")
{
    FixedSizeRowMajorMatrix<double, 2, 2> F;
    CHECK_THROWS_AS(allocate(F, 3, 3), std::runtime_error);
    allocate(F, 2, 2);                        // matching dimensions are accepted
}

/******************************** Externally owned memory **********************************/

TEST_CASE("Memory: many matrices can share one externally owned block")
{
    std::vector<double> block(3 * 16, -1.0);
    {
        RowMajorMatrix<double> A(4, 4, block.data()), B(4, 4, block.data() + 16), C(4, 4, block.data() + 32);
        A = 1.0; B = 2.0;
        C = A + B * 2.0;
        CHECK_MATRIX_EQ(C, Ref(4, 4, 5.0));
        CHECK(A.external_mem() && B.external_mem() && C.external_mem());
    }
    CHECK_EQ(block[0], 1.0);                  // the data stayed in the block
    CHECK_EQ(block[16], 2.0);
    CHECK_EQ(block[47], 5.0);                 // and the block was not freed (ASan would flag it)
}

TEST_CASE("Memory: moving a matrix over external memory keeps it external")
{
    double buffer[4] = {1, 2, 3, 4};
    RowMajorMatrix<double> A(2, 2, buffer);
    RowMajorMatrix<double> B(std::move(A));
    CHECK(B.external_mem());
    CHECK_EQ(&B(1), &buffer[0]);
    CHECK_EQ(A.size(), 0ul);
    CHECK(!A.external_mem());

    RowMajorMatrix<double> C;
    C = std::move(B);
    CHECK(C.external_mem());
    CHECK_EQ(&C(1), &buffer[0]);
    C(2,2) = 40.0;
    CHECK_EQ(buffer[3], 40.0);
}

TEST_CASE("Memory: assigning into a matrix over external memory writes into the buffer")
{
    double buffer[4] = {0, 0, 0, 0};
    RowMajorMatrix<double> E(2, 2, buffer);

    RowMajorMatrix<double> A(2, 2);
    A = 3.0;
    E = std::move(A);                         // data is copied: E keeps the caller's buffer
    CHECK(E.external_mem());
    CHECK_EQ(&E(1), &buffer[0]);
    CHECK_EQ(buffer[2], 3.0);

    RowMajorMatrix<double> B(2, 2);
    B = 5.0;
    E = B;
    CHECK_EQ(buffer[3], 5.0);
    E = B * 2.0;
    CHECK_EQ(buffer[0], 10.0);
}

TEST_CASE("Memory: a copy of a matrix over external memory owns its data")
{
    double buffer[4] = {1, 2, 3, 4};
    RowMajorMatrix<double> E(2, 2, buffer);
    RowMajorMatrix<double> C(E);
    CHECK(!C.external_mem());
    CHECK(&C(1) != &buffer[0]);
    C(1,1) = 100.0;
    CHECK_EQ(buffer[0], 1.0);

    deallocate(E);                            // stops using the buffer without freeing it
    CHECK_EQ(E.size(), 0ul);
    CHECK(!E.external_mem());
    CHECK_EQ(buffer[0], 1.0);
}
