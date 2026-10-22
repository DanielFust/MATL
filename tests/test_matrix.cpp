/**************************************************************************************************
 * File: test_matrix.cpp
 *
 * Description: Construction, indexing, copying, moving and scalar operations for the MATL
 *              matrix types (RowMajorMatrix and FixedSizeRowMajorMatrix).
 **************************************************************************************************/
#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

/*********************************** RowMajorMatrix ***************************************/

TEST_CASE("RowMajorMatrix: default constructor is empty")
{
    RowMajorMatrix<double> A;
    CHECK_EQ(A.rows(), 0ul);
    CHECK_EQ(A.cols(), 0ul);
    CHECK_EQ(A.size(), 0ul);
}

TEST_CASE("RowMajorMatrix: sized constructor sets dimensions")
{
    RowMajorMatrix<double> A(3, 5);
    CHECK_EQ(A.rows(), 3ul);
    CHECK_EQ(A.cols(), 5ul);
    CHECK_EQ(A.size(), 15ul);
    CHECK(!A.external_mem());
}

TEST_CASE("RowMajorMatrix: 1-based row-major indexing")
{
    RowMajorMatrix<double> A(2, 3);
    matl_test::fill_sequence(A);              // single-index writes 1..6 in storage order
    CHECK_EQ(A(1,1), 1.0); CHECK_EQ(A(1,2), 2.0); CHECK_EQ(A(1,3), 3.0);
    CHECK_EQ(A(2,1), 4.0); CHECK_EQ(A(2,2), 5.0); CHECK_EQ(A(2,3), 6.0);

    A(2,3) = 42.0;
    CHECK_EQ(A(6), 42.0);

    RowMajorMatrix<double> const& C = A;      // const getters
    CHECK_EQ(C(2,3), 42.0);
    CHECK_EQ(C(1), 1.0);
}

TEST_CASE("RowMajorMatrix: preallocated (external) memory is used in place")
{
    double buffer[6] = {1, 2, 3, 4, 5, 6};
    {
        RowMajorMatrix<double> A(2, 3, buffer);
        CHECK(A.external_mem());
        CHECK_EQ(A(2,1), 4.0);
        A(1,1) = 10.0;
        CHECK_EQ(buffer[0], 10.0);            // writes go straight to the caller's buffer
    }
    CHECK_EQ(buffer[5], 6.0);                 // destructor does not free or clobber the buffer
}

TEST_CASE("RowMajorMatrix: copy constructor (same type) makes a deep copy")
{
    RowMajorMatrix<double> A(2, 3);
    matl_test::fill_sequence(A);
    RowMajorMatrix<double> B(A);
    CHECK_MATRIX_EQ(B, Ref::from(A));
    B(1,1) = -1.0;
    CHECK_EQ(A(1,1), 1.0);
}

TEST_CASE("RowMajorMatrix: generic copy constructor handles non-square matrices")
{
    FixedSizeRowMajorMatrix<double, 2, 3> F;
    matl_test::fill_sequence(F);
    RowMajorMatrix<double> A(F);
    CHECK_MATRIX_EQ(A, Ref::from(F));
}

TEST_CASE("RowMajorMatrix: move constructor steals heap data")
{
    RowMajorMatrix<double> A(3, 4);
    matl_test::fill_sequence(A);
    Ref expected = Ref::from(A);
    RowMajorMatrix<double> B(std::move(A));
    CHECK_MATRIX_EQ(B, expected);
    CHECK_EQ(A.size(), 0ul);
}

TEST_CASE("RowMajorMatrix: move assignment")
{
    RowMajorMatrix<double> A(3, 4), B;
    matl_test::fill_sequence(A);
    Ref expected = Ref::from(A);
    B = std::move(A);
    CHECK_MATRIX_EQ(B, expected);
    CHECK_EQ(A.size(), 0ul);
}

TEST_CASE("RowMajorMatrix: copy assignment from matrix")
{
    RowMajorMatrix<double> A(2, 3), B(2, 3), C;
    matl_test::fill_sequence(A);
    B = A;
    CHECK_MATRIX_EQ(B, Ref::from(A));
    C = A;                                    // unallocated target is allocated on assignment
    CHECK_MATRIX_EQ(C, Ref::from(A));
}

TEST_CASE("RowMajorMatrix: copy() from another matrix type, including multi-block sizes")
{
    // The generic copy is cache-blocked in 32x32 tiles, so test sizes that are not multiples
    // of the block size, and sizes spanning several blocks.
    FixedSizeRowMajorMatrix<double, 2, 3> small;
    matl_test::fill_sequence(small);
    RowMajorMatrix<double> A(2, 3);
    A = 0.0;
    A.copy(small);
    CHECK_MATRIX_EQ(A, Ref::from(small));

    FixedSizeRowMajorMatrix<double, 40, 70> large;
    matl_test::fill_sequence(large);
    RowMajorMatrix<double> B(40, 70);
    B = 0.0;
    B.copy(large);
    CHECK_MATRIX_EQ(B, Ref::from(large));
}

TEST_CASE("RowMajorMatrix: scalar assignment and compound scalar operators")
{
    RowMajorMatrix<double> A(2, 2);
    A = 3.0;   CHECK_MATRIX_EQ(A, Ref(2, 2, 3.0));
    A += 2.0;  CHECK_MATRIX_EQ(A, Ref(2, 2, 5.0));
    A -= 1.0;  CHECK_MATRIX_EQ(A, Ref(2, 2, 4.0));
    A *= 2.5;  CHECK_MATRIX_EQ(A, Ref(2, 2, 10.0));
    A /= 4.0;  CHECK_MATRIX_EQ(A, Ref(2, 2, 2.5));
}

TEST_CASE("RowMajorMatrix: integer division assignment divides exactly")
{
    RowMajorMatrix<int> I(2, 2);
    I = 10;
    I /= 3;
    for (index_t k = 1; k <= 4; ++k) CHECK_EQ(I(k), 3);

    RowMajorMatrix<double> D(2, 2);
    D = 10.0;
    D /= 4;                                   // integer divisor on a floating-point matrix
    for (index_t k = 1; k <= 4; ++k) CHECK_EQ(D(k), 2.5);
}

TEST_CASE("RowMajorMatrix: compound assignment by matrix and expression")
{
    RowMajorMatrix<double> A(3, 3), B(3, 3);
    matl_test::fill_sequence(A);
    matl_test::fill_sequence(B, 10.0);
    Ref a = Ref::from(A), b = Ref::from(B);

    A += B;            a = a + b;             CHECK_MATRIX_EQ(A, a);
    A -= B;            a = a - b;             CHECK_MATRIX_EQ(A, a);
    A *= B;            a = a * b;             CHECK_MATRIX_EQ(A, a);
    A += B + B;        a = a + (b + b);       CHECK_MATRIX_EQ(A, a);
    A -= B * 2.0;      a = a - b * 2.0;       CHECK_MATRIX_EQ(A, a);
    A *= B - B * 0.5;  a = a * (b - b * 0.5); CHECK_MATRIX_EQ(A, a);
}

TEST_CASE("RowMajorMatrix: compound assignment with the target on the right-hand side")
{
    RowMajorMatrix<double> A(3, 3);
    matl_test::fill_sequence(A);
    Ref a = Ref::from(A);

    A += A;            a = a + a;             CHECK_MATRIX_EQ(A, a);
    A *= A;            a = a * a;             CHECK_MATRIX_EQ(A, a);
    A -= transpose(A); a = a - transpose(a);  CHECK_MATRIX_EQ(A, a);
}

TEST_CASE("RowMajorMatrix: mixed element types promote")
{
    RowMajorMatrix<int>    I(2, 2);
    RowMajorMatrix<double> D(2, 2), C(2, 2);
    I = 1;
    D = 0.5;
    C = I + D;
    CHECK_MATRIX_EQ(C, Ref(2, 2, 1.5));
}

/******************************* FixedSizeRowMajorMatrix ************************************/

TEST_CASE("FixedSizeRowMajorMatrix: compile-time dimensions")
{
    using F = FixedSizeRowMajorMatrix<double, 3, 4>;
    static_assert(F::ROWS == 3 && F::COLS == 4 && F::SIZE == 12);
    F A;
    CHECK_EQ(A.rows(), 3ul);
    CHECK_EQ(A.cols(), 4ul);
    CHECK_EQ(A.size(), 12ul);
}

TEST_CASE("FixedSizeRowMajorMatrix: 1-based row-major indexing")
{
    FixedSizeRowMajorMatrix<double, 2, 3> A;
    matl_test::fill_sequence(A);
    CHECK_EQ(A(1,3), 3.0);
    CHECK_EQ(A(2,1), 4.0);
    A(2,2) = -5.0;
    CHECK_EQ(A(5), -5.0);
}

TEST_CASE("FixedSizeRowMajorMatrix: scalar and matrix assignment")
{
    FixedSizeRowMajorMatrix<double, 2, 2> A, B;
    A = 7.0;
    CHECK_MATRIX_EQ(A, Ref(2, 2, 7.0));
    B = A;
    CHECK_MATRIX_EQ(B, Ref(2, 2, 7.0));
}

TEST_CASE("FixedSizeRowMajorMatrix: assignment from a RowMajorMatrix")
{
    RowMajorMatrix<double> D(2, 3);
    matl_test::fill_sequence(D);
    FixedSizeRowMajorMatrix<double, 2, 3> F;
    F = D;
    CHECK_MATRIX_EQ(F, Ref::from(D));
}

TEST_CASE("Matrix alias selects dynamic or fixed-size storage")
{
    static_assert(std::is_same_v<Matrix<double>,       RowMajorMatrix<double>>);
    static_assert(std::is_same_v<Matrix<float, 2, 3>,  FixedSizeRowMajorMatrix<float, 2, 3>>);
}

/************************************* API functions ****************************************/

TEST_CASE("eye: dynamic identity matrices")
{
    RowMajorMatrix<double> I = eye<double>(3);
    Ref expected(3, 3);
    for (index_t k = 1; k <= 3; ++k) expected(k, k) = 1.0;
    CHECK_MATRIX_EQ(I, expected);

    RowMajorMatrix<double> R = eye<double>(2, 3);
    Ref rect(2, 3);
    rect(1,1) = rect(2,2) = 1.0;
    CHECK_MATRIX_EQ(R, rect);
}

TEST_CASE("eye: fixed-size identity matrices")
{
    auto I = eye<double, 3>();
    static_assert(std::is_same_v<decltype(I), FixedSizeRowMajorMatrix<double, 3, 3>>);
    Ref expected(3, 3);
    for (index_t k = 1; k <= 3; ++k) expected(k, k) = 1.0;
    CHECK_MATRIX_EQ(I, expected);

    auto R = eye<double, 2, 3>();
    static_assert(std::is_same_v<decltype(R), FixedSizeRowMajorMatrix<double, 2, 3>>);
    Ref rect(2, 3);
    rect(1,1) = rect(2,2) = 1.0;
    CHECK_MATRIX_EQ(R, rect);
}
