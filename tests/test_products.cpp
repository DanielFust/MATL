/**************************************************************************************************
 * File: test_products.cpp
 *
 * Description: The cache-blocked matrix product kernel (see _matl_product_kernel_.hpp), compared
 *              with the reference implementation at sizes around its block boundaries: partial
 *              register tiles (MR x NR), several row blocks (MC), several depth blocks (KC, where
 *              partial sums are reloaded from the target) and several column blocks (NC). Data are
 *              small integers, so every result is exact.
 **************************************************************************************************/
#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    using blocking = detail::product_blocking<double>;

    template<typename T>
    RowMajorMatrix<T> random_matrix(index_t m, index_t n)
    {
        RowMajorMatrix<T> A(m, n);
        matl_test::fill_random(A);
        return A;
    }

    // C = A*B for an M x K and a K x N random matrix, compared with the reference
    template<typename T>
    void check_product(index_t M, index_t K, index_t N)
    {
        RowMajorMatrix<T> A = random_matrix<T>(M, K), B = random_matrix<T>(K, N), C(M, N);
        C = A * B;
        CHECK_MATRIX_EQ(C, Ref::from(A) * Ref::from(B));
    }
}

TEST_CASE("Products: the blocked kernel is used above the size threshold")
{
    static_assert(blocking::MIN_BLOCKED_WORK <= 16 * 16 * 16, "16x16 products use the blocked kernel");
    static_assert(blocking::NR % blocking::LANES == 0);
    CHECK(true);
}

TEST_CASE("Products: sizes that are not multiples of the register tile")
{
    for (index_t M : {9, 13, 17})
        for (index_t K : {10, 15})
            for (index_t N : {11, 18, 23})
                check_product<double>(M, K, N);
    check_product<double>(1, 1000, 1);     // single row and column
    check_product<double>(50, 1, 40);      // rank-one product (K = 1)
}

TEST_CASE("Products: several row and depth blocks (partial sums reloaded from the target)")
{
    check_product<double>(blocking::MC + 7, blocking::KC + 45, 21);
    check_product<double>(2 * blocking::MC + 3, 2 * blocking::KC + 1, 9);
}

TEST_CASE("Products: several column blocks")
{
    check_product<double>(5, 7, blocking::NC + 13);
    check_product<double>(3, blocking::KC + 2, blocking::NC + 1);
}

TEST_CASE("Products: element types")
{
    check_product<float>(37, 41, 43);               // SIMD path with more lanes per vector
    check_product<int>(33, 29, 31);                 // scalar path
    check_product<long>(20, 300, 12);

    // Mixed element types are promoted while packing
    RowMajorMatrix<int>    I = random_matrix<int>(24, 30);
    RowMajorMatrix<double> D = random_matrix<double>(30, 18), C(24, 18);
    C = I * D;
    CHECK_MATRIX_EQ(C, Ref::from(I) * Ref::from(D));
}

TEST_CASE("Products: a target of another element type uses the simple loop")
{
    RowMajorMatrix<double> A = random_matrix<double>(20, 20), B = random_matrix<double>(20, 20);
    RowMajorMatrix<float> C(20, 20);
    C = A * B;
    CHECK_MATRIX_EQ(C, Ref::from(A) * Ref::from(B));
}

TEST_CASE("Products: expression operands, nested products and aliasing at blocked sizes")
{
    const index_t n = 70;
    RowMajorMatrix<double> A = random_matrix<double>(n, n), B = random_matrix<double>(n, n);
    RowMajorMatrix<double> C = random_matrix<double>(n, n), R(n, n);
    Ref a = Ref::from(A), b = Ref::from(B), c = Ref::from(C);

    R = (A + B) * transpose(C);          CHECK_MATRIX_EQ(R, (a + b) * transpose(c));
    R = A * B + transpose(C * A) - C;    CHECK_MATRIX_EQ(R, a * b + transpose(c * a) - c);
    R = A * (B * C);                     CHECK_MATRIX_EQ(R, a * (b * c));
    R = C; Ref r = c;
    R = A * R;      r = a * r;           CHECK_MATRIX_EQ(R, r);     // target is an operand
    R = R * R + R;  r = r * r + r;       CHECK_MATRIX_EQ(R, r);
}

TEST_CASE("Products: fixed-size and mixed fixed-size operands at blocked sizes")
{
    FixedSizeRowMajorMatrix<double, 18, 20> F;
    FixedSizeRowMajorMatrix<double, 20, 14> G;
    FixedSizeRowMajorMatrix<double, 18, 14> H;
    matl_test::fill_random(F);
    matl_test::fill_random(G);
    H = F * G;                                   CHECK_MATRIX_EQ(H, Ref::from(F) * Ref::from(G));

    RowMajorMatrix<double> D = random_matrix<double>(20, 14), R(18, 14);
    R = F * D;                                   CHECK_MATRIX_EQ(R, Ref::from(F) * Ref::from(D));
}
