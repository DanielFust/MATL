/**************************************************************************************************
 * File: test_randomized.cpp
 *
 * Description: Randomized differential tests. Random matrix sizes and small-integer contents are
 *              generated from a fixed seed, MATL results are compared with the naive reference
 *              implementation in matl_test.hpp. Integer-valued data keeps results exact.
 **************************************************************************************************/
#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    const int TRIALS = 50;

    RowMajorMatrix<double> random_matrix(index_t m, index_t n)
    {
        RowMajorMatrix<double> A(m, n);
        matl_test::fill_random(A);
        return A;
    }

    void run_randomized_trials()
    {
        for (int trial = 0; trial < TRIALS; ++trial)
        {
            const index_t m = matl_test::random_dim(), n = matl_test::random_dim(), p = matl_test::random_dim();
            RowMajorMatrix<double> A = random_matrix(m, n), B = random_matrix(m, n);
            RowMajorMatrix<double> C = random_matrix(n, p), D = random_matrix(p, n);
            Ref a = Ref::from(A), b = Ref::from(B), c = Ref::from(C), d = Ref::from(D);

            RowMajorMatrix<double> R;
            R = A + B;                       CHECK_MATRIX_EQ(R, a + b);
            R = A - B * 2.0 + 3.0;           CHECK_MATRIX_EQ(R, a - b * 2.0 + 3.0);

            RowMajorMatrix<double> P;
            P = A * C;                       CHECK_MATRIX_EQ(P, a * c);
            P = (A + B) * C;                 CHECK_MATRIX_EQ(P, (a + b) * c);
            P = A * C + B * C;               CHECK_MATRIX_EQ(P, a * c + b * c);
            P = A * transpose(D);            CHECK_MATRIX_EQ(P, a * transpose(d));

            RowMajorMatrix<double> T;
            T = transpose(A);                CHECK_MATRIX_EQ(T, transpose(a));
            T = transpose(A - B);            CHECK_MATRIX_EQ(T, transpose(a - b));
            RowMajorMatrix<double> U;
            U = transpose(A * C);            CHECK_MATRIX_EQ(U, transpose(c) * transpose(a));
            U = transpose(C) * transpose(A); CHECK_MATRIX_EQ(U, transpose(a * c));

            RowMajorMatrix<double> Q;
            Q = A * C * D;                   CHECK_MATRIX_EQ(Q, a * c * d);

            // Square aliasing cases
            RowMajorMatrix<double> S = random_matrix(n, n);
            Ref s = Ref::from(S);
            S = S * S + S;                   s = s * s + s;            CHECK_MATRIX_EQ(S, s);
            S = transpose(S) - S * 2.0;      s = transpose(s) - s * 2.0; CHECK_MATRIX_EQ(S, s);
        }
    }
}

TEST_CASE("Randomized: expressions match the reference implementation (heap memory)")
{
    run_randomized_trials();
}

TEST_CASE("Randomized: expressions match the reference implementation (scratch memory)")
{
    matl_test::ScratchMode scratch;
    run_randomized_trials();
}
