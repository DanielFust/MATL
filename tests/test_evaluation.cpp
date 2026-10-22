/**************************************************************************************************
 * File: test_evaluation.cpp
 *
 * Description: Deferred evaluation. Expressions are pure formulas: building one computes and
 *              allocates nothing. All work, including matrix products nested inside larger
 *              expressions, happens when the expression is assigned to a matrix. Temporaries
 *              needed during an assignment are released when the assignment finishes.
 **************************************************************************************************/
#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    // Memory currently held by evaluation temporaries
    index_t temporary_bytes_in_use()
    {return MemoryManager::temporary_arena().bytes_in_use();}

    struct Operands
    {
        RowMajorMatrix<double> A{3, 3}, B{3, 3}, C{3, 3}, D{3, 3};
        Operands()
        {
            matl_test::fill_sequence(A, 1.0);
            matl_test::fill_sequence(B, -4.0);
            matl_test::fill_sequence(C, 2.0);
            matl_test::fill_random(D);
            C(1,3) = 9.0;
        }
    };
}

/****************************** Products are evaluated lazily *******************************/

TEST_CASE("Evaluation: a nested product reads the current values of its operands")
{
    Operands s;
    auto e = s.A * s.B + s.C;                 // nothing is computed here
    RowMajorMatrix<double> R(3, 3);

    R = e;  CHECK_MATRIX_EQ(R, Ref::from(s.A) * Ref::from(s.B) + Ref::from(s.C));
    s.A = 2.0;                                // change an operand of the nested product
    R = e;  CHECK_MATRIX_EQ(R, Ref::from(s.A) * Ref::from(s.B) + Ref::from(s.C));
    s.B(2,2) = 100.0;
    R = e;  CHECK_MATRIX_EQ(R, Ref::from(s.A) * Ref::from(s.B) + Ref::from(s.C));
}

TEST_CASE("Evaluation: products of expressions are evaluated lazily")
{
    Operands s;
    auto e = (s.A + s.B) * transpose(s.C) - 2.0 * (s.A * s.C);
    s.A = 1.0;
    s.C(1,1) = -3.0;
    RowMajorMatrix<double> R(3, 3);
    Ref a = Ref::from(s.A), b = Ref::from(s.B), c = Ref::from(s.C);
    R = e;  CHECK_MATRIX_EQ(R, (a + b) * transpose(c) - 2.0 * (a * c));
}

TEST_CASE("Evaluation: an expression gives the same result each time it is evaluated")
{
    Operands s;
    auto e = transpose(s.A * s.B) + s.C * s.D;
    RowMajorMatrix<double> R1(3, 3), R2(3, 3);
    R1 = e;
    R2 = e;
    CHECK_MATRIX_EQ(R1, Ref::from(R2));
    CHECK_MATRIX_EQ(R1, transpose(Ref::from(s.A) * Ref::from(s.B)) + Ref::from(s.C) * Ref::from(s.D));
}

TEST_CASE("Evaluation: products are stored as formulas, not results")
{
    using M    = RowMajorMatrix<double>;
    using Prod = expr::Product<M, M>;
    static_assert(std::is_same_v<detail::stored_t<Prod>, Prod>,
                  "a product nested in an expression is stored as a product expression");
    CHECK(true);
}

TEST_CASE("Evaluation: deeply nested products")
{
    Operands s;
    Ref a = Ref::from(s.A), b = Ref::from(s.B), c = Ref::from(s.C), d = Ref::from(s.D);
    RowMajorMatrix<double> R(3, 3);
    R = (s.A * s.B) * (s.C * s.D) + s.A * s.D;    CHECK_MATRIX_EQ(R, (a * b) * (c * d) + a * d);
    R = s.A * (s.B * (s.C * s.D));                CHECK_MATRIX_EQ(R, a * (b * (c * d)));
    R = transpose(transpose(s.A * s.B) * s.C);    CHECK_MATRIX_EQ(R, transpose(transpose(a * b) * c));
    R = ((s.A + s.B) * 2.0) * (s.C - s.D) / 4.0;  CHECK_MATRIX_EQ(R, ((a + b) * 2.0) * (c - d) / 4.0);
}

TEST_CASE("Evaluation: fixed-size nested products are evaluated lazily")
{
    FixedSizeRowMajorMatrix<double, 2, 3> A;
    FixedSizeRowMajorMatrix<double, 3, 2> B;
    FixedSizeRowMajorMatrix<double, 2, 2> C, R;
    matl_test::fill_sequence(A);
    matl_test::fill_sequence(B, -2.0);
    C = 1.0;
    auto e = A * B + C;
    A = 2.0;
    R = e;
    CHECK_MATRIX_EQ(R, Ref::from(A) * Ref::from(B) + Ref::from(C));
}

/************************* Aliasing with deferred products ***********************************/

TEST_CASE("Evaluation: target used as an operand of products at the top level")
{
    Operands s;
    Ref a = Ref::from(s.A), b = Ref::from(s.B), c = Ref::from(s.C), r;
    RowMajorMatrix<double> R(3, 3);
    R = s.C; r = c;  R = s.A * R;               r = a * r;               CHECK_MATRIX_EQ(R, r);
    R = s.C; r = c;  R = R * R;                 r = r * r;               CHECK_MATRIX_EQ(R, r);
    R = s.C; r = c;  R = (s.A * R) * s.B;       r = (a * r) * b;         CHECK_MATRIX_EQ(R, r);
    R = s.C; r = c;  R = s.A * (s.B * R);       r = a * (b * r);         CHECK_MATRIX_EQ(R, r);
    R = s.C; r = c;  R = transpose(R) * R;      r = transpose(r) * r;    CHECK_MATRIX_EQ(R, r);
    R = s.C; r = c;  R = (R + s.A) * R;         r = (r + a) * r;         CHECK_MATRIX_EQ(R, r);
}

TEST_CASE("Evaluation: target used inside products nested in elementwise expressions")
{
    Operands s;
    Ref a = Ref::from(s.A), b = Ref::from(s.B), c = Ref::from(s.C), r;
    RowMajorMatrix<double> R(3, 3);
    R = s.C; r = c;  R = s.A * R + R;               r = a * r + r;               CHECK_MATRIX_EQ(R, r);
    R = s.C; r = c;  R = R - R * s.B;               r = r - r * b;               CHECK_MATRIX_EQ(R, r);
    R = s.C; r = c;  R = transpose(R * s.A);        r = transpose(r * a);        CHECK_MATRIX_EQ(R, r);
    R = s.C; r = c;  R = transpose(R * s.A) + R;    r = transpose(r * a) + r;    CHECK_MATRIX_EQ(R, r);
    R = s.C; r = c;  R = 2.0 * (R * R) - transpose(R);
                     r = 2.0 * (r * r) - transpose(r);                           CHECK_MATRIX_EQ(R, r);
}

/************************ No work or memory until assignment *********************************/

TEST_CASE("Evaluation: building an expression allocates no memory")
{
    Operands s;
    MemoryManager::temporary_arena().reset_high_water_mark();
    auto e = s.A * s.B + transpose(s.C * s.A) - (s.A + s.D) * s.B;
    CHECK_EQ(MemoryManager::temporary_arena().high_water_mark(), 0ul);   // no temporaries

    RowMajorMatrix<double> R(3, 3);
    R = e;
    Ref a = Ref::from(s.A), b = Ref::from(s.B), c = Ref::from(s.C), d = Ref::from(s.D);
    CHECK_MATRIX_EQ(R, a * b + transpose(c * a) - (a + d) * b);
}

TEST_CASE("Evaluation: temporaries are released when an assignment finishes")
{
    Operands s;
    RowMajorMatrix<double> R(3, 3);
    MemoryManager::temporary_arena().reset_high_water_mark();

    R = s.A * s.B + s.C * s.D;                    CHECK_EQ(temporary_bytes_in_use(), 0ul);
    R = (s.A * s.B) * (s.C * s.D) - s.A;          CHECK_EQ(temporary_bytes_in_use(), 0ul);
    R = transpose(s.A * R) + (s.B - s.C) * R;     CHECK_EQ(temporary_bytes_in_use(), 0ul);
    R = R * R;                                    CHECK_EQ(temporary_bytes_in_use(), 0ul);
    CHECK(MemoryManager::temporary_arena().high_water_mark() > 0);   // temporaries were used
}

TEST_CASE("Evaluation: initializing a matrix from a product expression")
{
    Operands s;
    Matrix<double> P = s.A * s.B + s.C;
    CHECK_MATRIX_EQ(P, Ref::from(s.A) * Ref::from(s.B) + Ref::from(s.C));
    Matrix<double, 3, 3> F = s.A * transpose(s.B);
    CHECK_MATRIX_EQ(F, Ref::from(s.A) * transpose(Ref::from(s.B)));
}
