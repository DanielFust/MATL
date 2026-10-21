/**************************************************************************************************
 * File: test_expressions.cpp
 *
 * Description: Lazy expression templates: elementwise arithmetic, scalar arithmetic, matrix
 *              multiplication, transposition, nested expressions and aliasing (the target matrix
 *              appearing on the right-hand side).
 **************************************************************************************************/
#include "matl_test.hpp"

using namespace matl;
using matl_test::Ref;

namespace
{
    // Three 3x3 matrices with distinct, non-symmetric contents and their references
    struct Square3
    {
        RowMajorMatrix<double> A{3, 3}, B{3, 3}, C{3, 3};
        Ref a, b, c;
        Square3()
        {
            matl_test::fill_sequence(A, 1.0);
            matl_test::fill_sequence(B, -4.0);
            matl_test::fill_sequence(C, 2.0);
            C(1,3) = 9.0;                     // break the arithmetic-sequence structure
            a = Ref::from(A); b = Ref::from(B); c = Ref::from(C);
        }
    };
}

/************************************ Elementwise ********************************************/

TEST_CASE("Expressions: matrix + matrix and matrix - matrix")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = s.A + s.B;           CHECK_MATRIX_EQ(R, s.a + s.b);
    R = s.A - s.B;           CHECK_MATRIX_EQ(R, s.a - s.b);
    R = s.A + s.B - s.C;     CHECK_MATRIX_EQ(R, s.a + s.b - s.c);
    R = s.A - (s.B - s.C);   CHECK_MATRIX_EQ(R, s.a - (s.b - s.c));
    R = (s.A + s.B) + (s.B + s.C); CHECK_MATRIX_EQ(R, (s.a + s.b) + (s.b + s.c));
}

TEST_CASE("Expressions: scalar arithmetic on matrices")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = s.A + 2.0;           CHECK_MATRIX_EQ(R, s.a + 2.0);
    R = 2.0 + s.A;           CHECK_MATRIX_EQ(R, 2.0 + s.a);
    R = s.A - 2.0;           CHECK_MATRIX_EQ(R, s.a - 2.0);
    R = 2.0 - s.A;           CHECK_MATRIX_EQ(R, 2.0 - s.a);
    R = s.A * 3.0;           CHECK_MATRIX_EQ(R, s.a * 3.0);
    R = 3.0 * s.A;           CHECK_MATRIX_EQ(R, 3.0 * s.a);
    R = s.A / 4.0;           CHECK_MATRIX_EQ(R, s.a / 4.0);
    R = s.A / 4;             CHECK_MATRIX_EQ(R, s.a / 4.0);
}

TEST_CASE("Expressions: scalar arithmetic on expressions")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = (s.A + s.B) + 2.0;   CHECK_MATRIX_EQ(R, (s.a + s.b) + 2.0);
    R = 2.0 + (s.A + s.B);   CHECK_MATRIX_EQ(R, 2.0 + (s.a + s.b));
    R = (s.A + s.B) - 2.0;   CHECK_MATRIX_EQ(R, (s.a + s.b) - 2.0);
    R = 2.0 - (s.A + s.B);   CHECK_MATRIX_EQ(R, 2.0 - (s.a + s.b));
    R = (s.A + s.B) * 3.0;   CHECK_MATRIX_EQ(R, (s.a + s.b) * 3.0);
    R = 3.0 * (s.A + s.B);   CHECK_MATRIX_EQ(R, 3.0 * (s.a + s.b));
    R = (s.A + s.B) / 4.0;   CHECK_MATRIX_EQ(R, (s.a + s.b) / 4.0);
    R = (s.A + s.B) / 4;     CHECK_MATRIX_EQ(R, (s.a + s.b) / 4.0);
}

TEST_CASE("Expressions: division by a scalar of a different floating-point type")
{
    RowMajorMatrix<double> A(2, 2), R(2, 2);
    A = 3.0;
    R = A / 2.0f;
    CHECK_MATRIX_EQ(R, Ref(2, 2, 1.5));
    R = (A + A) / 4.0f;
    CHECK_MATRIX_EQ(R, Ref(2, 2, 1.5));
}

TEST_CASE("Expressions: integer matrices divided by integers use integer division")
{
    RowMajorMatrix<int> I(2, 2), R(2, 2);
    I = 7;
    R = I / 2;
    for (index_t k = 1; k <= 4; ++k) CHECK_EQ(R(k), 3);
}

TEST_CASE("Expressions: long mixed elementwise chain")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = 1.0 + s.A * 2.0 + s.B - 3.0 - s.B / 2.0 + 4.0 * s.C;
    CHECK_MATRIX_EQ(R, 1.0 + s.a * 2.0 + s.b - 3.0 - s.b / 2.0 + 4.0 * s.c);
}

TEST_CASE("Expressions: assignment allocates an empty target")
{
    Square3 s;
    RowMajorMatrix<double> R;
    R = s.A + s.B;
    CHECK_MATRIX_EQ(R, s.a + s.b);
}

/******************************** Matrix multiplication ***************************************/

TEST_CASE("Expressions: square matrix product")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = s.A * s.B;           CHECK_MATRIX_EQ(R, s.a * s.b);
    R = s.B * s.A;           CHECK_MATRIX_EQ(R, s.b * s.a);
}

TEST_CASE("Expressions: non-square matrix product")
{
    RowMajorMatrix<double> A(2, 3), B(3, 4), R(2, 4);
    matl_test::fill_sequence(A);
    matl_test::fill_sequence(B, -3.0);
    R = A * B;
    CHECK_MATRIX_EQ(R, Ref::from(A) * Ref::from(B));

    RowMajorMatrix<double> col(3, 1), row(1, 3), outer(3, 3), inner(1, 1);
    matl_test::fill_sequence(col);
    matl_test::fill_sequence(row, 4.0);
    outer = col * row;       CHECK_MATRIX_EQ(outer, Ref::from(col) * Ref::from(row));
    inner = row * col;       CHECK_MATRIX_EQ(inner, Ref::from(row) * Ref::from(col));
}

TEST_CASE("Expressions: products combined with elementwise operations")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = s.A + s.B * s.C;           CHECK_MATRIX_EQ(R, s.a + s.b * s.c);
    R = s.A * s.B + s.C;           CHECK_MATRIX_EQ(R, s.a * s.b + s.c);
    R = s.A * s.B - s.C;           CHECK_MATRIX_EQ(R, s.a * s.b - s.c);
    R = s.C - s.A * s.B;           CHECK_MATRIX_EQ(R, s.c - s.a * s.b);
    R = s.A * s.B + s.B * s.C;     CHECK_MATRIX_EQ(R, s.a * s.b + s.b * s.c);
    R = s.A * s.B - s.B * s.C;     CHECK_MATRIX_EQ(R, s.a * s.b - s.b * s.c);
    R = (s.A * s.B) * 2.0;         CHECK_MATRIX_EQ(R, (s.a * s.b) * 2.0);
    R = 2.0 * (s.A * s.B);         CHECK_MATRIX_EQ(R, 2.0 * (s.a * s.b));
    R = (s.A * s.B) + 1.0;         CHECK_MATRIX_EQ(R, (s.a * s.b) + 1.0);
    R = (s.A * s.B) - 1.0;         CHECK_MATRIX_EQ(R, (s.a * s.b) - 1.0);
    R = 1.0 - (s.A * s.B);         CHECK_MATRIX_EQ(R, 1.0 - (s.a * s.b));
    R = (s.A * s.B) / 2.0;         CHECK_MATRIX_EQ(R, (s.a * s.b) / 2.0);
    R = (s.A * s.B) / 2;           CHECK_MATRIX_EQ(R, (s.a * s.b) / 2.0);
}

TEST_CASE("Expressions: products of expressions")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = s.A * (s.B + s.C);         CHECK_MATRIX_EQ(R, s.a * (s.b + s.c));
    R = (s.A + s.B) * s.C;         CHECK_MATRIX_EQ(R, (s.a + s.b) * s.c);
    R = (s.A + s.B) * (s.A - s.C); CHECK_MATRIX_EQ(R, (s.a + s.b) * (s.a - s.c));
    R = s.A * s.B * s.C;           CHECK_MATRIX_EQ(R, s.a * s.b * s.c);
    R = (s.A * s.B) * (s.B * s.C); CHECK_MATRIX_EQ(R, (s.a * s.b) * (s.b * s.c));
}

/************************************** Transpose *********************************************/

TEST_CASE("Transpose: square and non-square matrices")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = transpose(s.A);            CHECK_MATRIX_EQ(R, transpose(s.a));

    RowMajorMatrix<double> W(2, 3), T(3, 2);
    matl_test::fill_sequence(W);
    T = transpose(W);              CHECK_MATRIX_EQ(T, transpose(Ref::from(W)));
    W = transpose(T);              CHECK_MATRIX_EQ(W, transpose(Ref::from(T)));
}

TEST_CASE("Transpose: mixed with untransposed operands")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = s.A + transpose(s.B);      CHECK_MATRIX_EQ(R, s.a + transpose(s.b));
    R = s.A - transpose(s.B);      CHECK_MATRIX_EQ(R, s.a - transpose(s.b));
    R = transpose(s.A) + s.B;      CHECK_MATRIX_EQ(R, transpose(s.a) + s.b);
    R = transpose(s.A) - s.B;      CHECK_MATRIX_EQ(R, transpose(s.a) - s.b);
    R = transpose(s.A) * 2.0;      CHECK_MATRIX_EQ(R, transpose(s.a) * 2.0);
    R = transpose(s.A) / 2;        CHECK_MATRIX_EQ(R, transpose(s.a) / 2.0);
    R = 2.0 - transpose(s.A);      CHECK_MATRIX_EQ(R, 2.0 - transpose(s.a));
    R = transpose(s.A) * s.B;      CHECK_MATRIX_EQ(R, transpose(s.a) * s.b);
    R = s.A * transpose(s.B);      CHECK_MATRIX_EQ(R, s.a * transpose(s.b));
}

TEST_CASE("Transpose: of expressions")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = transpose(s.A + s.B);          CHECK_MATRIX_EQ(R, transpose(s.a + s.b));
    R = transpose(s.A * s.B);          CHECK_MATRIX_EQ(R, transpose(s.a * s.b));
    R = transpose(transpose(s.A));     CHECK_MATRIX_EQ(R, s.a);
    R = s.B + s.A + transpose(s.A * s.B + s.C);
    CHECK_MATRIX_EQ(R, s.b + s.a + transpose(s.a * s.b + s.c));

    RowMajorMatrix<double> W(2, 3), V(2, 3), T(3, 2);
    matl_test::fill_sequence(W);
    matl_test::fill_sequence(V, 5.0);
    T = transpose(W + V);              CHECK_MATRIX_EQ(T, transpose(Ref::from(W) + Ref::from(V)));
}

/********************************* Aliasing (target on RHS) ***********************************/

TEST_CASE("Aliasing: target appears in elementwise expressions")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = s.A;  Ref r = s.a;
    R = R + s.B;           r = r + s.b;           CHECK_MATRIX_EQ(R, r);
    R = s.B - R;           r = s.b - r;           CHECK_MATRIX_EQ(R, r);
    R = R * 2.0 + R;       r = r * 2.0 + r;       CHECK_MATRIX_EQ(R, r);
}

TEST_CASE("Aliasing: target appears in a matrix product")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    R = s.C;  Ref r = s.c;
    R = s.A * R;           r = s.a * r;           CHECK_MATRIX_EQ(R, r);
    R = R * s.B;           r = r * s.b;           CHECK_MATRIX_EQ(R, r);
    R = R * R;             r = r * r;             CHECK_MATRIX_EQ(R, r);
    R = s.A + s.B * R - s.B; r = s.a + s.b * r - s.b; CHECK_MATRIX_EQ(R, r);
}

TEST_CASE("Aliasing: target appears transposed")
{
    Square3 s;
    RowMajorMatrix<double> R(3, 3);
    Ref r;
    R = s.A; r = s.a;  R = transpose(R);           r = transpose(r);           CHECK_MATRIX_EQ(R, r);
    R = s.A; r = s.a;  R = transpose(R + R);       r = transpose(r + r);       CHECK_MATRIX_EQ(R, r);
    R = s.A; r = s.a;  R = R + transpose(R);       r = r + transpose(r);       CHECK_MATRIX_EQ(R, r);
    R = s.A; r = s.a;  R = transpose(R) - R;       r = transpose(r) - r;       CHECK_MATRIX_EQ(R, r);
    R = s.A; r = s.a;  R = transpose(R) * 0.5;     r = transpose(r) * 0.5;     CHECK_MATRIX_EQ(R, r);
    R = s.A; r = s.a;  R = transpose(R) / 2;       r = transpose(r) / 2.0;     CHECK_MATRIX_EQ(R, r);
    R = s.A; r = s.a;  R = transpose(R) + 1.0;     r = transpose(r) + 1.0;     CHECK_MATRIX_EQ(R, r);
    R = s.A; r = s.a;  R = transpose(R) - 1.0;     r = transpose(r) - 1.0;     CHECK_MATRIX_EQ(R, r);
    R = s.A; r = s.a;  R = 1.0 - transpose(R);     r = 1.0 - transpose(r);     CHECK_MATRIX_EQ(R, r);
    R = s.A; r = s.a;  R = transpose(R) * R;       r = transpose(r) * r;       CHECK_MATRIX_EQ(R, r);
}

/****************************** Fixed-size and mixed operands *********************************/

TEST_CASE("Expressions: fixed-size operands")
{
    FixedSizeRowMajorMatrix<double, 3, 3> A, B, R;
    matl_test::fill_sequence(A);
    matl_test::fill_sequence(B, -2.0);
    Ref a = Ref::from(A), b = Ref::from(B);
    R = A + B;              CHECK_MATRIX_EQ(R, a + b);
    R = A * B;              CHECK_MATRIX_EQ(R, a * b);
    R = A * B + A;          CHECK_MATRIX_EQ(R, a * b + a);
    R = transpose(A) - B;   CHECK_MATRIX_EQ(R, transpose(a) - b);
    R = A * R;              CHECK_MATRIX_EQ(R, a * (transpose(a) - b));

    FixedSizeRowMajorMatrix<double, 2, 3> W;
    FixedSizeRowMajorMatrix<double, 3, 2> T;
    matl_test::fill_sequence(W);
    T = transpose(W);       CHECK_MATRIX_EQ(T, transpose(Ref::from(W)));
}

TEST_CASE("Expressions: mixed fixed-size and dynamic operands")
{
    FixedSizeRowMajorMatrix<double, 2, 3> F;
    RowMajorMatrix<double> D(3, 2), R(2, 2), S(2, 3);
    matl_test::fill_sequence(F);
    matl_test::fill_sequence(D, 7.0);
    R = F * D;              CHECK_MATRIX_EQ(R, Ref::from(F) * Ref::from(D));
    S = F + transpose(D);   CHECK_MATRIX_EQ(S, Ref::from(F) + transpose(Ref::from(D)));
}
