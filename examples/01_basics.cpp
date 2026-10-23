/**************************************************************************************************
 * Example 1: the basics
 *
 * Creating matrices, reading and writing elements, arithmetic, matrix products, elementwise
 * functions and reductions.
 *
 * Build and run from the repository root:
 *
 *     make examples && ./build/examples/01_basics
 *
 * or compile it directly:  g++ -std=c++20 -I include examples/01_basics.cpp -o basics
 **************************************************************************************************/
#include <iostream>
#include <matl.hpp>

using namespace matl;

int main()
{
    /******************************** Creating matrices ****************************************/

    // Matrix<T> has dimensions chosen at run time. An initializer list gives the rows.
    Matrix<double> A = {{1, 2, 3},
                        {4, 5, 6}};                 // 2x3

    // A sized matrix. Its elements are not initialized, so give them values before reading.
    Matrix<double> B(2, 3);
    B = 0.5;                                        // assigning a scalar sets every element

    // zeros, ones and eye (the identity) make common matrices
    Matrix<double> Z = zeros<double>(2, 3);
    Matrix<double> I = eye<double>(3);

    // Matrix<T, M, N> has dimensions fixed at compile time: it needs no heap memory, and
    // combining fixed-size matrices of the wrong dimensions does not compile.
    Matrix<double, 2, 2> F = {{2, 1},
                              {1, 3}};

    // operator<< writes a matrix with its columns aligned; print() writes it to std::cout
    std::cout << "A =\n" << A << "\n\n";
    std::cout << "Z =\n" << Z << "\n\nI =\n" << I << "\n\n";

    /*********************************** Elements ***********************************************/

    // Indices start at 1, as in MATLAB and Fortran. 'last' is the last row or column.
    A(2, 3) = 10.0;
    std::cout << "A(1,2) = " << A(1, 2) << ",  A(last, last) = " << A(last, last) << "\n";
    std::cout << "A is " << A.rows() << " x " << A.cols() << " (" << A.size() << " elements)\n\n";

    /*********************************** Arithmetic *********************************************/

    // An arithmetic expression is computed when it is assigned to a matrix, in a single pass
    // over the elements (2*A + B - 1 creates no intermediate matrices).
    Matrix<double> C = 2.0 * A + B - 1.0;
    std::cout << "2A + B - 1 =\n" << C << "\n\n";

    // '*' between two matrices is the matrix product; transpose(A) is the transpose.
    Matrix<double> P = A * transpose(A);           // (2x3)(3x2) = 2x2
    std::cout << "A * transpose(A) =\n" << P << "\n\n";

    // Fixed-size and dynamic matrices mix freely
    Matrix<double> Q = F * P + F;
    std::cout << "F * P + F =\n" << Q << "\n\n";

    // Compound assignment updates a matrix in place
    C += B;                                         // C = C + B
    C *= 2.0;                                       // every element times 2
    P *= F;                                         // P = P * F (matrix product)

    /***************************** Elementwise operations ***************************************/

    // C++ has no .* or ./ operators, so MATL uses MATLAB's function names for them
    Matrix<double> W = {{1, -1,  2},
                        {0,  3, -2}};
    Matrix<double> E = times(A, W);                // elementwise product,  MATLAB A .* W
    Matrix<double> R = rdivide(W, A);              // elementwise quotient, MATLAB W ./ A
    Matrix<double> S = pow(A, 2);                  // elementwise power,    MATLAB A .^ 2

    // Math functions apply to every element (abs, sqrt, exp, log, sin, cos, ...)
    Matrix<double> G = sqrt(abs(W)) + exp(-A);
    std::cout << "A .* W =\n" << E << "\n\nW ./ A =\n" << R << "\n\nA .^ 2 =\n" << S
              << "\n\nsqrt(|W|) + exp(-A) =\n" << G << "\n\n";

    /*********************************** Reductions *********************************************/

    // Each reduction gives one number for the whole matrix (unlike MATLAB's sum(A), which gives
    // the column sums)
    std::cout << "sum(A) = " << sum(A) << ",  mean(A) = " << mean(A)
              << ",  min(A) = " << min(A) << ",  max(A) = " << max(A) << "\n";
    std::cout << "norm_fro(A) = " << norm_fro(A) << ",  norm_1(A) = " << norm_1(A)
              << ",  norm_inf(A) = " << norm_inf(A) << ",  trace(F) = " << trace(F) << "\n";

    // Expressions can be reduced directly, without storing them in a matrix first
    std::cout << "max(|A - C|) = " << max(abs(A - C)) << ",  dot(A, B) = " << dot(A, B) << "\n\n";

    /**************************** A note on 'auto' ************************************************/

    // An expression stored with 'auto' is a formula, not a matrix: it reads A and B whenever it
    // is evaluated. Name the type to get a matrix holding the values.
    auto formula = A + B;                           // a formula over A and B
    Matrix<double> value = A + B;                   // a matrix
    A = 0.0;
    Matrix<double> later = formula;                 // evaluated now: uses the new A
    std::cout << "value(1,1) = " << value(1, 1) << ",  formula evaluated after A = 0: "
              << later(1, 1) << "\n";

    return 0;
}
