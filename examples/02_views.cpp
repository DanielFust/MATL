/**************************************************************************************************
 * Example 2: views
 *
 * Selecting rows, columns and blocks of a matrix with range, all and last (like MATLAB's
 * A(1:3, :) or A(end, 2:end)), reading and writing through them, and row operations.
 *
 * A view does not copy anything: it refers to the elements of the matrix, so writing to a view
 * writes the matrix.
 *
 * Build and run from the repository root:
 *
 *     make examples && ./build/examples/02_views
 **************************************************************************************************/
#include <iostream>
#include <matl.hpp>

using namespace matl;

int main()
{
    Matrix<double> A = {{ 1,  2,  3,  4},
                        { 5,  6,  7,  8},
                        { 9, 10, 11, 12},
                        {13, 14, 15, 16}};

    /********************************* Selecting elements ***************************************/

    //   MATL                          MATLAB
    //   A(range(1, 2), all)           A(1:2, :)        rows 1 to 2
    //   A(all, 3)                     A(:, 3)          column 3
    //   A(last, range(2, last))       A(end, 2:end)    the last row, from column 2
    //   A(range(1, 2, last), all)     A(1:2:end, :)    every other row
    //   A(range(last, -1, 1), all)    A(end:-1:1, :)   the rows in reverse order
    std::cout << "rows 1-2:\n"          << A(range(1, 2), all)        << "\n\n";
    std::cout << "column 3:\n"          << A(all, 3)                  << "\n\n";
    std::cout << "last row, 2:end:\n"   << A(last, range(2, last))    << "\n\n";
    std::cout << "every other row:\n"   << A(range(1, 2, last), all)  << "\n\n";
    std::cout << "rows reversed:\n"     << A(range(last, -1, 1), all) << "\n\n";

    // A view can be used anywhere a matrix can: in expressions, products and reductions
    std::cout << "sum of column 2 = " << sum(A(all, 2)) << "\n";
    std::cout << "block * block =\n"
              << A(range(1, 2), range(1, 2)) * A(range(3, 4), range(3, 4)) << "\n\n";

    /********************************* Writing through views ************************************/

    Matrix<double> B = A;

    B(all, 1) = 0.0;                                        // set a column
    B(1, all) = transpose(A(all, 4));                       // a row from a column of A
    B(range(3, 4), range(3, 4)) = eye<double>(2) * 100.0;   // a block from an expression
    B(2, all) += 2.0 * B(3, all);                           // a row operation
    swap(B(3, all), B(4, all));                             // exchange two rows
    std::cout << "B after writing through views:\n" << B << "\n\n";

    // An assignment may read the elements it overwrites; MATL notices and gives the result as if
    // the right-hand side were computed first. Here every row moves down by one:
    B(range(2, last), all) = B(range(1, last - 1), all);
    std::cout << "B with its rows shifted down:\n" << B << "\n\n";

    /********************************* Views and copies *****************************************/

    // 'auto' keeps a view: v refers to column 2 of A, and writing to v writes A
    auto v = A(all, 2);
    v = -1.0;
    std::cout << "A after v = -1:\n" << A << "\n\n";

    // A Matrix made from a view is a copy: changing it does not change A
    Matrix<double> column = A(all, 3);
    column = 0.0;
    std::cout << "A(1,3) is still " << A(1, 3) << "\n";

    // A view is valid only while its matrix exists and keeps its memory (it must not be
    // destroyed, moved from or resized). A view of a temporary keeps the temporary alive:
    auto diagonal_block = eye<double>(4)(range(2, 3), range(2, 3));
    std::cout << "a view of a temporary identity:\n" << diagonal_block << "\n\n";

    /****************************** Example: Gaussian elimination ********************************/

    // Views make row operations read like the algorithm. Eliminate below the diagonal of M
    // (without pivoting; see 03_linear_systems.cpp for the library's LU factorization):
    Matrix<double> M = {{ 2,  1, 1},
                        { 4,  3, 3},
                        { 8,  7, 9}};
    for (index_t k = 1; k < M.rows(); ++k)
        for (index_t i = k + 1; i <= M.rows(); ++i)
            M(i, range(k, last)) -= (M(i, k) / M(k, k)) * M(k, range(k, last));
    std::cout << "M reduced to upper triangular form:\n" << M << "\n";

    return 0;
}
