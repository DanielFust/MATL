/**************************************************************************************************
 * Example 4: storage order and memory
 *
 * Row-major and column-major matrices, matrices over memory you own (a buffer, or an array
 * filled by Fortran or C code), passing a matrix's elements to other code, and fixed-size
 * matrices.
 *
 * Build and run from the repository root:
 *
 *     make examples && ./build/examples/04_storage_and_memory
 **************************************************************************************************/
#include <iostream>
#include <vector>
#include <matl.hpp>

using namespace matl;

// Stands for a routine of a C or Fortran library that works on a plain array
void scale_array(double* data, std::size_t n, double factor)
{
    for (std::size_t k = 0; k < n; ++k) data[k] *= factor;
}

int main()
{
    /********************************** Storage order *******************************************/

    // Matrix<T> stores its elements row by row (row-major). ColMajorMatrix<T> stores them column
    // by column, as Fortran, LAPACK and MATLAB do. Both are used the same way: A(i,j) is row i,
    // column j, and initializer lists give the rows.
    Matrix<double>         R = {{1, 2, 3}, {4, 5, 6}};
    ColMajorMatrix<double> C = {{1, 2, 3}, {4, 5, 6}};

    // Only the order in memory differs, which the single index A(k) and data() follow:
    std::cout << "row-major elements in memory:    ";
    for (index_t k = 1; k <= R.size(); ++k) std::cout << R(k) << ' ';
    std::cout << "\ncolumn-major elements in memory: ";
    for (index_t k = 1; k <= C.size(); ++k) std::cout << C(k) << ' ';
    std::cout << "\n\n";

    // Matrices of both orders mix freely in expressions and assignments
    Matrix<double> D = R + C;
    C = 2.0 * R;
    std::cout << "R + C =\n" << D << "\n\n";

    /*************************** Matrices over memory you own ***********************************/

    // A matrix can use memory it does not own, for example a buffer received through MPI or one
    // block holding many matrices. The matrix reads and writes the buffer in place and never
    // frees it, so the buffer must outlive the matrix.
    std::vector<double> buffer(2 * 4, 0.0);
    Matrix<double> P(2, 2, buffer.data());            // uses buffer[0..3]
    Matrix<double> Q(2, 2, buffer.data() + 4);        // uses buffer[4..7]
    P = eye<double>(2);
    Q = P * 3.0 + 1.0;                                 // the result is written into the buffer
    std::cout << "buffer: ";
    for (double value : buffer) std::cout << value << ' ';
    std::cout << "\n\n";

    // An array filled by Fortran (or LAPACK) is column by column; wrap it in a ColMajorMatrix to
    // use it without copying
    double fortran_array[6] = {1, 4,                   // column 1
                               2, 5,                   // column 2
                               3, 6};                  // column 3
    ColMajorMatrix<double> F(2, 3, fortran_array);
    std::cout << "the Fortran array as a matrix:\n" << F << "\n\n";

    // data() points to a matrix's elements in storage order, to pass them to such code
    scale_array(F.data(), F.size(), 10.0);
    std::cout << "after scale_array: F(2,3) = " << F(2, 3) << "\n\n";

    /******************************** Fixed-size matrices ***************************************/

    // Matrix<T, M, N> holds its M*N elements inside the object (on the stack), which suits many
    // small matrices: 3x3 rotations, 4x4 transformations, small blocks. Mismatched dimensions are
    // compile errors, e.g. adding a Matrix<double, 2, 3> to a Matrix<double, 3, 2>.
    Matrix<double, 3, 3> rotation = {{0, -1, 0},
                                     {1,  0, 0},
                                     {0,  0, 1}};       // 90 degrees about z
    Matrix<double, 3, 1> point = {{1}, {0}, {0}};
    Matrix<double, 3, 1> rotated = rotation * point;
    std::cout << "rotated point:\n" << rotated << "\n\n";

    // Converting between fixed-size and dynamic matrices copies the elements
    Matrix<double> dynamic = rotation;
    std::cout << "the rotation as a dynamic matrix is " << dynamic.rows() << " x " << dynamic.cols() << "\n";

    return 0;
}
