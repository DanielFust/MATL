/**************************************************************************************************
 * Example 7: sparse matrices
 *
 * Assembling sparse matrices element by element (as a finite-element code does), products with
 * vectors, changing values in a fixed sparsity pattern, conversions, and the CSR arrays for other
 * libraries.
 *
 * The problem: the 1D Poisson equation -u'' = f on [0, 1] with linear finite elements. Each element
 * [x_e, x_e+1] contributes a 2x2 block to the stiffness matrix K and the mass matrix M; blocks of
 * neighboring elements overlap at their shared node, and SparseBuilder sums the overlapping
 * entries.
 *
 * Build and run from the repository root:
 *
 *     make examples && ./build/examples/07_sparse_matrices
 **************************************************************************************************/
#include <cstdint>
#include <iostream>
#include <matl.hpp>

using namespace matl;

int main()
{
    const index_t elements = 5;
    const index_t nodes = elements + 1;
    const double h = 1.0 / double(elements);               // element length

    /************************************* Assembly *********************************************/

    // Element matrices of a linear element of length h
    const Matrix<double, 2, 2> k_e = Matrix<double, 2, 2>{{ 1, -1},
                                                          {-1,  1}} * (1.0 / h);
    const Matrix<double, 2, 2> m_e = Matrix<double, 2, 2>{{2, 1},
                                                          {1, 2}} * (h / 6.0);

    // Add each element's block at its two nodes (e, e+1); overlapping entries are summed
    SparseBuilder<double> K_builder(nodes, nodes), M_builder(nodes, nodes);
    for (index_t e = 1; e <= elements; ++e)
        for (index_t a = 1; a <= 2; ++a)
            for (index_t b = 1; b <= 2; ++b)
            {
                K_builder.add(e + a - 1, e + b - 1, k_e(a, b));
                M_builder.add(e + a - 1, e + b - 1, m_e(a, b));
            }
    CSRMatrix<double> K = K_builder.build();
    CSRMatrix<double> M = M_builder.build();

    // A sparse matrix prints its stored entries; dense(K) shows the whole (small) matrix
    std::cout << "K has " << K.nnz() << " stored entries of " << K.rows() * K.cols() << ":\n" << dense(K) << "\n\n";

    /*********************************** Products ***********************************************/

    // Sparse times dense is an expression with a dense result, like a dense product
    Matrix<double> u(nodes, 1);
    for (index_t i = 1; i <= nodes; ++i)
    {
        const double x = double(i - 1) * h;
        u(i) = x * (1.0 - x);                              // u = x(1 - x), so -u'' = 2
    }
    Matrix<double> Ku = K * u;
    std::cout << "K u =\n" << Ku << "\n\n";

    // Interior rows of K u equal the load M f for f = 2 (the boundary rows carry the boundary terms)
    Matrix<double> f = 2.0 * ones<double>(nodes, 1);
    Matrix<double> r = K * u - M * f;
    std::cout << "residual K u - M f at the interior nodes: " << max(abs(r(range(2, last - 1), all))) << "\n";

    // transpose(K) * x multiplies by the transpose without forming it (K is symmetric here)
    std::cout << "max |K^T u - K u| = " << max(abs(transpose(K) * u - Ku)) << "\n\n";

    /******************************** Changing values *******************************************/

    // The pattern is fixed once built, but the values can change, e.g. to form the time-stepping
    // matrix M + dt K, or to re-assemble with new coefficients
    const double dt = 0.01;
    CSRMatrix<double> A = M * 1.0;                         // a copy with M's pattern (the same as K's)
    for (index_t i = 1; i <= nodes; ++i)
        for (index_t j = 1; j <= nodes; ++j)
            if (A.is_stored(i, j)) A.entry(i, j) += dt * K(i, j);
    std::cout << "M + dt K (row 2): " << A(2, 1) << "  " << A(2, 2) << "  " << A(2, 3) << "\n\n";

    /********************************** The CSR arrays ******************************************/

    // The arrays are 0-based, as other libraries (SciPy, MKL, PETSc, SuiteSparse) expect; element
    // access is 1-based, like the rest of MATL
    std::cout << "row_offsets:    ";
    for (auto p : K.row_offsets()) std::cout << p << ' ';
    std::cout << "\ncolumn_indices: ";
    for (auto c : K.column_indices()) std::cout << c << ' ';
    std::cout << "\n\n";

    // 32-bit indices use less memory and make products faster, for up to 2^31 stored entries
    CSRMatrix<double, std::int32_t> K32 = sparse<std::int32_t>(dense(K));
    std::cout << "K with 32-bit indices: " << K32.nnz() << " entries, K32(2, 2) = " << K32(2, 2) << "\n";

    // Dense algorithms do not accept sparse matrices; lu(K) would not compile, and say why.
    // For a small matrix, convert explicitly:
    auto F = lu(dense(K) + eye<double>(nodes));
    std::cout << "rcond(K + I) = " << F.rcond() << "\n";

    return 0;
}
