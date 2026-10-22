#ifndef MATL_MATRIX_TYPES_DECLARATIONS_HPP
#define MATL_MATRIX_TYPES_DECLARATIONS_HPP
namespace matl
{
    // Storage order of a matrix: RowMajor (the rows are contiguous in memory; the default) or
    // ColMajor (the columns are contiguous, as in Fortran, LAPACK and MATLAB). For expressions,
    // Mixed means operands with different storage orders, and Any is a scalar (no order of its own).
    enum class StorageOrder
    {
        RowMajor,
        ColMajor,
        Mixed,
        Any
    };

    // Dynamic (heap-allocated) matrix with the given storage order
    template<class T, StorageOrder ORDER> class DynamicMatrix;

    // Fixed-size M x N matrix with the given storage order
    template<class T, index_t M, index_t N, StorageOrder ORDER> class FixedSizeMatrix;

    // The matrix types by storage order
    template<class T> using RowMajorMatrix = DynamicMatrix<T, StorageOrder::RowMajor>;
    template<class T> using ColMajorMatrix = DynamicMatrix<T, StorageOrder::ColMajor>;
    template<class T, index_t M, index_t N> using FixedSizeRowMajorMatrix = FixedSizeMatrix<T, M, N, StorageOrder::RowMajor>;
    template<class T, index_t M, index_t N> using FixedSizeColMajorMatrix = FixedSizeMatrix<T, M, N, StorageOrder::ColMajor>;

    // View of rows and columns of a matrix (T is const for views of const matrices)
    template<class T, StorageOrder ORDER> class MatrixView;

    // Rows of a matrix, as a nested initializer list (e.g. {{1, 2}, {3, 4}})
    template<class T>
    using matrix_rows = std::initializer_list< std::initializer_list<T> >;
}
#endif
