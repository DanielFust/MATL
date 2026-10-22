#ifndef MATL_MATRIX_TYPES_DECLARATIONS_HPP
#define MATL_MATRIX_TYPES_DECLARATIONS_HPP
namespace matl
{
    // Row Major Matrix
    template<class T> class RowMajorMatrix;

    // Fixed-size Row Major Matrix
    template<class T, index_t M, index_t N> class FixedSizeRowMajorMatrix;

    // Rows of a matrix, as a nested initializer list (e.g. {{1, 2}, {3, 4}})
    template<class T>
    using matrix_rows = std::initializer_list< std::initializer_list<T> >;
}
#endif
