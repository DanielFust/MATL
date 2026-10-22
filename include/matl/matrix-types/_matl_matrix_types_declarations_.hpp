#ifndef MATL_MATRIX_TYPES_DECLARATIONS_HPP
#define MATL_MATRIX_TYPES_DECLARATIONS_HPP
namespace matl
{
    // Row Major Matrix
    template<class T> class RowMajorMatrix;

    // Fixed-size Row Major Matrix
    template<class T, index_t M, index_t N> class FixedSizeRowMajorMatrix;
}
#endif
