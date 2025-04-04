#ifndef _MATL_MATRIX_TYPES_DECLARATIONS__HPP
#define _MATL_MATRIX_TYPES_DECLARATIONS__HPP
namespace matl
{
    // Row Major Matrix
    template<class T> class RowMajorMatrix;

    // Fixed-size Row Major Matrix
    template<class T, index_t M, index_t N> class FixedSizeRowMajorMatrix;
}
#endif
