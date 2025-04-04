//#ifndef _MATL_SUBTRACTION_EXPRESSIONS_DECLARATIONS__HPP
//#define _MATL_SUBTRACTION_EXPRESSIONS_DECLARATIONS__HPP

namespace matl
{
    // Generic Matrix Subtraction
    template<typename AE_T1, typename AE_T2> class _GENERIC_MAT_SUB_;

    // RowMatrix - RowMatrix Subtraction
    template<typename RMAE1, typename RMAE2> class _ROWMAT_SUB_;

    // Generic Matrix - Scalar Subtraction
    template<typename EXPR,   class T>        class _GENERIC_MAT_SCALAR_SUB_; // Matrix - Scalar
    template<class T,        typename EXPR>  class _GENERIC_SCALAR_MAT_SUB_; // Scalar - Matrix
}

//#endif