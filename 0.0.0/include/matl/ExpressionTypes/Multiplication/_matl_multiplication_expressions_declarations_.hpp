//#ifndef _MATL_MULTIPLICATION_EXPRESSIONS_DECLARATIONS__HPP
//#define _MATL_MULTIPLICATION_EXPRESSIONS_DECLARATIONS__HPP

namespace matl
{
    // Generic Matrix - Scalar Multiplication
    template<typename AE, class T>           class _GENERIC_MAT_SCALAR_PROD_;

    // RowMatrix - RowMatrix Multiplication
    template<typename RMAE1, typename RMAE2> class _ROWMAT_PROD_;

    // RowMatrix - Scalar Multiplication
    template<typename RMAE, class T>         class _ROWMAT_SCALAR_PROD_; 

    // Generic Matrix - Matrix Multiplication
    template<typename MAT1, typename MAT2>   class _GENERIC_MAT_PROD_;
}

//#endif