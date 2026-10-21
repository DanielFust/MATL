#ifndef _MATL_INTERFACES_DECLARATIONS__HPP
#define _MATL_INTERFACES_DECLARATIONS__HPP

namespace matl 
{
    // Matrix Algebra Object
    template<typename MAO> class _MATRIX_ALGEBRA_OBJECT_;

    // Matrix Base Class
    template<typename MAT_T> class _MATRIX_;

    // Base class for matrix expressions
    template<typename EXPR_T> class _EXPRESSION_;

    // Interface for elementwise matrix data access
    template<typename EW_T> class _ELEMENTWISE_;

    // Interface for fixed-size matrices and expressions
    template<typename FS> class _FIXED_SIZE_;

    // Interface for fixed-size matrices
    template<typename FSM> class _FIXED_SIZE_MATRIX_;
}
#endif