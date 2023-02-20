#ifndef _MATL_EXPRESSION_TYPES_DECLARATIONS__HPP
#define _MATL_EXPRESSION_TYPES_DECLARATIONS__HPP

namespace matl
{
    // Matrix Expression Base class
    template<typename EXPR_T> class _MATEXPRESSION_;

    // Atomic Expressions
    template<typename EXPR_T> class _ATOMIC_EXPR_;

    // Row Matrix/Major Atomic Expressions
    template<typename RMAE_T> class _ROWMAT_ATOMIC_EXPR_;
}

#include "./Addition/_matl_addition_expressions_declarations_.hpp"
#include "./Subtraction/_matl_subtraction_expressions_declarations_.hpp"
#include "./Multiplication/_matl_multiplication_expressions_declarations_.hpp"


#endif