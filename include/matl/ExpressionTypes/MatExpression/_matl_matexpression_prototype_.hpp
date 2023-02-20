//#ifndef _MATL_MATEXPRESSION_PROTOTYPE__HPP
//#define _MATL_MATEXPRESSION_PROTOTYPE__HPP

//#include <iostream>
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"

namespace matl
{
   /**********************************************************************************************
    * class _MATEXPRESSION_ (Prototype)
    * 
    * Base class for expression templates based around the Matrix static interface in the 'matl'
    * namespace. This enables syntactically concise expressions that support nuanced behavior
    * (such as lazy execution, memory managment, cache blocking, etc.) specific to the expression
    * 
    * Template Parameters:
    *     EXPR_T         (_MATEXPRESSION_) Expression type
    **********************************************************************************************/
    template<typename EXPR_T>
    class _MATEXPRESSION_
    {
      public:

        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename TypeTraits<EXPR_T>::dtype dtype;

        // Preferred matrix return type
        typedef typename TypeTraits<EXPR_T>::matrix_type matrix_type;

        // Evaluate Expression to target matrix
        template<typename MAT_T> void eval(Matrix<MAT_T>& target) const;

        // Evaluate Expression to new matrix
        Matrix<matrix_type> eval() const;

        // Get number of rows of resultant matrix
        index_t rows() const;

        // Get number of columns of resultant matrix
        index_t cols() const;

        // Get number of elements in resultant matrix
        index_t size() const;

        // Check if there is memory conflict with target matrix
        template<typename MAT_T>
        bool memory_conflict(Matrix<MAT_T> const&) const;
              
    };
}

//#endif