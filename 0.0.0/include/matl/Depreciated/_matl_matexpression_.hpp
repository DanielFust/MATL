#ifndef _MATL_MATEXPRESSION__HPP
#define _MATL_MATEXPRESSION__HPP

#include <iostream>
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"

namespace matl
{
   /**********************************************************************************************
    * class class _MATEXPRESSION_
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

        // Resultant matrix type (e.g. row matrix, compressed sparse column, etc.)
        //typedef typename TypeTraits<EXPR_T>::mat_type mat_type;

        // Evaluate expression
        ///mat_type eval() const
        //{return static_cast<EXPR_T const&>(*this).eval();}

        // Evaluate Expression to target matrix
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const
        {static_cast<EXPR_T const&>(*this).eval(static_cast<MAT_T&>(target));}

        // Get number of rows of resultant matrix
        index_t rows() const
        {return static_cast<EXPR_T const&>(*this).rows();}

        // Get number of columns of resultant matrix
        index_t cols() const
        {return static_cast<EXPR_T const&>(*this).cols();}

        // Get number of elements in resultant matrix
        index_t size() const
        {return static_cast<EXPR_T const&>(*this).size();}
              
    };
}

#endif