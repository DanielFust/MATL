#ifndef _MATL_ROWMAT_ATOMIC_EXPR__HPP
#define _MATL_ROWMAT_ATOMIC_EXPR__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method declarations and implementations for atomic row matrix 
 * expressions. Atomic expressions are operations such as matrix addition and scalar 
 * multiplication, in which each element of the result depends only on elements that share the 
 * same index. Atomic expressions are garaunteed to support lazy execution  
 **************************************************************************************************/

#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"

namespace matl
{
   /**********************************************************************************************
    * class _ROWMAT_ATOMIC_EXPR_ 
    * 
    * Atomic (element-wise) expressions for row major matrices. Atomic expressions support 
    * lazy execution without requiring memory allocations. Each element of the resultant matrix
    * can be calculated via indexing the _ROWMAT_ATOMIC_EXPR_.
    *    
    * Template Arguments:
    *     EXPR_T       (_ROWMAT_ATOMIC_EXPR_)
    **********************************************************************************************/
    template<typename EXPR_T>
    class _ROWMAT_ATOMIC_EXPR_ : public _ATOMIC_EXPR_< _ROWMAT_ATOMIC_EXPR_<EXPR_T> >
    {
      public:
        // Expression data type
        typedef typename TypeTraits<EXPR_T>::dtype    dtype;

        // Matrix return type
        //typedef typename TypeTraits<EXPR_T>::mat_type mat_type; 

        // Single index data getter
        dtype operator() (index_t i) const
        {return static_cast<EXPR_T const&>(*this)(i);}
        
        
        // Double index data getter
        dtype operator() (index_t i, index_t j) const
        {return static_cast<EXPR_T const&>(*this)(i,j);}
        
        
        // Evaluate expression to target matrix
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const
        {
            static_cast<EXPR_T const&>(*this).eval(static_cast<MAT_T&>(target));
            return;
        }


        // Number of rows in Matrix
        index_t rows() const
        {return static_cast<EXPR_T const&>(*this).rows();}
        

        // Number of columns in Matrix
        index_t cols() const
        {return static_cast<EXPR_T const&>(*this).cols();}
        

        // Number of elements in Matrix
        index_t size() const
        {return static_cast<EXPR_T const&>(*this).size();}
        
    };
}


#endif