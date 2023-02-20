//#ifndef _MATL_ROWMAT_ATOMIC_EXPR_PROTOTYPE__HPP
//#define _MATL_ROWMAT_ATOMIC_EXPR_PROTOTYPE__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method declarations and implementations for atomic row matrix 
 * expressions. Atomic expressions are operations such as matrix addition and scalar 
 * multiplication, in which each element of the result depends only on elements that share the 
 * same index. Atomic expressions are garaunteed to support lazy execution  
 **************************************************************************************************/

//#include "_matl_forward_declaraions_.hpp"
//#include "_matl_type_traits_.hpp"

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

        // Single index data getter
        dtype operator() (index_t i) const;
        
        // Double index data getter
        dtype operator() (index_t i, index_t j) const;
        
        // Evaluate expression to target matrix
        template<typename MAT_T> void eval(Matrix<MAT_T>& target) const;

        // Number of rows in Matrix
        index_t rows() const;

        // Number of columns in Matrix
        index_t cols() const;

        // Number of elements in Matrix
        index_t size() const;

        // Check if there is memory conflict with target matrix
        template<typename MAT_T>
        bool memory_conflict(Matrix<MAT_T> const&) const;
        
    };
}


//#endif