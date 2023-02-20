#ifndef _MATL_ATOMIC_EXPR__HPP
#define _MATL_ATOMIC_EXPR__HPP

#include <iostream>
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"

namespace matl
{
   /**********************************************************************************************
    * class _ATOMIC_EXPR_
    * 
    * Base class for atomic matrix expression templates. Atomic expressions are operations that
    * can be performed independently on all elements of a matrix concurrently.
    * 
    * Template Parameters:
    *     EXPR_T         (_ATOMIC_EXPR_) Expression type
    **********************************************************************************************/
    template<typename EXPR_T>
    class _ATOMIC_EXPR_ : public _MATEXPRESSION_< _ATOMIC_EXPR_<EXPR_T> >
    {
      public:

        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename TypeTraits<EXPR_T>::dtype dtype;

        // Resultant matrix type (e.g. row matrix, compressed sparse column, etc.)
        //typedef typename TypeTraits<EXPR_T>::mat_type mat_type;

        // Single index data getter
        dtype operator() (index_t i) const
        {return static_cast<EXPR_T const&>(*this)(i);}

        // Double index data getter
        dtype operator() (index_t i, index_t j) const
        {return static_cast<EXPR_T const&>(*this)(i,j);}

        // Evaluate expression
        //mat_type eval() const
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



   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Single index getter. Downcasts to derived type.
    * 
    *  Template Parameters
    *      EXPR_T           Atomic expression type
    * 
    * Parameters
    *     index_t i         index  
    ***********************************************************************************************/
    template<typename EXPR_T>
    _ATOMIC_EXPR_<EXPR_T>::dtype _ATOMIC_EXPR_::operator() (index_t i) const
    {return static_cast<EXPR_T const&>(*this)(i);}



   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Double index getter. Downcasts to derived type.
    * 
    *  Template Parameters
    *      EXPR_T           Atomic expression type
    * 
    * Parameters
    *     index_t i         Row index  
    *     index_t j         Column index  
    ***********************************************************************************************/
    template<typename EXPR_T>
    _ATOMIC_EXPR_<EXPR_T>::dtype _ATOMIC_EXPR_::(index_t i, index_t j) const
    {return static_cast<EXPR_T const&>(*this)(i,j);}




   /***********************************************************************************************
    * void eval(Matrix& target) const
    * 
    * Evaluate expression to target matrix. Downcasts to derived type.
    * 
    *  Template Parameters
    *      EXPR_T           Atomic expression type
    *      MAT_T            Matrix type
    * 
    * Parameters
    *     Matrix target     Target matrix  
    ***********************************************************************************************/
    template<typename EXPR_T, typename MAT_T>
    void _ATOMIC_EXPR_TYPE_<EXPR_T>::eval(Matrix<MAT_T>& target) const
    {static_cast<EXPR_T const&>(*this).eval(static_cast<MAT_T&>(target));}



   /***********************************************************************************************
    * index_t rows() const
    * 
    * Number of rows in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Atomic expression type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename EXPR_T>
    index_t _ATOMIC_EXPR_<EXPR_T>::rows() const
    {return static_cast<EXPR_T const&>(*this).rows();}

        // Get number of columns of resultant matrix
        index_t cols() const
        {return static_cast<EXPR_T const&>(*this).cols();}

        // Get number of elements in resultant matrix
        index_t size() const
        {return static_cast<EXPR_T const&>(*this).size();}

}

#endif