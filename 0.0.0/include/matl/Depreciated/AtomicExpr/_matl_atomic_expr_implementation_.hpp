//#ifndef _MATL_ATOMIC_EXPR__HPP
//#define _MATL_ATOMIC_EXPR__HPP

#include <iostream>
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"

namespace matl
{

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
    typename _ATOMIC_EXPR_<EXPR_T>::dtype _ATOMIC_EXPR_<EXPR_T>::operator() (index_t i) const
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
    typename _ATOMIC_EXPR_<EXPR_T>::dtype _ATOMIC_EXPR_<EXPR_T>::operator() 
        (index_t i, index_t j) const
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
    template<typename EXPR_T>
    template<typename MAT_T>
    void _ATOMIC_EXPR_<EXPR_T>::eval(Matrix<MAT_T>& target) const
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



   /***********************************************************************************************
    * index_t cols() const
    * 
    * Number of columns in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Atomic expression type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename EXPR_T>
    index_t _ATOMIC_EXPR_<EXPR_T>::cols() const
    {return static_cast<EXPR_T const&>(*this).cols();}




   /***********************************************************************************************
    * index_t size() const
    * 
    * Number of elements in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Atomic expression type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename EXPR_T>
    index_t _ATOMIC_EXPR_<EXPR_T>::size() const
    {return static_cast<EXPR_T const&>(*this).size();}



   /***********************************************************************************************
    * bool memory_conflict(Matrix const&) const
    * 
    * Check if there is memory conflict between expression and target matrix
    * 
    * Template Parameters:
    *     EXPR_T           Expression type
    *     MAT_T            Target matrix type
    * 
    * Parameters:
    *     target           Target matrix
    * 
    * Returns:
    *     'true' if direct evaluation corrupts memory 'false' otherwise
    ***********************************************************************************************/  
    template<typename EXPR_T>
    template<typename MAT_T>
    bool _ATOMIC_EXPR_<EXPR_T>::memory_conflict(Matrix<MAT_T> const& target) const {return false;} 

}

//#endif