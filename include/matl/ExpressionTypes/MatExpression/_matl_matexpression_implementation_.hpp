#ifndef _MATL_MATEXPRESSION_IMPLEMENTATION__HPP
#define _MATL_MATEXPRESSION_IMPLEMENTATION__HPP

#include <iostream>
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"

namespace matl
{
  
   /***********************************************************************************************
    * void eval(Matrix& target) const
    * 
    * Evaluate expression to target matrix. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Atomic expression type
    *      MAT_T            Matrix type
    * 
    * Parameters:
    *     Matrix target     Target matrix
    * Returns:  
    ***********************************************************************************************/
    template<typename EXPR_T>
    template<typename MAT_T> 
    void _MATEXPRESSION_<EXPR_T>::eval(Matrix<MAT_T>& target) const
    {static_cast<EXPR_T const&>(*this).eval(target);}




   /***********************************************************************************************
    * Matrix eval() const
    * 
    * Evaluate expression to new matrix. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Atomic expression type
    * 
    * Parameters:
    * Returns:
    *     Matrix of expression's preferred type  
    ***********************************************************************************************/
    template<typename EXPR_T>
    Matrix<typename _MATEXPRESSION_<EXPR_T>::matrix_type> _MATEXPRESSION_<EXPR_T>::eval() const
    {return static_cast<EXPR_T const&>(*this).eval();}




   /***********************************************************************************************
    * index_t rows() const
    * 
    * Number of rows in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Expression type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename EXPR_T>
    index_t _MATEXPRESSION_<EXPR_T>::rows() const
    {return static_cast<EXPR_T const&>(*this).rows();}





   /***********************************************************************************************
    * index_t cols() const
    * 
    * Number of columns in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Expression type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename EXPR_T>
    index_t _MATEXPRESSION_<EXPR_T>::cols() const
    {return static_cast<EXPR_T const&>(*this).cols();}




   /***********************************************************************************************
    * index_t size() const
    * 
    * Number of elements in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Expression type
    * 
    * Parameters:
    ***********************************************************************************************/  
    template<typename EXPR_T>
    index_t _MATEXPRESSION_<EXPR_T>::size() const
    {return static_cast<EXPR_T const&>(*this).size();}




   /***********************************************************************************************
    * bool memory_conflict(Matrix const&) const
    * 
    * Check if there is memory conflict between expression and target matrix
    * 
    * Template Parameters:
    *     EXPR_T           Expression type
    * 
    * Parameters:
    *     target           Target matrix
    * 
    * Returns:
    *     'true' if direct evaluation corrupts memory 'false' otherwise
    ***********************************************************************************************/  
    template<typename EXPR_T>
    template<typename MAT_T>
    bool _MATEXPRESSION_<EXPR_T>::memory_conflict(Matrix<MAT_T> const& target) const
    {return static_cast<EXPR_T const&>(*this).memory_conflict(target);}
}

#endif