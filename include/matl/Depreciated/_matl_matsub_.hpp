#ifndef _MATL_MATSUB__HPP
#define _MATL_MATSUB__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains template specializations of the MATSUB() function used to build matrix 
 * subtraction expression templates.
 **************************************************************************************************/
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"
#include "_matl_generic_mat_sub_.hpp"
#include "_matl_rowmat_sub_.hpp"

namespace matl
{

   /**********************************************************************************************
    * _ROWMAT_SUB_ MATSUB(RowMatrix const&, RowMatrix const&)
    * 
    * MATSUB template specialization for RowMatrix inputs
    * 
    * Template Parameters:
    *     T1 (implicit)      First RowMatrix data type
    *     T2 (implicit)      Second RowMatrix data type
    * 
    * Parameters:
    *     mat1               First RowMatrix
    *     mat2               Second RowMatrix
    * 
    * Returns:
    *     _ROWMAT_SUB_       Row matrix subtraction expression template
    **********************************************************************************************/ 
    template<class T1, class T2>
    _ROWMAT_SUB_<RowMatrix<T1>, RowMatrix<T2> > 
        MATSUB(RowMatrix<T1> const& mat1, RowMatrix<T2> const& mat2)
    {return _ROWMAT_SUB_<RowMatrix<T1>, RowMatrix<T2> > (mat1, mat2);}





   /**********************************************************************************************
    * _ROWMAT_SUB_ MATSUB(_ROWMAT_ATOMIC_EXPR_ const&, RowMatrix const&)
    * 
    * MATSUB template specialization for RowMatrix inputs
    * 
    * Template Parameters:
    *     RMAE (implicit)     Row matrix atomic expression type data type
    *     T    (implicit)     RowMatrix data type
    * 
    * Parameters:
    *     expr               Row matrix atomic expression
    *     mat                RowMatrix
    * 
    * Returns:
    *     _ROWMAT_SUB_       Row matrix subtraction expression template
    **********************************************************************************************/ 
    template<class RMAE, class T>
    _ROWMAT_SUB_< _ROWMAT_ATOMIC_EXPR_<RMAE>, RowMatrix<T> > 
        MATSUB(_ROWMAT_ATOMIC_EXPR_<RMAE>  const& expr, RowMatrix<T> const& mat)
    {return _ROWMAT_SUB_< _ROWMAT_ATOMIC_EXPR_<RMAE>, RowMatrix<T> > (expr, mat);} 




   /**********************************************************************************************
    * _ROWMAT_SUB_ MATSUB(RowMatrix const&, _ROWMAT_ATOMIC_EXPR_ const&)
    * 
    * MATSUB template specialization for RowMatrix inputs
    * 
    * Template Parameters:
    *     T    (implicit)     RowMatrix data type
    *     RMAE (implicit)     Row matrix atomic expression type data type
    * 
    * Parameters:
    *     mat                RowMatrix
    *     expr               Row matrix atomic expression
    * 
    * Returns:
    *     _ROWMAT_SUB_       Row matrix subtraction expression template
    **********************************************************************************************/ 
    template<class T, class RMAE>
    _ROWMAT_SUB_< RowMatrix<T>, _ROWMAT_ATOMIC_EXPR_<RMAE> > 
        MATSUB(RowMatrix<T> const& mat, _ROWMAT_ATOMIC_EXPR_<RMAE>  const& expr)
    {return _ROWMAT_SUB_< RowMatrix<T>, _ROWMAT_ATOMIC_EXPR_<RMAE> > (mat, expr);}  




   /**********************************************************************************************
    * _ROWMAT_SUB_ MATSUB(_ROWMAT_ATOMIC_EXPR_ const&, _ROWMAT_ATOMIC_EXPR_ const&)
    * 
    * MATSUB template specialization for RowMatrix inputs
    * 
    * Template Parameters:
    *     RMAE1 (implicit)     First row matrix atomic expression type data type
    *     RMAE2 (implicit)     Second row matrix atomic expression type data type
    * 
    * Parameters:
    *     expr1                Row matrix atomic expression
    *     expr2                Row matrix atomic expression
    * 
    * Returns:
    *     _ROWMAT_SUB_         Row matrix subtraction expression template
    **********************************************************************************************/ 
    template<class RMAE1, class RMAE2>
    _ROWMAT_SUB_< _ROWMAT_ATOMIC_EXPR_<RMAE1>, _ROWMAT_ATOMIC_EXPR_<RMAE2> > 
        MATSUB(_ROWMAT_ATOMIC_EXPR_<RMAE1>  const& expr1, _ROWMAT_ATOMIC_EXPR_<RMAE2>const& expr2)
    {
        return _ROWMAT_SUB_< _ROWMAT_ATOMIC_EXPR_<RMAE1>, 
                             _ROWMAT_ATOMIC_EXPR_<RMAE2> > (expr1, expr2);
    }   
    

    
}
#endif