#ifndef _MATL_OPERATOR_MINUS__HPP
#define _MATL_OPERATOR_MINUS__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains implementation of the operator(+) for the BLAS library
 **************************************************************************************************/
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"
#include "_matl_matsub_.hpp"

namespace matl
{

   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator- (Matrix const&, _MATEXPRESSION_ const&);
    * 
    * General subtraction of matrix and matrix expression
    * 
    * Template Parameters:
    *     EXPR_T        (implicit) Expression type
    *     MAT_T         (implicit) matrix type
    * 
    * Parameters:
    *     expr           Matrix expression
    *     mat            Matrix
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename EXPR_T, typename MAT_T>
    auto operator- (_MATEXPRESSION_<EXPR_T> const& expr,
                    Matrix<MAT_T>           const& mat)  
    {return MATSUB(downcast(expr),downcast(mat));}



   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator- (_MATEXPRESSION_ const&, Matrix const&,);
    * 
    * General subtraction of matrix expression and matrix
    * 
    * Template Parameters:
    *     MAT_T         (implicit) matrix type
    *     EXPR_T        (implicit) Expression type
    * 
    * Parameters:
    *     mat            Matrix
    *     expr           Matrix expression
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename MAT_T, typename EXPR_T>
    auto operator- (Matrix<MAT_T>           const& mat,
                    _MATEXPRESSION_<EXPR_T> const& expr)
    {return MATSUB(downcast(mat),downcast(expr));}



   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator- (Matrix const&, Matrix const&,);
    * 
    * General subtraction of matrices
    * 
    * Template Parameters:
    *     MAT_T1         (implicit) first matrix type
    *     MAT_T2         (implicit) second matrix type
    * 
    * Parameters:
    *     mat1            First matrix
    *     mat2            Second matrix
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename MAT_T1, typename MAT_T2>
    auto operator- (Matrix<MAT_T1>           const& mat1,
                    Matrix<MAT_T2>           const& mat2)
    {return MATSUB(downcast(mat1),downcast(mat2));}  







   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator- (_MATEXPRESSION_ const&, const&, _MATEXPRESSION_ const&);
    * 
    * General subtraction of expressions
    * 
    * Template Parameters:
    *     EXPR_T1       (implicit) First expression type
    *     EXPR_T2       (implicit) Second expression type
    * 
    * Parameters:
    *     expr1          First matrix expression
    *     expr2          Second matrix expression
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename EXPR_T1, typename EXPR_T2>
    decltype(auto) operator- (_MATEXPRESSION_<EXPR_T1> const& expr1,
                              _MATEXPRESSION_<EXPR_T2> const& expr2) 
    {return MATSUB(downcast(expr1), downcast(expr2));} 



}
#endif