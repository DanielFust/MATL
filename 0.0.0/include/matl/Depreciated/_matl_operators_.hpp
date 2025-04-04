#ifndef _BLAS_OPERATORS__HPP
#define _BLAS_OPERATORS__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains BLAS operators and bindings to BLAS operations
 **************************************************************************************************/
#include "_blas_forward_declr_.hpp"
#include "_blas_type_traits_.hpp"
#include "_blas_matsum_.hpp"

namespace blas
{

   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator+ (_MATEXPRESSION_ const&, _MATEXPRESSION_ const&);
    * 
    * General summation of matrix expressions
    * 
    * Template Parameters:
    *     ATOM_T          (implicit) Matrix atomic expression return type
    *     EXPR_T1         (implicit) First matrix expression
    *     EXPR_T2         (implicit) Second matrix expression
    * 
    * Parameters:
    *     expr1          First matrix expression
    *     expr2          Second matrix expression
    ***********************************************************************************************/ 
    template<typename ATOM_T, typename EXPR_T1, typename EXPR_T2>
    _ATOMIC_EXPR_<ATOM_T> operator+ (_MATEXPRESSION_<EXPR_T1> const& expr1, 
                                     _MATEXPRESSION_<EXPR_T2> const& expr2)
    {return MATSUM(expr1, expr2);}


   /***********************************************************************************************
    * _ATOMIC_EXPR_<ATOM_T> operator+ (Matrix<MAT_T> const&, _MATEXPRESSION_ const&);
    * 
    * General summation of matrix and matrix expression
    * 
    * Template Parameters:
    *     ATOM_T        (implicit) Matrix atomic expression return type
    *     MAT_T         (implicit) Matrix type
    *     EXPR_T        (implicit) Matrix expression type
    * 
    * Parameters:
    *     mat            Generic matrix
    *     expr           Matrix expression
    ***********************************************************************************************/ 
    template<typename ATOM_T, typename MAT_T, typename EXPR_T>
    _ATOMIC_EXPR_<ATOM_T> operator+ (Matrix<MAT_T>           const& mat, 
                                     _MATEXPRESSION_<EXPR_T> const& expr)
    {return MATSUM(mat, expr);}


   /***********************************************************************************************
    * _ATOMIC_EXPR_<ATOM_T> operator+ (_MATEXPRESSION_ const&, Matrix<MAT_T> const&);
    * 
    * General summation of matrix and matrix expression
    * 
    * Template Parameters:
    *     ATOM_T        (implicit) Matrix atomic expression return type
    *     EXPR_T        (implicit) Matrix expression type
    *     MAT_T         (implicit) Matrix type
    * 
    * Parameters:
    *     expr           Matrix expression
    *     mat            Generic matrix
    ***********************************************************************************************/ 
    template<typename ATOM_T, typename EXPR_T, typename MAT_T>
    _ATOMIC_EXPR_<ATOM_T> operator+ (_MATEXPRESSION_<EXPR_T> const& expr,
                                     Matrix<MAT_T>           const& mat)
    {return MATSUM(expr,mat);}


   /***********************************************************************************************
    * _ATOMIC_EXPR_<ATOM_T> operator+ (_MATEXPRESSION_ const&, Matrix<MAT_T> const&);
    * 
    * General summation of matrix and matrix expression
    * 
    * Template Parameters:
    *     ATOM_T        (implicit) Matrix atomic expression return type
    *     MAT_T1        (implicit) First matrix type
    *     MAT_T2        (implicit) Second matrix type
    * 
    * Parameters:
    *     mat1           First matrix
    *     mat2           Second matrix
    ***********************************************************************************************/ 
    template<typename ATOM_T, typename MAT_T1, typename MAT_T2>
    _ATOMIC_EXPR_<ATOM_T> operator+ (Matrix<MAT_T1>           const& mat1,
                                     Matrix<MAT_T2>           const& mat2)
    {return MATSUM(mat1,mat2);}  

}

#endif




