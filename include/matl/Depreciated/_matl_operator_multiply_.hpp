#ifndef _MATL_OPERATOR_MULTIPLY__HPP
#define _MATL_OPERATOR_MULTIPLY__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains implementation of the operator(*) for the BLAS library
 **************************************************************************************************/
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"
#include "_matl_rowmat_scalar_prod_.hpp"
#include "_matl_generic_mat_scalar_prod_.hpp"

namespace matl
{



   /***********************************************************************************************
    *                           RowMatrix Scalar Multiplication
    ***********************************************************************************************/ 


    
   /***********************************************************************************************
    * _ROWMAT_SCALAR_PROD_ operator* (_ROWMAT_ATOMIC_EXPR_ const&, T);
    * 
    * Multiplication of Row Matrix Atomic Expression and scalar
    * 
    * Template Parameters:
    *     RMAE        (implicit) Row Matrix Atomic Expression type
    *     T           (implicit) Scalar data type
    * 
    * Parameters:
    *     expr         Matrix expression
    *     val          Scalar
    * 
    * Returns:
    *     _ROWMAT_SCALAR_PROD_ Row matrix scalar product
    ***********************************************************************************************/ 
    template<typename RMAE, class T, enable_if_scalar<T> = true>
    _ROWMAT_SCALAR_PROD_<RMAE,T> 
         operator* (_ROWMAT_ATOMIC_EXPR_<RMAE> const& expr, T val)  
    {return _ROWMAT_SCALAR_PROD_<RMAE,T>(downcast(expr),val);} 




   /***********************************************************************************************
    * _ROWMAT_SCALAR_PROD_ operator* (T, _ROWMAT_ATOMIC_EXPR_ const&);
    * 
    * Multiplication of Row Matrix Atomic Expression and scalar
    * 
    * Template Parameters:
    *     RMAE        (implicit) Row Matrix Atomic Expression type
    *     T           (implicit) Scalar data type
    * 
    * Parameters:
    *     expr         Matrix expression
    *     val          Scalar
    * 
    * Returns:
    *     _ROWMAT_SCALAR_PROD_ Row matrix scalar product
    ***********************************************************************************************/ 
    template<typename RMAE, class T, enable_if_scalar<T> = true>
    _ROWMAT_SCALAR_PROD_<RMAE,T> 
         operator* (T val,_ROWMAT_ATOMIC_EXPR_<RMAE> const& expr)  
    {return _ROWMAT_SCALAR_PROD_<RMAE,T>(downcast(expr),val);} 


  
   /***********************************************************************************************
    * _ROWMAT_SCALAR_PROD_ operator* (RowMatrix<T1> const&, T2);
    * 
    * Multiplication of Row Matrix Atomic Expression and scalar
    * 
    * Template Parameters:
    *     T1          (implicit) Row Matrix data type
    *     T2          (implicit) Scalar data type
    * 
    * Parameters:
    *     mat          Row Matrix
    *     val          Scalar
    * 
    * Returns:
    *     _ROWMAT_SCALAR_PROD_ Row matrix scalar product
    ***********************************************************************************************/ 
    template<typename T1, class T2, enable_if_scalar<T2> = true>
    _ROWMAT_SCALAR_PROD_<RowMatrix<T1>,T2> operator* (RowMatrix<T1> const& mat, T2 val)  
    {return _ROWMAT_SCALAR_PROD_<RowMatrix<T1>,T2> (mat,val);}




   /***********************************************************************************************
    * _ROWMAT_SCALAR_PROD_ operator* (T2, RowMatrix<T1> const&);
    * 
    * Multiplication of Row Matrix Atomic Expression and scalar
    * 
    * Template Parameters:
    *     T1          (implicit) Row Matrix data type
    *     T2          (implicit) Scalar data type
    * 
    * Parameters:
    *     val          Scalar
    *     mat          Row Matrix
    * 
    * Returns:
    *     _ROWMAT_SCALAR_PROD_ Row matrix scalar product
    ***********************************************************************************************/ 
    template<typename T1, class T2, enable_if_scalar<T2> = true>
    _ROWMAT_SCALAR_PROD_<RowMatrix<T1>,T2> operator* ( T2 val, RowMatrix<T1> const& mat)  
    {return _ROWMAT_SCALAR_PROD_<RowMatrix<T1>,T2>(mat,val);} 














   /***********************************************************************************************
    *                         Generic Matrix Scalar Multiplication
    ***********************************************************************************************/ 


    
   /***********************************************************************************************
    * _GENERIC_MAT_SCALAR _PROD_ operator* (_ATOMIC_EXPR_ const&, T);
    * 
    * Multiplication of a Generic Matrix Atomic Expression and scalar
    * 
    * Template Parameters:
    *     AE          (implicit) Matrix Atomic Expression type
    *     T           (implicit) Scalar data type
    * 
    * Parameters:
    *     expr         Matrix expression
    *     val          Scalar
    * 
    * Returns:
    *     _GENERIC_MAT_SCALAR_PROD_ Matrix scalar product
    ***********************************************************************************************/ 
    template<typename AE, class T, enable_if_scalar<T> = true>
    _GENERIC_MAT_SCALAR_PROD_<AE,T> 
         operator* (_ATOMIC_EXPR_<AE> const& expr, T val)  
    {return _GENERIC_MAT_SCALAR_PROD_<AE,T>(downcast(expr),val);} 





   /***********************************************************************************************
    * _GENERIC_MAT_SCALAR _PROD_ operator* (T, _ATOMIC_EXPR_ const&);
    * 
    * Multiplication of a Generic Matrix Atomic Expression and scalar
    * 
    * Template Parameters:
    *     T           (implicit) Scalar data type
    *     AE          (implicit) Matrix Atomic Expression type
    * 
    * Parameters:
    *     val          Scalar
    *     expr         Matrix expression
    * 
    * Returns:
    *     _GENERIC_MAT_SCALAR_PROD_ Matrix scalar product
    ***********************************************************************************************/ 
    template<typename AE, class T, enable_if_scalar<T> = true>
    _GENERIC_MAT_SCALAR_PROD_<AE,T> 
         operator* (T val, _ATOMIC_EXPR_<AE> const& expr)  
    {return _GENERIC_MAT_SCALAR_PROD_<AE,T>(downcast(expr),val);} 


  

   /***********************************************************************************************
    * _GENERIC_MAT_SCALAR _PROD_ operator* (Matrix const&, T);
    * 
    * Multiplication of a Generic Matrix and scalar
    * 
    * Template Parameters:
    *     T           (implicit) Scalar data type
    *     MAT         (implicit) Matrix type
    * 
    * Parameters:
    *     mat          Matrix
    *     val          Scalar
    * 
    * Returns:
    *     _GENERIC_MAT_SCALAR_PROD_ Matrix scalar product
    ***********************************************************************************************/ 
    template<typename MAT, class T, enable_if_scalar<T> = true>
    _GENERIC_MAT_SCALAR_PROD_<Matrix<MAT>,T> 
         operator* (Matrix<MAT> const& mat, T val)  
    {return _GENERIC_MAT_SCALAR_PROD_<Matrix<MAT>,T>(downcast(mat),val);}

   




  /***********************************************************************************************
    * _GENERIC_MAT_SCALAR _PROD_ operator* (T, Matrix const&);
    * 
    * Multiplication of a Generic Matrix and scalar
    * 
    * Template Parameters:
    *     T           (implicit) Scalar data type
    *     MAT         (implicit) Matrix type
    * 
    * Parameters:
    *     val          Scalar
    *     mat          Matrix
    * 
    * Returns:
    *     _GENERIC_MAT_SCALAR_PROD_ Matrix scalar product
    ***********************************************************************************************/ 
    template<typename MAT, class T, enable_if_scalar<T> = true>
    _GENERIC_MAT_SCALAR_PROD_<Matrix<MAT>,T> 
         operator* (T val, Matrix<MAT> const& mat)  
    {return _GENERIC_MAT_SCALAR_PROD_<Matrix<MAT>,T>(downcast(mat),val);}




















































   


}
#endif