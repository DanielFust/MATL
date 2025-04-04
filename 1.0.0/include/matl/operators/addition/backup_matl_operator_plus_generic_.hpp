/**************************************************************************************************
 * _matl_operator_plus_generic_.hpp
 * 
 * This file contains implementations for the (+) operator for constructing _EXPRESSION_ objects
 * with the _MATRIX_ and _EXPRESSION_ base classes. 
 **************************************************************************************************/

namespace matl
{
   /***********************************************************************************************
    * operator+ (_MATRIX_ const&, MATRIX_ const&)
    * 
    * Instantiate a _MATRIX_PLUS_MATRIX_ object from two _MATRIX_ objects
    * 
    * Template Parameters:
    *     MAT_T1: Matrix type
    *     MAT_T2: Matrix type
    * Returns:
    *    _MATRIX_PLUS_MATRIX_<_MATRIX_, _MATRIX_>    Matrix plus matrix expression
    ***********************************************************************************************/
    template<typename MAT_T1, typename MAT_T2>
    _MATRIX_PLUS_MATRIX_<_MATRIX_<MAT_T1>, _MATRIX_<MAT_T2> > 
    operator+ (_MATRIX_<MAT_T1> const& mat1, 
               _MATRIX_<MAT_T2> const& mat2)
    {return _MATRIX_PLUS_MATRIX_<_MATRIX_<MAT_T1>, _MATRIX_<MAT_T2> >(mat1, mat2);} 

   /***********************************************************************************************
    * operator+ (_MATRIX_ const&, EXPR_T const&)
    * 
    * Instantiate a _MATRIX_PLUS_MATRIX_ object from a _MATRIX_ and _EXPRESSION_ object
    * Template Parameters:
    *    MAT_T:  Matrix type
    *    EXPR_T: Expression type 
    * Returns:
    *   _MATRIX_PLUS_MATRIX_<_MATRIX_, EXPR_T>    Matrix plus expression object
    ***********************************************************************************************/
    template<typename MAT_T, typename EXPR_T>
    _MATRIX_PLUS_MATRIX_<_MATRIX_<MAT_T>, _EXPRESSION_<EXPR_T> > operator+ (_MATRIX_<MAT_T>      const& mat, 
                                                                            _EXPRESSION_<EXPR_T> const& expr)
    {return _MATRIX_PLUS_MATRIX_<_MATRIX_<MAT_T>, _EXPRESSION_<EXPR_T> >(mat, expr);}

   /************************************************************************************************
    * operator+ (EXPR_T const&, _MATRIX_ const&)
    * 
    * Instantiate a _MATRIX_PLUS_MATRIX_ object from an _EXPRESSION_ and _MATRIX_ object
    * 
    * Template Parameters:
    *    EXPR_T: Expression type
    *    MAT_T:  Matrix type
    * Returns:
    *   _MATRIX_PLUS_MATRIX_<EXPR_T, _MATRIX_>    Expression plus matrix object
    ************************************************************************************************/
    template<typename EXPR_T, typename MAT_T>
    _MATRIX_PLUS_MATRIX_<_EXPRESSION_<EXPR_T>, _MATRIX_<MAT_T> > operator+ (_EXPRESSION_<EXPR_T> const& expr, 
                                                                            _MATRIX_<MAT_T>      const& mat)
    {return _MATRIX_PLUS_MATRIX_<_EXPRESSION_<EXPR_T>, _MATRIX_<MAT_T> >(expr, mat);}
    
   /************************************************************************************************
    * operator+ (EXPR_T const&, EXPR_T const&)
    * 
    * Instantiate a _MATRIX_PLUS_MATRIX_ object from two _EXPRESSION_ objects
    * 
    * Template Parameters:
    *    EXPR_T1: Expression type
    *    EXPR_T2: Expression type
    * Returns:
    *  _MATRIX_PLUS_MATRIX_<EXPR_T1, EXPR_T2>    Expression plus expression object
    ************************************************************************************************/  
    template<typename EXPR_T1, typename EXPR_T2>
    _MATRIX_PLUS_MATRIX_<_EXPRESSION_<EXPR_T1>, _EXPRESSION_<EXPR_T2> > operator+ (_EXPRESSION_<EXPR_T1> const& expr1, 
                                                                                   _EXPRESSION_<EXPR_T2> const& expr2)
    {return _MATRIX_PLUS_MATRIX_<_EXPRESSION_<EXPR_T1>, _EXPRESSION_<EXPR_T2> >(expr1, expr2);}

   /***********************************************************************************************
    * operator+ (_MATRIX_ const&, T const&)
    * 
    * Instantiate a _MATRIX_PLUS_SCALAR_ object from a _MATRIX_ and scalar
    * 
    * Template Parameters:
    *   MAT_T: Matrix type
    *   T:     Scalar type
    * Returns:
    *  _MATRIX_PLUS_SCALAR_<_MATRIX_, T>    Matrix plus scalar object
    ***********************************************************************************************/
    template<typename MAT_T, class T, enable_if_scalar<T> = true>
    _MATRIX_PLUS_SCALAR_<_MATRIX_<MAT_T>, T> operator+ (_MATRIX_<MAT_T> const& mat, T scalar)
    {return _MATRIX_PLUS_SCALAR_<_MATRIX_<MAT_T>, T>(mat, scalar);}

   /***********************************************************************************************
    * operator+ (T const&, _MATRIX_ const&)
    * 
    * Instantiate a _MATRIX_PLUS_SCALAR_ object from a scalar and _MATRIX_ object
    * 
    * Template Parameters:
    *     T:     Scalar type
    *     MAT_T: Matrix type
    * Returns:
    *   _MATRIX_PLUS_SCALAR_<T, _MATRIX_>    Scalar plus matrix object
    ***********************************************************************************************/
    template<class T, typename MAT_T, enable_if_scalar<T> = true>
    _MATRIX_PLUS_SCALAR_<_MATRIX_<MAT_T>,T> operator+ (T scalar, _MATRIX_<MAT_T> const& mat)
    {return _MATRIX_PLUS_SCALAR_<_MATRIX_<MAT_T>,T>(mat, scalar);}

   /***********************************************************************************************
    * operator+ (EXPR_T const&, T const&)
    * 
    * Instantiate a _MATRIX_PLUS_SCALAR_ object from an _EXPRESSION_ and scalar
    * 
    * Template Parameters:
    *     EXPR_T: Expression type
    *     T:      Scalar type
    * Returns:
    *     _MATRIX_PLUS_SCALAR_<EXPR_T, T>    Expression plus scalar object
    ************************************************************************************************/
    template<typename EXPR_T, class T, enable_if_scalar<T> = true>
    _MATRIX_PLUS_SCALAR_<_EXPRESSION_<EXPR_T>, T> operator+ (_EXPRESSION_<EXPR_T> const& expr, T scalar)
    {return _MATRIX_PLUS_SCALAR_<_EXPRESSION_<EXPR_T>, T>(expr, scalar);}
    
    
   /************************************************************************************************
    * operator+ (T const&, EXPR_T const&)
    * 
    * Instantiate a _MATRIX_PLUS_SCALAR_ object from a scalar and _EXPRESSION_ object
    * 
    * Template Parameters:
    *     T:      Scalar type
    *     EXPR_T: Expression type
    * Returns:
    *     _MATRIX_PLUS_SCALAR_<T, EXPR_T>    Scalar plus expression object
    ************************************************************************************************/
    template<class T, typename EXPR_T, enable_if_scalar<T> = true>  
    _MATRIX_PLUS_SCALAR_<_EXPRESSION_<EXPR_T>, T> operator+ (T scalar, _EXPRESSION_<EXPR_T> const& expr)
    {return _MATRIX_PLUS_SCALAR_<_EXPRESSION_<EXPR_T>, T>(expr, scalar);}
}   





