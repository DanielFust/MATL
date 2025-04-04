/**************************************************************************************************
 * File: _matl_operator_minus_generic_.hpp
 * 
 * Description: Implementations for the (-) operator for constructing _EXPRESSION_ objects
 **************************************************************************************************/

namespace matl
{
   /***********************************************************************************************
    * operator- (_MATRIX_ const&, MATRIX_ const&)
    * 
    * Instantiate a _MATRIX_MINUS_MATRIX_ object from two _MATRIX_ objects
    * 
    * Template Parameters:
    *     MAT_T1: Matrix type
    *     MAT_T2: Matrix type
    * Returns:
    *    _MATRIX_MINUS_MATRIX_<_MATRIX_, _MATRIX_>    Matrix minus matrix expression
    ***********************************************************************************************/
    template<typename MAT_T1, typename MAT_T2>
    _MATRIX_MINUS_MATRIX_<_MATRIX_<MAT_T1>, _MATRIX_<MAT_T2> > 
    operator- (_MATRIX_<MAT_T1> const& mat1, 
               _MATRIX_<MAT_T2> const& mat2)
    {return _MATRIX_MINUS_MATRIX_<_MATRIX_<MAT_T1>, _MATRIX_<MAT_T2> >(mat1, mat2);} 
    
   /***********************************************************************************************
    * operator- (_MATRIX_ const&, EXPR const&)
    * 
    * Instantiate a _MATRIX_MINUS_MATRIX_ object from a _MATRIX_ and _EXPRESSION_ object
    * Template Parameters:
    *    MAT_T:  Matrix type
    *    EXPR: Expression type 
    * Returns:
    *   _MATRIX_MINUS_MATRIX_<_MATRIX_, EXPR>    Matrix minus expression object
    ***********************************************************************************************/
    template<typename MAT_T, typename EXPR>
    _MATRIX_MINUS_MATRIX_<_MATRIX_<MAT_T>, _EXPRESSION_<EXPR> > 
    operator- (_MATRIX_<MAT_T>      const& mat, 
               _EXPRESSION_<EXPR> const& expr)
    {return _MATRIX_MINUS_MATRIX_<_MATRIX_<MAT_T>, _EXPRESSION_<EXPR> >(mat, expr);}
    
   /***********************************************************************************************
    * operator- (EXPR const&, _MATRIX_ const&)
    * 
    * Instantiate a _MATRIX_MINUS_MATRIX_ object from an _EXPRESSION_ and _MATRIX_ object
    * 
    * Template Parameters:
    *     EXPR: Expression type
    *     MAT_T:  Matrix type
    * Returns:
    *     _MATRIX_MINUS_MATRIX_<EXPR, _MATRIX_>    Expression minus matrix object
    ***********************************************************************************************/
    template<typename EXPR, typename MAT_T>
    _MATRIX_MINUS_MATRIX_<_EXPRESSION_<EXPR>, _MATRIX_<MAT_T> > 
    operator- (_EXPRESSION_<EXPR> const& expr, 
                _MATRIX_<MAT_T>      const& mat)
    {return _MATRIX_MINUS_MATRIX_<_EXPRESSION_<EXPR>, _MATRIX_<MAT_T> >(expr, mat);}



   /***********************************************************************************************
    * operator- (EXPR const&, EXPR const&) 
    * 
    * Instantiate a _MATRIX_MINUS_MATRIX_ object from two _EXPRESSION_ objects
    * 
    * Template Parameters:
    *     EXPR1: Expression type
    *     EXPR2: Expression type
    * Returns:
    *     _MATRIX_MINUS_MATRIX_<EXPR1, EXPR2>    Expression minus expression object
    ***********************************************************************************************/
    template<typename EXPR1, typename EXPR2>
    _MATRIX_MINUS_MATRIX_<_EXPRESSION_<EXPR1>, _EXPRESSION_<EXPR2> >
    operator- (_EXPRESSION_<EXPR1> const& expr1, 
               _EXPRESSION_<EXPR2> const& expr2)
    {return _MATRIX_MINUS_MATRIX_<_EXPRESSION_<EXPR1>, _EXPRESSION_<EXPR2> >(expr1, expr2);} 


   /************************************************************************************************
    * operator- (_MATRIX_ const&, T)
    * 
    * Instantiate a _MATRIX_MINUS_SCALAR_ object from a _MATRIX_ and scalar
    * 
    * Template Parameters:
    *     MAT_T: Matrix type
    *     T:     Scalar type 
    * Returns:
    *     _MATRIX_MINUS_SCALAR_<_MATRIX_, T>    Matrix minus scalar object 
    ***********************************************************************************************/
    template<typename MAT_T, class T, enable_if_scalar<T> = true>
    _MATRIX_MINUS_SCALAR_<_MATRIX_<MAT_T>, T> operator- (_MATRIX_<MAT_T> const& mat, T scalar)
    {return _MATRIX_MINUS_SCALAR_<_MATRIX_<MAT_T>, T>(mat, scalar);}



   /************************************************************************************************
    * operator- (_EXPRESSION_ const&, T)
    *  
    * Instantiate a _MATRIX_MINUS_SCALAR_ object from an _EXPRESSION_ and scalar
    * 
    * Template Parameters:
    *     EXPR: Expression type
    *     T:     Scalar type
    * Returns:
    *    _MATRIX_MINUS_SCALAR_<EXPR, T>    Expression minus scalar object
    ***********************************************************************************************/
    template<typename EXPR, class T, enable_if_scalar<T> = true>  
    _MATRIX_MINUS_SCALAR_<_EXPRESSION_<EXPR>, T> 
    operator- (_EXPRESSION_<EXPR> const& expr, T scalar)
    {return _MATRIX_MINUS_SCALAR_<_EXPRESSION_<EXPR>, T>(expr, scalar);}


   /************************************************************************************************
    * operator- (T, _MATRIX_ const&)
    * 
    * Instantiate a _SCALAR_MINUS_MATRIX_ object from a scalar and _MATRIX_ object
    * 
    * Template Parameters:
    *     T:     Scalar type
    *     MAT_T: Matrix type
    * Returns:
    *     _SCALAR_MINUS_MATRIX_<MAT_T, T>    Scalar minus matrix object
    ***********************************************************************************************/
    template<class T, typename MAT_T, enable_if_scalar<T> = true>
    _SCALAR_MINUS_MATRIX_<T, _MATRIX_<MAT_T> > operator- (T scalar, _MATRIX_<MAT_T> const& mat)
    {return _SCALAR_MINUS_MATRIX_<T, _MATRIX_<MAT_T> >(scalar, mat);}

   /************************************************************************************************
    * operator- (T, EXPR const&)
    * 
    * Instantiate a _SCALAR_MINUS_MATRIX_ object from a scalar and _EXPRESSION_ object
    * 
    * Template Parameters:
    *     T:      Scalar type
    *     EXPR: Expression type
    * Returns:
    *     _SCALAR_MINUS_MATRIX_<T, EXPR>    Scalar minus expression object
    ***********************************************************************************************/
    template<class T, typename EXPR, enable_if_scalar<T> = true>
    _SCALAR_MINUS_MATRIX_<T, _EXPRESSION_<EXPR> > 
    operator- (T scalar, _EXPRESSION_<EXPR> const& expr)
    {return _SCALAR_MINUS_MATRIX_<T, _EXPRESSION_<EXPR> >(scalar, expr);}

}