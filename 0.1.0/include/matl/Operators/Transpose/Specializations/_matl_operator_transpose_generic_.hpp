
/**************************************************************************************************
 * Tranpose operations for Generic Matrices and Matrix Expressions
 **************************************************************************************************/


namespace matl
{
   /***********************************************************************************************
    * _GENERIC_MAT_TRANSPOSE_ transpose(Matrix const& mat) 
    * 
    * General transpose of generic matrix
    * 
    * Template Parameters:
    *     MAT_T              (implicit) matrix type
    * 
    * Parameters:
    *     mat                Matrix
    * 
    * Returns:
    *     _GENERIC_MAT_TRANSPOSE_  Generic Matrix Transpose
    ***********************************************************************************************/ 
    template<typename MAT_T>
    auto transpose(Matrix<MAT_T> const& mat)  
    {return _GENERIC_MAT_TRANSPOSE_(mat);}

   /***********************************************************************************************
    * _GENERIC_MAT_TRANSPOSE_ transpose(_MATEXPRESSION_ const&)
    * 
    * General transpose of generic matrix expression
    * 
    * Template Parameters:
    *     EXPR_T              (implicit) matrix type
    * 
    * Parameters:
    *     expr                Matrix expression
    * 
    * Returns:
    *     _GENERIC_MAT_TRANSPOSE_  Generic Matrix Transpose
    ***********************************************************************************************/ 
    template<typename EXPR_T>
    auto transpose(_MATEXPRESSION_<EXPR_T> const& expr)  
    {return _GENERIC_MAT_TRANSPOSE_(expr);} 


}
