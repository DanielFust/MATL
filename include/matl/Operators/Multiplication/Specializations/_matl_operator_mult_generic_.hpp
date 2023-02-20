
/**************************************************************************************************
 * Multiplication operations for Generic Matrices and Matrix Expressions
 **************************************************************************************************/


//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"

namespace matl
{

   /***********************************************************************************************
    * _MATEXPRESSION_ operator* (_MATEXPRESSION_ const&, Matrix const&);
    * 
    * General multiplication of matrix and matrix expression
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
    *     _MATEXPRESSION_    Matrix - Matrix Product expression (specializations possible)
    ***********************************************************************************************/ 
    template<typename EXPR_T, typename MAT_T>
    auto operator* (_MATEXPRESSION_<EXPR_T> const& expr,
                    Matrix<MAT_T>           const& mat)  
    {
        //std::cout << "'operator* (_MATEXPRESSION_ const& expr, Matrix const& mat)' called" << std::endl;
        typename TypeTraits<EXPR_T>::matrix_type  mat2;
        mat2.resize(expr.rows(), expr.cols());
        expr.eval(mat2);
        //std::cout << "mat2=" << std::endl;
        //mat2.println();
        //auto mat2 = expr.eval();
        return mat2*mat;
        //mat2.println();
        //return expr.eval()*mat;
    }


   /***********************************************************************************************
    * _MATEXPRESSION_ operator* (Matrix const&, _MATEXPRESSION_ const&);
    * 
    * General multiplication of matrix and matrix expression
    * 
    * Template Parameters:
    *     EXPR_T        (implicit) Expression type
    *     MAT_T         (implicit) matrix type
    * 
    * Parameters:
    *     mat            Matrix
    *     expr           Matrix expression
    * 
    * Returns:
    *     _MATEXPRESSION_    Matrix - Matrix Product expression (specializations possible)
    ***********************************************************************************************/ 
    template<typename EXPR_T, typename MAT_T>
    auto operator* (Matrix<MAT_T>           const& mat,
                    _MATEXPRESSION_<EXPR_T> const& expr)  
    {
        //std::cout << "'operator* (Matrix const& mat, _MATEXPRESSION_ const& expr)' called" << std::endl;
        typename TypeTraits<EXPR_T>::matrix_type  mat2;
        mat2.resize(expr.rows(), expr.cols());
        expr.eval(mat2);
        //std::cout << "mat2=" << std::endl;
        //mat2.println();
        return mat*mat2;
    }


   /***********************************************************************************************
    * _MATEXPRESSION_ operator* (_MATEXPRESSION_ const&, _MATEXPRESSION_ const&);
    * 
    * General multiplication of matrix expressions
    * 
    * Template Parameters:
    *     EXPR_T1       (implicit) First expression type
    *     EXPR_T2       (implicit) First expression type
    * 
    * Parameters:
    *     expr1          First matrix expression
    *     expr2          Second matrix expression
    * 
    * Returns:
    *     _MATEXPRESSION_    Matrix - Matrix Product expression (specializations possible)
    ***********************************************************************************************/ 
    template<typename EXPR_T1, typename EXPR_T2>
    auto operator* (_MATEXPRESSION_<EXPR_T1> const& expr1,
                    _MATEXPRESSION_<EXPR_T2> const& expr2)  
    {
        //std::cout << "'operator* (_MATEXPRESSION_ const& expr1, _MATEXPRESSION_ const& expr2)' called" << std::endl;
        typename TypeTraits<EXPR_T1>::matrix_type mat1;
        typename TypeTraits<EXPR_T2>::matrix_type mat2;
    
        mat1.resize(expr1.rows(), expr1.cols());
        mat2.resize(expr2.rows(), expr2.cols());
        
        expr1.eval(mat1);
        expr2.eval(mat2);

        return mat1 * mat2;
    }


   /***********************************************************************************************
    * _GENERIC_MAT_PROD_ operator* (Matrix const&, Matrix const&);
    * 
    * General multiplication of matrices
    * 
    * Template Parameters:
    *     MAT_T1         (implicit) First matrix type
    *     MAT_T2         (implicit) Second matrix type
    * 
    * Parameters:
    *     mat1           First matrix
    *     mat2           Second matrix
    * 
    * Returns:
    *     _GENERIC_MAT_PROD_    Matrix - Matrix Product expression 
    ***********************************************************************************************/ 
    template<typename MAT_T1, typename MAT_T2>
    auto operator* (Matrix<MAT_T1>  const& mat1,
                    Matrix<MAT_T2>  const& mat2)  
    {
        //std::cout << "'operator* (Matrix const& mat1, Matrix const& mat2)' called" << std::endl;
        
        return _GENERIC_MAT_PROD_(mat1,mat2);
    }

















   /***********************************************************************************************
    * _GENERIC_MAT_PROD_ operator* (_ATOMIC_EXPR_ const&, Matrix const&);
    * 
    * General multiplication of matrix and atomic matrix expression
    * 
    * Template Parameters:
    *     AE_T                 (implicit) Atomic expression type
    *     MAT_T                (implicit) Matrix type
    * 
    * Parameters:
    *     expr                 Atomic matrix expression
    *     mat                  Matrix
    * 
    * Returns:
    *     _GENERIC_MAT_PROD_    Matrix - Matrix Product expression 
    ***********************************************************************************************/ 
    template<typename AE_T, typename MAT_T>
    auto operator* (_ATOMIC_EXPR_<AE_T> const& expr,
                    Matrix<MAT_T>       const& mat)  
    {return _GENERIC_MAT_PROD_(expr,mat);}






   /***********************************************************************************************
    * _GENERIC_MAT_PROD_ operator* (Matrix const&, _ATOMIC_EXPR_ const&);
    * 
    * General multiplication of matrix and atomic matrix expression
    * 
    * Template Parameters:
    *     MAT_T                (implicit) Matrix type
    *     AE_T                 (implicit) Atomic expression type
    * 
    * Parameters:
    *     mat                  Matrix
    *     expr                 Atomic matrix expression
    * 
    * Returns:
    *     _GENERIC_MAT_PROD_    Matrix - Matrix Product expression 
    ***********************************************************************************************/ 
    template<typename AE_T, typename MAT_T>
    auto operator* (Matrix<MAT_T>       const& mat,
                    _ATOMIC_EXPR_<AE_T> const& expr)  
    {return _GENERIC_MAT_PROD_(mat,expr);}






   /***********************************************************************************************
    * _GENERIC_MAT_SUM_ operator* (_ATOMIC_EXPR_ const&, _ATOMIC_EXPR_ const&);
    * 
    * General mutiplication of atomic matrix expressions
    * 
    * Template Parameters:
    *     AE_T1                (implicit) First atomic expression type
    *     AE_T2                (implicit) Second atomic expression type
    * 
    * Parameters:
    *     expr1                First atomic matrix expression
    *     expr2                Second atomic matrix expression
    * 
    * Returns:
    *     _GENERIC_MAT_PROD_    Matrix - Matrix Product expression 
    ***********************************************************************************************/ 
    template<typename AE_T1, typename AE_T2>
    auto operator* (_ATOMIC_EXPR_<AE_T1> const& expr1,
                    _ATOMIC_EXPR_<AE_T2> const& expr2)  
    {return _GENERIC_MAT_PROD_(expr1,expr2);}








   /***********************************************************************************************
    * _GENERIC_MAT_SCALAR_PROD_ operator* (Matrix const&, T);
    * 
    * General multiplication of matrix and scalar
    * 
    * Template Parameters:
    *     MAT_T                        Matrix type
    *     T                            Scalar data type
    * 
    * Parameters:
    *     mat                          Matrix
    *     val                          Scalar
    * 
    * Returns:
    *     _GENERIC_MAT_SCALAR_PROD_    Generic matrix-scalar multiplication
    ***********************************************************************************************/ 
    template<typename MAT_T, class T, enable_if_scalar<T> = true>
    auto operator* (Matrix<MAT_T> const& mat,
                    T                    val)  
    {
        //std::cout << "mat.rows()="<<mat.rows()<<", mat.cols()="<<mat.cols()<<std::endl;
        return _GENERIC_MAT_SCALAR_PROD_(mat,val);
    }





   /***********************************************************************************************
    * _GENERIC_MAT_SCALAR_PROD_ operator* (T, Matrix const&);
    * 
    * General multiplication of a matrix and scalar
    * 
    * Template Parameters:
    *     MAT_T                        Matrix type
    *     T                            Scalar data type
    * 
    * Parameters:
    *     val                          Scalar
    *     mat                          Matrix
    * 
    * Returns:
    *     _GENERIC_MAT_SCALAR_PROD_    Generic matrix-scalar multiplication
    ***********************************************************************************************/ 
    template<typename MAT_T, class T, enable_if_scalar<T> = true>
    auto operator* (T                    val,
                    Matrix<MAT_T> const& mat)  
    {return _GENERIC_MAT_SCALAR_PROD_(mat,val);}     









   /***********************************************************************************************
    * _GENERIC_MAT_SCALAR_PROD_ operator* (_ATOMIC_EXPR_ const&, T);
    * 
    * General multiplication of an atomic matrix expression and a scalar
    * 
    * Template Parameters:
    *     AE_T                        Atomic expression type
    *     T                           Scalar data type
    * 
    * Parameters:
    *     expr                        Atomic matrix expression
    *     val                         Scalar
    * 
    * Returns:
    *     _GENERIC_MAT_SCALAR_PROD_    Generic matrix-scalar multiplication
    ***********************************************************************************************/ 
    template<typename AE_T, class T, enable_if_scalar<T> = true>
    auto operator* (_ATOMIC_EXPR_<AE_T> const& expr,
                    T                          val)  
    {return _GENERIC_MAT_SCALAR_PROD_(expr,val);}








   /***********************************************************************************************
    * _GENERIC_MAT_SCALAR_PROD_ operator* (T, _ATOMIC_EXPR_ const&);
    * 
    * General multiplication of an atomic matrix expression and scalar
    * 
    * Template Parameters:
    *     AE_T                        Atomic expression type
    *     T                           Scalar data type
    * 
    * Parameters:
    *     val                         Scalar
    *     expr                        Atomic matrix expression
    * 
    * Returns:
    *     _GENERIC_MAT_SCALAR_PROD_    Generic matrix-scalar multiplication
    ***********************************************************************************************/ 
    template<typename AE_T, class T, enable_if_scalar<T> = true>
    auto operator* (T                          val,
                    _ATOMIC_EXPR_<AE_T> const& expr)  
    {return _GENERIC_MAT_SCALAR_PROD_(expr,val);}













   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator* (_MATEXPRESSION_ const&, T);
    * 
    * General multiplication of matrix expression and scalar
    * 
    * Template Parameters:
    *     EXPR_T                      Expression type
    *     T                           Scalar data type
    * 
    * Parameters:
    *     expr                        Expression
    *     val                         Scalar
    * 
    * Returns:
    *     _ATOMIC_EXPR_    Generic matrix-scalar multiplication
    ***********************************************************************************************/ 
    template<typename EXPR_T, class T, enable_if_scalar<T> = true>
    auto operator* (_MATEXPRESSION_<EXPR_T> const& expr,
                    T                              val)  
    {
        typename TypeTraits<EXPR_T>::matrix_type mat;
        mat.resize(expr.rows(), expr.cols());
        expr.eval(mat);
        
        return mat*val;
    }








   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator* (T, _MATEXPRESSION_ const&);
    * 
    * General multiplication of matrix expression and scalar
    * 
    * Template Parameters:
    *     EXPR_T                      Expression type
    *     T                           Scalar data type
    * 
    * Parameters:
    *     val                         Scalar
    *     expr                        Expression
    * 
    * Returns:
    *     _ATOMIC_EXPR_               Generic matrix-scalar multiplication
    ***********************************************************************************************/ 
    template<typename EXPR_T, class T, enable_if_scalar<T> = true>
    auto operator* (T                              val,
                    _MATEXPRESSION_<EXPR_T> const& expr)  
    {
        typename TypeTraits<EXPR_T>::matrix_type mat;
        mat.resize(expr.rows(), expr.cols());
        expr.eval(mat);
        return mat*val;
    }  







}
//#endif