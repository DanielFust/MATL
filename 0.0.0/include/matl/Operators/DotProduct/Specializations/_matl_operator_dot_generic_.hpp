
/**************************************************************************************************
 * Dot product operations for Generic Matrices and Matrix Expressions
 **************************************************************************************************/


namespace matl
{
   /***********************************************************************************************
    * dtype dot(Matrix const&, Matrix const&) 
    * 
    * General transpose of generic matrix
    * 
    * Template Parameters:
    *     MAT1               (implicit) matrix type
    *     MAT2               (implicit) matrix type
    * 
    * Parameters:
    *     mat1                Matrix
    *     mat2                Matrix
    * 
    * Returns:
    *     dtype      scalar product result
    ***********************************************************************************************/ 
    template<typename MAT1, typename MAT2>
    auto dot(Matrix<MAT1> const& mat1, Matrix<MAT2> const& mat2)  
    {
        typename promote<typename TypeTraits<MAT1>::dtype, 
                         typename TypeTraits<MAT2>::dtype>::type val = 0;

        #ifdef MATL_DEBUG
        if ((mat1.rows()!=mat2.rows())||(mat1.cols() != mat2.cols()))
        {
            std::cerr << "Error in 'dtype dot(Matrix const&, Matrix const&)'\n" 
                      << "Matrix with dimensions " << mat1.rows() << "x" << mat1.cols() 
                      << " cannot be dotted with Matrix with dimensions "
                      << mat2.rows() << "x" << mat2.cols() << "."
                      << std::endl; 
            abort();          
        }
        #endif

        for (index_t i=1; i<=mat1.rows(); ++i)
            for (index_t j=1; j<=mat1.cols(); ++j)
                val += mat1(i,j)*mat2(i,j);
        return val;        
    }




   /***********************************************************************************************
    * dtype dot(Matrix const&, _MATEXPRESSION_ const&) 
    * 
    * General transpose of generic matrix
    * 
    * Template Parameters:
    *     MAT                (implicit) matrix type
    *     EXPR               (implicit) matrix expression type
    * 
    * Parameters:
    *     mat                 Matrix
    *     expr                Matrix expression
    * 
    * Returns:
    *     dtype      scalar product result
    ***********************************************************************************************/ 
    template<typename MAT, typename EXPR>
    auto dot(Matrix<MAT> const& mat, _MATEXPRESSION_<EXPR> const& expr)  
    {
        typename promote<typename TypeTraits<MAT>::dtype, 
                         typename TypeTraits<EXPR>::dtype>::type val = 0;

        #ifdef MATL_DEBUG
        if ((mat.rows()!=expr.rows())||(mat.cols() != expr.cols()))
        {
            std::cerr << "Error in 'dtype dot(Matrix const&, _MATEXPRESSION_ const&)'\n" 
                      << "Matrix with dimensions " << mat.rows() << "x" << mat.cols() 
                      << " cannot be dotted with Matrix Expression with dimensions "
                      << expr.rows() << "x" << expr.cols() << "."
                      << std::endl; 
            abort();          
        }
        #endif

        for (index_t i=1; i<=mat.rows(); ++i)
            for (index_t j=1; j<=mat.cols(); ++j)
                val += mat(i,j)*expr(i,j);
        return val;        
    }



   /***********************************************************************************************
    * dtype dot(_MATEXPRESSION_ const&, Matrix const&) 
    * 
    * General Matrix-Matrix Expression dot product
    * 
    * Template Parameters:
    *     EXPR               (implicit) matrix expression type
    *     MAT                (implicit) matrix type
    * 
    * Parameters:
    *     expr                Matrix expression
    *     mat                 Matrix
    * 
    * Returns:
    *     dtype      scalar product result
    ***********************************************************************************************/ 
    template<typename EXPR, typename MAT>
    auto dot(_MATEXPRESSION_<EXPR> const& expr, Matrix<MAT> const& mat)  
    {
        typename promote<typename TypeTraits<MAT>::dtype, 
                         typename TypeTraits<EXPR>::dtype>::type val = 0;

        #ifdef MATL_DEBUG
        if ((mat.rows()!=expr.rows())||(mat.cols() != expr.cols()))
        {
            std::cerr << "Error in 'dtype dot(_MATEXPRESSION_ const&, Matrix const&)'\n" 
                      << "Matrix with dimensions " << mat.rows() << "x" << mat.cols() 
                      << " cannot be dotted with Matrix Expression with dimensions "
                      << expr.rows() << "x" << expr.cols() << "."
                      << std::endl; 
            abort();          
        }
        #endif

        for (index_t i=1; i<=mat.rows(); ++i)
            for (index_t j=1; j<=mat.cols(); ++j)
                val += mat(i,j)*expr(i,j);
        return val;        
    }




   /***********************************************************************************************
    * dtype dot(_MATEXPRESSION_ const&, _MATEXPRESSION_ const&) 
    * 
    * General Matrix-Matrix Expression dot product
    * 
    * Template Parameters:
    *     EXPR1              (implicit) matrix expression type
    *     EXPR2              (implicit) matrix expression type
    * 
    * Parameters:
    *     expr1               Matrix expression
    *     expr1               Matrix expression
    * 
    * Returns:
    *     dtype      scalar product result
    ***********************************************************************************************/ 
    template<typename EXPR1, typename EXPR2>
    auto dot(_MATEXPRESSION_<EXPR1> const& expr1, _MATEXPRESSION_<EXPR2> const& expr2)  
    {
        typename promote<typename TypeTraits<EXPR1>::dtype, 
                         typename TypeTraits<EXPR2>::dtype>::type val = 0;

        #ifdef MATL_DEBUG
        if ((expr1.rows()!=expr2.rows())||(expr1.cols() != expr2.cols()))
        {
            std::cerr << "Error in 'dtype dot(_MATEXPRESSION_ const&, _MATEXPRESSION_ const&)'\n" 
                      << "Matrix Expression with dimensions " 
                      << expr1.rows() << "x" << expr1.cols() 
                      << " cannot be dotted with Matrix Expression with dimensions "
                      << expr2.rows() << "x" << expr2.cols() << "."
                      << std::endl; 
            abort();          
        }
        #endif

        for (index_t i=1; i<=expr1.rows(); ++i)
            for (index_t j=1; j<=expr1.cols(); ++j)
                val += expr1(i,j)*expr2(i,j);
        return val;        
    }   



}
