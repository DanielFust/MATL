#ifndef _BLAS_ROWMAT_PROD__HPP
#define _BLAS_ROWMAT_PROD__HPP

#include "_blas_forward_declr_.hpp"
#include "_blas_type_traits_.hpp"

namespace blas
{
   /**********************************************************************************************
    * class _ROWMAT_PROD_ : public _MATEXPRESSION_
    * 
    * Expression template for matrix products between two dense row-major matricies (RowMatrix).
    * 
    * Template Parameters:
    *     T1 (implicit)        Data type of first RowMatrix
    *     T2 (implicit)        Data type of second RowMatrix     
    **********************************************************************************************/
    template<class T1, class T2>
    class _ROWMAT_PROD_ : public _MATEXPRESSION_< _ROWMAT_PROD_<T1,T2> >
    {
      private:
        RowMatrix<T1> const& mat1_;
        RowMatrix<T2> const& mat2_;
      
      protected:

        // Return Type
        typedef typename promote<T1, T2>::type dtype;

        // Constructor
        _ROWMAT_PROD_(RowMatrix<T1> const& mat1, RowMatrix<T2> const& mat2);

        // Evaluate expression 
        RowMatrix<dtype> eval() const;

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(Matrix<MAT_T>* target) const;

        // Number of rows in resultant matrix
        index_t rows() const {return mat1_.rows();}
        
        // Number of columns in resultant matrix
        index_t cols() const {return mat2_.cols();}

        // Total size of resultant matrix
        index_t size() const {return this->rows() * this->cols();}

    };




   /**********************************************************************************************
    * _ROWMAT_PROD_(RowMatrix const& mat1, RowMatrix const& mat2)
    * 
    * Constructor for RowMatrix-RowMatrix multiplication expression templates. 
    * 
    * Template Parameters:
    *     T1 (implicit)    First RowMatrix data type
    *     T2 (implicit)    Second RowMatrix data type
    * 
    * Parameters:
    *     mat1    (RowMatrix const&) Left hand side matrix
    *     mat2    (RowMatrix const&) Right hand side matrix
    * 
    * Returns:
    *     this    (_ROWMAT_PROD_) Expression template for row major matrix multiplication
    **********************************************************************************************/
    template<class T1, class T2>
    _ROWMAT_PROD_<T1,T2>::_ROWMAT_PROD_(RowMatrix<T1> const& mat1, 
                                        RowMatrix<T2> const& mat2) : mat1_(mat1), 
                                                                     mat2_(mat2) 
    {
        // Check for dimensional compatibility
        #ifdef DIM_CHECK 
        if ((mat1.rows() <= 0) || (mat1.cols() <= 0))
        {
            std::cerr << "Error in _ROWMAT_PROD_(RowMatrix<T1> const& mat1, " 
                      << "RowMatrix<T2> const& mat2)"
                      << "\n    RowMatrix 'mat1' with dimensions " 
                      << mat1.rows() << "x" << mat1.cols()
                      << " unsuitable for matrix multiplication."
                      << std::endl;
            abort();          
        }

        if ((mat2.rows() <= 0) || (mat2.cols() <= 0))
        {
            std::cerr << "Error in _ROWMAT_PROD_(RowMatrix<T1> const& mat1, " 
                      << "RowMatrix<T2> const& mat2)"
                      << "\n    RowMatrix 'mat2' with dimensions " 
                      << mat2.rows() << "x" << mat2.cols()
                      << " unsuitable for matrix multiplication."
                      << std::endl;
            abort();          
        }

        if (mat1.cols() != mat2.rows())
        {
            std::cerr << "Error in _ROWMAT_PROD_(RowMatrix<T1> const& mat1, " 
                      << "RowMatrix<T2> const& mat2)"
                      << "\n    RowMatrix 'mat1' with dimensions " 
                      << mat1.rows() << "x" << mat1.cols()
                      << " and RowMatrix 'mat2' with dimensions "
                      << mat2.rows() << "x" << mat2.cols()
                      << " are incompatible for matrix multiplication"
                      << std::endl; 
            abort();           
        }
        #endif
    }



   /***********************************************************************************************
    * RowMatrix<dtype> _ROWMAT_PROD_::eval() const
    *
    * Evaluate the expression template, creating a new RowMatrix. 
    * 
    * Template Parameters:
    *     T1 (implicit)    First RowMatrix data type
    *     T2 (implicit)    Second RowMatrix data type
    * 
    * Parameters:
    * 
    * Returns:
    *     RowMatrix    New RowMatrix object resulting from the evaluation of the multiplication
    ***********************************************************************************************/ 
    template<class T1, class T2>
    RowMatrix<typename _ROWMAT_PROD_<T1,T2>::dtype> _ROWMAT_PROD_<T1,T2>::eval() const
    {
        // Build RowMatrix
        RowMatrix<_ROWMAT_PROD_<T1,T2>::dtype> matrix 
            = RowMatrix<typename _ROWMAT_PROD_<T1,T2>::dtype>(this->rows(), this->cols());

        // Evaluate expression to newly allocated matrix
        this->eval(matrix);

        return matrix;
    }




   /**********************************************************************************************
    * void _ROWMAT_PROD_::eval(Matrix* target) const
    * 
    * Evaluates the expression template to target Matrix.
    * 
    * NOTE: Parallelization and cache blocking are not yet implemented. THIS MUST BE ADDED WHEN
    *       POSSIBLE!!!!
    * 
    * Template Parameters:
    *     T1 (implicit)        First RowMatrix data type
    *     T2 (implicit)        Second RowMatrix data type
    *     MAT_T (implicit)     (Matrix) Matrix type (default: RowMatrix)
    * 
    * Parameters:
    *     target    (Matrix*) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<class T1, class T2>
    template<typename MAT_T>
    void _ROWMAT_PROD_<T1, T2>::eval(Matrix<MAT_T>* target) const
    {
        if (target != nullptr) // Check that target matrix is allocated
        {
            #ifdef DIM_CHECK // Check that target matrix is correctly sized
            if (  ((this->rows_ != target->rows()) || (this->rows_ != 0))
                ||((this->cols_ != target->cols()) || (this->cols_ != 0)) )
            {
                std::cerr << "Error in 'void _ROWMAT_PROD_<T1,T2>::eval(Matrix<MAT_T>* target)"
                          << "\n    target Matrix with dimensions " << target->rows() << "x" 
                          << target->cols() << " incompatible with _ROWMAT_PROD_ with dimensions " 
                          << this->rows() << "x" << this->cols() 
                          << std::endl;
                abort();          
            }    
            #endif
            // Fall through to assignment
        }
        else // Allocate new RowMatrix as target
        {
            target = new RowMatrix<_ROWMAT_PROD_<T1,T2>::dtype>(this->rows(), this->cols());
        }

        // Evaluate expression to target matrix
        for (index_t i = 0; i < this->rows(); ++i)
        {
            for (index_t j = 0; j < this->cols(); ++j) 
            {
                target(i,j) = 0.0;
                for (index_t k = 0; k < this->mat1_.cols(); ++k)
                {
                    target(i,j) += this->mat1_(i,k)*this->mat2_(k,j);
                }
            }
        }
        return;
    }







   /***********************************************************************************************
    * _ROWMAT_PROD_ operator* (RowMatrix const& mat1, RowMatrix const& mat2)
    * 
    * Maps the _ROWMAT_PROD_ constructor to the (*) operator for mathematical syntax
    * 
    * Template Parameters:
    *     T1 (implicit)        First RowMatrix data type
    *     T2 (implicit)        Second RowMatrix data type
    * 
    * Parameters:
    *     mat1    (RowMatrix const&) Left hand side matrix
    *     mat2    (RowMatrix const&) Right hand side matrix
    * 
    * Returns:
    *     _ROWMAT_PROD_    Expression template for row major matrix multiplication
    ***********************************************************************************************/
    template<class T1, class T2>
    _ROWMAT_PROD_<T1, T2> operator* (RowMatrix<T1> const& mat1, RowMatrix<T2> const& mat2)
    {return _ROWMAT_PROD_<T1,T2>(mat1, mat2);} 




   /***********************************************************************************************
    * _ROWMAT_PROD_ operator* (RowMatrix const& mat, _ROWMAT_ATOMIC_EXPR_ const& expr)
    * 
    * Evaluates row-major atomic expressions and maps the resultant matrices to the expression
    * template constructor to an overload of the (*) operator for mathematical syntax
    * 
    * Template Parameters:
    *     T (implicit)         RowMatrix data type
    *     RMAE_T (implicit)    Row major atomic expression subtype
    * 
    * Parameters:
    *     mat     (RowMatrix const&) Left hand side matrix
    *     expr    (_ROWMAT_ATOMIC_EXPR_ const&) Expression for right hand side matrix
    * 
    * Returns:
    *     _ROWMAT_PROD_    Expression template for row major matrix multiplication
    ***********************************************************************************************/
    template<class T, typename RMAE_T>
    _ROWMAT_PROD_<T, typename _ROWMAT_ATOMIC_EXPR_<RMAE_T>::dtype> 
        operator* (RowMatrix<T> const& mat, _ROWMAT_ATOMIC_EXPR_<RMAE_T> const& expr)
    {return _ROWMAT_PROD_<T, typename TypeTraits<RMAE_T>::dtype>(mat, expr.eval());}


   
   /***********************************************************************************************
    * _ROWMAT_PROD_ operator* (_ROWMAT_ATOMIC_EXPR_ const& expr, RowMatrix const& mat)
    * 
    * Evaluates row-major atomic expressions and maps the resultant matrices to the expression
    * template constructor to an overload of the (*) operator for mathematical syntax
    * 
    * Template Parameters:
    *     RMAE_T (implicit)    Row major atomic expression subtype
    *     T (implicit)         RowMatrix data type
    *     
    * 
    * Parameters:
    *     expr    (_ROWMAT_ATOMIC_EXPR_ const&) Expression for left hand side matrix
    *     mat     (RowMatrix const&) Right hand side matrix
    *     
    * Returns:
    *     _ROWMAT_PROD_    Expression template for row major matrix multiplication
    ***********************************************************************************************/
    template<typename RMAE_T, class T>
    _ROWMAT_PROD_<typename _ROWMAT_ATOMIC_EXPR_<RMAE_T>::dtype, T> 
        operator* (_ROWMAT_ATOMIC_EXPR_<RMAE_T> const& expr, RowMatrix<T> const& mat)
    {return _ROWMAT_PROD_<typename TypeTraits<RMAE_T>::dtype,T>(expr.eval(), mat);}




   /***********************************************************************************************
    * _ROWMAT_PROD_ operator* (_ROWMAT_ATOMIC_EXPR_ const& expr, _ROWMAT_ATOMIC_EXPR_ const& expr2)
    * 
    * Evaluates row-major atomic expressions and maps the resultant matrices to the expression
    * template constructor to an overload of the (*) operator for mathematical syntax
    * 
    * Template Parameters:
    *     RMAE_T1 (implicit)    Row major atomic expression subtype of LHS expression
    *     RMAE_T2 (implicit)    Row major atomic expression subtype of RHS expression
    * 
    * Parameters:
    *     expr1    (_ROWMAT_ATOMIC_EXPR_ const&) Expression for left hand side matrix
    *     expr2    (_ROWMAT_ATOMIC_EXPR_ const&) Expression for right hand side matrix
    *     
    * Returns:
    *     _ROWMAT_PROD_    Expression template for row major matrix multiplication
    ***********************************************************************************************/
    template<typename RMAE_T1, typename RMAE_T2>
    _ROWMAT_PROD_<typename _ROWMAT_ATOMIC_EXPR_<RMAE_T1>::dtype, 
                  typename _ROWMAT_ATOMIC_EXPR_<RMAE_T2>::dtype> 
        operator* (_ROWMAT_ATOMIC_EXPR_<RMAE_T1> const& expr1, 
                   _ROWMAT_ATOMIC_EXPR_<RMAE_T1> const& expr2)
    {return _ROWMAT_PROD_<typename TypeTraits<RMAE_T1>::dtype,
                          typename TypeTraits<RMAE_T2>::dtype>(expr1.eval(), expr1.eval());}



}

#endif