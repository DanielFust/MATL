
namespace matl
{
   /**********************************************************************************************
    * _ROWMAT_PROD_(RowMatrix const& mat1, RowMatrix const& mat2)
    * 
    * Constructor for RowMatrix-RowMatrix multiplication expression templates. 
    * 
    * Template Parameters:
    *     RMAE1 (implicit)    RowMatrix or row major atomic expression type
    *     RMAE2 (implicit)    RowMatrix or row major atomic expression type
    * 
    * Parameters:
    *     mat1    (RowMatrix const&) RowMatrix or row major atomic expression
    *     mat2    (RowMatrix const&) RowMatrix or row major atomic expression
    * 
    * Returns:
    **********************************************************************************************/
    template<typename RMAE1, typename RMAE2>
    _ROWMAT_PROD_<RMAE1,RMAE2>::_ROWMAT_PROD_(RMAE1 const& mat1, 
                                              RMAE2 const& mat2) : mat1_(mat1), 
                                                                   mat2_(mat2) 
    {
        // Check for dimensional compatibility
        #ifdef MATL_DIM_CHECK 
        if ((mat1.rows() <= 0) || (mat1.cols() <= 0))
        {
            std::cerr << "Error in _ROWMAT_PROD_(RMAE1 const& mat1, " 
                      << "RMAE2 const& mat2)"
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
            std::cerr << "Error in _ROWMAT_PROD_(RMAE1 const& mat1, " 
                      << "RMAE2 const& mat2)"
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






   /**********************************************************************************************
    * void _ROWMAT_PROD_::eval(Matrix& target) const
    * 
    * Evaluates the expression template to target Matrix.
    * 
    * NOTE: Parallelization and cache blocking are not yet implemented. THIS MUST BE ADDED WHEN
    *       POSSIBLE!!!!
    * 
    * Template Parameters:
    *     RMAE1 (implicit)    RowMatrix or row major atomic expression type
    *     RMAE2 (implicit)    RowMatrix or row major atomic expression type
    *     MAT_T (implicit)     (Matrix) Matrix type (default: RowMatrix)
    * 
    * Parameters:
    *     target    (Matrix*) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    template<typename MAT_T>
    void _ROWMAT_PROD_<RMAE1,RMAE2>::eval(Matrix<MAT_T>& target) const
    {
        
        #ifdef MATL_DIM_CHECK // Check that target matrix is correctly sized
        if (  ((this->rows_ != target->rows()) || (this->rows_ != 0))
            ||((this->cols_ != target->cols()) || (this->cols_ != 0)) )
        {
            std::cerr << "Error in 'void _ROWMAT_PROD_<RMAE1,RMAE2>::eval(Matrix<MAT_T>* target)"
                        << "\n    target Matrix with dimensions " << target->rows() << "x" 
                        << target->cols() << " incompatible with _ROWMAT_PROD_ with dimensions " 
                        << this->rows() << "x" << this->cols() 
                        << std::endl;
            abort();          
        }    
        #endif    

        // Evaluate expression to target matrix
        for (index_t i = 0; i != this->rows()+1; ++i)
        {
            for (index_t j = 0; j != this->cols()+1; ++j) 
            {
                target(i,j) = 0.0;
                for (index_t k = 0; k != this->mat1_.cols()+1; ++k)
                {
                    target(i,j) += this->mat1_(i,k)*this->mat2_(k,j);
                }
            }
        }
        return;
    }





   /**********************************************************************************************
    * index_t rows() const
    * 
    * Gets number of rows in expression.
    * 
    * Template Parameters:
    *     RMAE1 (implicit)    RowMatrix or row major atomic expression type
    *     RMAE2 (implicit)    RowMatrix or row major atomic expression type
    * Parameters:
    * Returns:
    *     index_t             Number of rows in expression
    **********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    index_t _ROWMAT_PROD_<RMAE1,RMAE2>::rows() const
    {return this->mat1_.rows();}



   /**********************************************************************************************
    * index_t cols() const
    * 
    * Gets number of columns in expression.
    * 
    * Template Parameters:
    *     RMAE1 (implicit)    RowMatrix or row major atomic expression type
    *     RMAE2 (implicit)    RowMatrix or row major atomic expression type
    * Parameters:
    * Returns:
    *     index_t             Number of columns in expression
    **********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    index_t _ROWMAT_PROD_<RMAE1,RMAE2>::cols() const
    {return this->mat2_.cols();}





   /**********************************************************************************************
    * index_t size() const
    * 
    * Gets size of expression.
    * 
    * Template Parameters:
    *     RMAE1 (implicit)    RowMatrix or row major atomic expression type
    *     RMAE2 (implicit)    RowMatrix or row major atomic expression type
    * Parameters:
    * Returns:
    *     index_t             Number of elements in expression
    **********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    index_t _ROWMAT_PROD_<RMAE1,RMAE2>::size() const
    {return this->rows() * this->cols();}  







   /***********************************************************************************************
    * bool memory_conflict(Matrix const&) const
    * 
    * Check if there is memory conflict between expression and target matrix
    * 
    * Template Parameters:
    *     RMAE1 (implicit)    RowMatrix or row major atomic expression type
    *     RMAE2 (implicit)    RowMatrix or row major atomic expression type
    * 
    * Parameters:
    *     target             Target matrix
    * 
    * Returns:
    *     'true' if direct evaluation corrupts memory 'false' otherwise
    ***********************************************************************************************/  
    template<typename RMAE1, typename RMAE2>
    template<typename MAT_T>
    bool _ROWMAT_PROD_<RMAE1,RMAE2>::memory_conflict(Matrix<MAT_T> const& target) const
    {return (memory_conflict(this->mat1_, target) || memory_conflict(this->mat2_, target));} 












}
