//#ifndef _MATL_ROWMAT_PROD__HPP
//#define _MATL_ROWMAT_PROD__HPP

//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"

namespace matl
{
  
   /**********************************************************************************************
    * _GENERIC_MAT_PROD_(RowMatrix const& mat1, RowMatrix const& mat2)
    * 
    * Constructor for RowMatrix-RowMatrix multiplication expression templates. 
    * 
    * Template Parameters:
    *     MAT1 (implicit)    Matrix or atomic matrix expression type
    *     MAT2 (implicit)    Matrix or atomic matrix expression type
    * 
    * Parameters:
    *     mat1               Matrix or atomic matrix expression
    *     mat2               Matrix or atomic matrix expression
    * 
    * Returns:
    **********************************************************************************************/
    template<typename MAT1, typename MAT2>
    _GENERIC_MAT_PROD_<MAT1,MAT2>::_GENERIC_MAT_PROD_(MAT1 const& mat1, 
                                                      MAT2 const& mat2) : mat1_(mat1), 
                                                                          mat2_(mat2),
                                                                          rows_(mat1.rows()),
                                                                          cols_(mat2.cols())                                                                       
    {
        // Check for dimensional compatibility
        #ifdef MATL_DIM_CHECK 
        if ((mat1.rows() <= 0) || (mat1.cols() <= 0))
        {
            std::cerr << "Error in _GENERIC_MAT_PROD_(MAT1 const& mat1, " 
                      << "MAT2 const& mat2)"
                      << "\n    Matrix 'mat1' with dimensions " 
                      << mat1.rows() << "x" << mat1.cols()
                      << std::endl;
            abort();          
        }

        if ((mat2.rows() <= 0) || (mat2.cols() <= 0))
        {
            std::cerr << "Error in _GENERIC_MAT_PROD_(MAT1 const& mat1, " 
                      << "MAT2 const& mat2)"
                      << "\n    Matrix 'mat2' with dimensions " 
                      << mat2.rows() << "x" << mat2.cols()
                      << " unsuitable for matrix multiplication."
                      << std::endl;
            abort();          
        }

        if (mat1.cols() != mat2.rows())
        {
            std::cerr << "Error in _GENERIC_MAT_PROD_(MAT1 const& mat1, " 
                      << "MAT2 const& mat2)"
                      << "\n    Matrix 'mat1' with dimensions " 
                      << mat1.rows() << "x" << mat1.cols()
                      << " and Matrix 'mat2' with dimensions "
                      << mat2.rows() << "x" << mat2.cols()
                      << " are incompatible for matrix multiplication"
                      << std::endl; 
            abort();           
        }
        #endif
    }





   /**********************************************************************************************
    * void eval(Matrix& target) const
    * 
    * Evaluates the expression template to target Matrix.
    * 
    * NOTE: Parallelization and cache blocking are not yet implemented. THIS MUST BE ADDED WHEN
    *       POSSIBLE!!!!
    * 
    * Template Parameters:
    *     MAT1 (implicit)    Matrix or atomic matrix expression type
    *     MAT2 (implicit)    Matrix or atomic matrix expression type
    *     MAT  (implicit)    (Matrix) Target matrix type 
    * 
    * Parameters:
    *     target             (Matrix&) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename MAT1, typename MAT2>
    template<typename MAT>
    void _GENERIC_MAT_PROD_<MAT1,MAT2>::eval(Matrix<MAT>& target) const
    {
        
        #ifdef MATL_DIM_CHECK // Check that target matrix is correctly sized
        if (((this->rows() != target.rows() )||(this->cols() != target.cols())) 
            && (target.size()!=0))
        {
            std::cerr << "Error in 'void _GENERIC_MAT_PROD_<MAT1,MAT2>::eval(Matrix<MAT>& target)"
                        << "\n    target Matrix with dimensions " << target.rows() << "x" 
                        << target.cols() << " incompatible with _GENERIC_MAT_PROD_ with dimensions " 
                        << this->rows() << "x" << this->cols() 
                        << std::endl;
            abort();          
        }    
        #endif    

        if (target.size()==0) target.resize(this->rows(), this->cols());

        // Evaluate expression to target matrix
        for (index_t i = 1; i != this->rows()+1; ++i)
        {
            for (index_t j = 1; j != this->cols()+1; ++j) 
            {
                target(i,j) = 0.0;
                for (index_t k = 1; k != this->mat1_.cols()+1; ++k)
                {
                    target(i,j) += this->mat1_(i,k)*this->mat2_(k,j);
                }
            }
        }
        return;
    }




   /**********************************************************************************************
    * Matrix eval() const
    * 
    * Evaluates the expression template to a new Matrix.
    * 
    * 
    * Template Parameters:
    *     MAT1 (implicit)    Matrix or atomic matrix expression type
    *     MAT2 (implicit)    Matrix or atomic matrix expression type
    * 
    * Parameters:
    * Returns:
    *     Matrix with first matrix type
    **********************************************************************************************/ 
    template<typename MAT1, typename MAT2> 
    MAT1 _GENERIC_MAT_PROD_<MAT1,MAT2>::eval() const
    {
        typename TypeTraits<MAT1>::matrix_type mat;
        this->eval(mat);
        return mat;
    }





   /**********************************************************************************************
    * index_t rows() const
    * 
    * Gets number of rows in expression.
    * 
    * Template Parameters:
    *     MAT1 (implicit)     Matrix or atomic matrix expression type
    *     MAT2 (implicit)     Matrix or atomic matrix expression type
    * Parameters:
    * Returns:
    *     index_t             Number of rows in expression
    **********************************************************************************************/ 
    template<typename MAT1, typename MAT2>
    index_t _GENERIC_MAT_PROD_<MAT1,MAT2>::rows() const
    {return this->rows_;}



   /**********************************************************************************************
    * index_t cols() const
    * 
    * Gets number of columns in expression.
    * 
    * Template Parameters:
    *     MAT1 (implicit)    Matrix or atomic matrix expression type
    *     MAT2 (implicit)    Matrix or atomic matrix expression type
    * Parameters:
    * Returns:
    *     index_t             Number of columns in expression
    **********************************************************************************************/ 
    template<typename MAT1, typename MAT2>
    index_t _GENERIC_MAT_PROD_<MAT1,MAT2>::cols() const
    {return this->cols_;}





   /**********************************************************************************************
    * index_t cols() const
    * 
    * Gets size of expression.
    * 
    * Template Parameters:
    *     MAT1 (implicit)    Matrix or atomic matrix expression type
    *     MAT2 (implicit)    Matrix or atomic matrix expression type
    * Parameters:
    * Returns:
    *     index_t             Number of elements in expression
    **********************************************************************************************/ 
    template<typename MAT1, typename MAT2>
    index_t _GENERIC_MAT_PROD_<MAT1,MAT2>::size() const
    {return this->rows() * this->cols();}  




   /***********************************************************************************************
    * bool memory_conflict(Matrix const&) const
    * 
    * Check if there is memory conflict between expression and target matrix
    * 
    * Template Parameters:
    *     MAT1 (implicit)    Matrix or atomic matrix expression type
    *     MAT2 (implicit)    Matrix or atomic matrix expression type
    * 
    * Parameters:
    *     target             Target matrix
    * 
    * Returns:
    *     'true' if direct evaluation corrupts memory 'false' otherwise
    ***********************************************************************************************/  
    template<typename MAT1, typename MAT2>
    template<typename MAT_T>
    bool _GENERIC_MAT_PROD_<MAT1,MAT2>::memory_conflict(Matrix<MAT_T> const& target) const
    {return (address_overlap(this->mat1_, target) || address_overlap(this->mat2_, target));} 





   /***********************************************************************************************
    * index_t block_size()
    *
    * Computes a cache block size based on current thread's chunk size
    * 
    * Template Parameters:
    *     MAT1 (implicit)    Matrix or atomic matrix expression type
    *     MAT2 (implicit)    Matrix or atomic matrix expression type
    * Parameters:
    * Returns:
    *     index_t             Number of elements in a block
    ***********************************************************************************************/ 
    template<typename MAT1, typename MAT2>
    index_t _GENERIC_MAT_PROD_<MAT1,MAT2>::block_size() const
    {
        index_t type_size = sizeof(_GENERIC_MAT_PROD_<MAT1,MAT2>::dtype);
        return chunk_size()/type_size;
    }


}

//#endif