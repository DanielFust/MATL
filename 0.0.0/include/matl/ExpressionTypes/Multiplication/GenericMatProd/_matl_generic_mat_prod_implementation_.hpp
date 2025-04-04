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
    *     MAT1 (implicit)    Matrix or matrix expression type
    *     MAT2 (implicit)    Matrix or matrix expression type
    * 
    * Parameters:
    *     mat1               Matrix or matrix expression
    *     mat2               Matrix or matrix expression
    * 
    * Returns:
    **********************************************************************************************/
    template<typename MAT1, typename MAT2>
    _GENERIC_MAT_PROD_<MAT1,MAT2>::_GENERIC_MAT_PROD_(MAT1 const& mat1, 
                                                      MAT2 const& mat2)                                                                      
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


        this->mat_.resize(mat1.rows(),mat2.cols());

        // Evaluate expression to member matrix
        for (index_t i = 1; i != mat1.rows()+1; ++i)
        {
            for (index_t j = 1; j != mat2.cols()+1; ++j) 
            {
                this->mat_(i,j) = 0.0;
                for (index_t k = 1; k != mat1.cols()+1; ++k)
                {
                    this->mat_(i,j) += mat1(i,k)*mat2(k,j);
                }
            }
        }
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
        /*
        for (index_t i = 1; i != this->rows()+1; ++i)
        {
            for (index_t j = 1; j != this->cols()+1; ++j) 
            {
                target(i,j) = 0.0;
                for (index_t k = 1; k != mat1.cols()+1; ++k)
                {
                    target(i,j) += this->mat1_(i,k)*this->mat2_(k,j);
                }
            }
        }
        */
        
        for (index_t i = 1; i != this->rows()+1; ++i)
            for (index_t j = 1; j != this->cols()+1; ++j) 
                target(i,j) = (*this)(i,j);
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




   /***********************************************************************************************
    * dtype _GENERIC_MAT_PROD_::operator() (index_t i, index_t j) const
    *
    * Retrieve evaluated (double) indexed value 
    * 
    * Template Parameters:
    *     AE_T1     (implicit) (_ATOMIC_EXPR_) Generic atomic matrix expression
    *     AE_T2     (implicit) (_ATOMIC_EXPR_) Generic atomic matrix expression
    * 
    * Parameters:
    *     i       (index_t) Row index
    *     j       (index_t) Column index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix  
    ***********************************************************************************************/ 
    template<typename AE_T1, typename AE_T2>
    typename _GENERIC_MAT_PROD_<AE_T1,AE_T2>::dtype 
        _GENERIC_MAT_PROD_<AE_T1,AE_T2>::operator() (index_t i, index_t j) const
    {return this->mat_(i,j);} 









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
    {return this->mat_.rows();}



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
    {return this->mat_.cols();}





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
    {return false;}
    //{return (address_overlap(this->mat_, target));} 





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