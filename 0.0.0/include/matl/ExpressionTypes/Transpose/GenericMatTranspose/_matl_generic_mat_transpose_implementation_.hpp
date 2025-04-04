//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"

namespace matl
{
  
   /**********************************************************************************************
    * _GENERIC_MAT_TRANSPOSE_(MAT const& mat)
    * 
    * Constructor for generic matrix transposition. 
    * 
    * Template Parameters:
    *     MAT (implicit)    Matrix or matrix expression type
    * 
    * Parameters:
    *     mat               Matrix or matrix expression
    * 
    * Returns:
    **********************************************************************************************/
    template<typename MAT>
    _GENERIC_MAT_TRANSPOSE_<MAT>::_GENERIC_MAT_TRANSPOSE_(MAT const& mat) : mat_(mat) {}                                                                      
    




   /**********************************************************************************************
    * void eval(Matrix& target) const
    * 
    * Evaluates the expression template to target Matrix.
    * 
    * NOTE: Parallelization and cache blocking are not yet implemented. THIS MUST BE ADDED WHEN
    *       POSSIBLE!!!!
    * 
    * Template Parameters:
    *     MAT    (implicit)    Matrix or atomic matrix expression type
    *     MAT_T  (implicit)    (Matrix) Target matrix type 
    * 
    * Parameters:
    *     target             (Matrix&) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename MAT>
    template<typename MAT_T>
    void _GENERIC_MAT_TRANSPOSE_<MAT>::eval(Matrix<MAT_T>& target) const
    {   
        #ifdef MATL_DIM_CHECK // Check that target matrix is correctly sized
        if (((this->rows() != target.rows() )||(this->cols() != target.cols())) 
            && (target.size()!=0))
        {
            std::cerr << "Error in 'void _GENERIC_MAT_TRANSPOSE_<MAT>::eval(Matrix<MAT_T>& target)"
                        << "\n    target Matrix with dimensions " << target.rows() << "x" 
                        << target.cols() << " incompatible with _GENERIC_MAT_TRANSPOSE_ with dimensions " 
                        << this->rows() << "x" << this->cols() 
                        << std::endl;
            abort();          
        }    
        #endif    

        // Self-target
        if (this->memory_conflict(target))
        {
            for (index_t i=1; i != this->rows()+1; ++i)
                for (index_t j=i+1; j <= this->cols(); ++j)
                {
                    typename TypeTraits<MAT_T>::dtype val = target(i,j);
                    target(i,j) =  (*this)(i,j);
                    target(j,i) = val;
                }    
        }
        else
        {
            if (target.size()==0) target.resize(this->rows(), this->cols());
        
            for (index_t i = 1; i != this->rows()+1; ++i)
                for (index_t j = 1; j != this->cols()+1; ++j) 
                    target(i,j) = (*this)(i,j);
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
    *     MAT (implicit)    Matrix or atomic matrix expression type
    * 
    * Parameters:
    * Returns:
    *     Matrix with deduced matrix type
    **********************************************************************************************/ 
    template<typename MAT> 
    MAT _GENERIC_MAT_TRANSPOSE_<MAT>::eval() const
    {
        typename TypeTraits<MAT>::matrix_type mat;
        this->eval(mat);
        return mat;
    }




   /***********************************************************************************************
    * dtype _GENERIC_MAT_TRANSPOSE_::operator() (index_t i, index_t j) const
    *
    * Retrieve evaluated (double) indexed value 
    * 
    * Template Parameters:
    *     MAT (implicit)    Matrix or atomic matrix expression type
    * 
    * Parameters:
    *     i       (index_t) Row index
    *     j       (index_t) Column index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix  
    ***********************************************************************************************/ 
    template<typename MAT>
    typename _GENERIC_MAT_TRANSPOSE_<MAT>::dtype 
        _GENERIC_MAT_TRANSPOSE_<MAT>::operator() (index_t i, index_t j) const
    {return this->mat_(j,i);} 









   /**********************************************************************************************
    * index_t rows() const
    * 
    * Gets number of rows in expression.
    * 
    * Template Parameters:
    *     MAT (implicit)    Matrix or atomic matrix expression type
    * Parameters:
    * Returns:
    *     index_t             Number of rows in expression
    **********************************************************************************************/ 
    template<typename MAT>
    index_t _GENERIC_MAT_TRANSPOSE_<MAT>::rows() const
    {return this->mat_.cols();}



   /**********************************************************************************************
    * index_t cols() const
    * 
    * Gets number of columns in expression.
    * 
    * Template Parameters:
    *     MAT (implicit)    Matrix or atomic matrix expression type
    * Parameters:
    * Returns:
    *     index_t             Number of columns in expression
    **********************************************************************************************/ 
    template<typename MAT>
    index_t _GENERIC_MAT_TRANSPOSE_<MAT>::cols() const
    {return this->mat_.rows();}





   /**********************************************************************************************
    * index_t cols() const
    * 
    * Gets size of expression.
    * 
    * Template Parameters:
    *     MAT (implicit)    Matrix or atomic matrix expression type
    * Parameters:
    * Returns:
    *     index_t             Number of elements in expression
    **********************************************************************************************/ 
    template<typename MAT>
    index_t _GENERIC_MAT_TRANSPOSE_<MAT>::size() const
    {return this->rows() * this->cols();}  




   /***********************************************************************************************
    * bool memory_conflict(Matrix const&) const
    * 
    * Check if there is memory conflict between expression and target matrix
    * 
    * Template Parameters:
    *     MAT (implicit)    Matrix or atomic matrix expression type
    * 
    * Parameters:
    *     target             Target matrix
    * 
    * Returns:
    *     'true' if direct evaluation corrupts memory 'false' otherwise
    ***********************************************************************************************/  
    template<typename MAT>
    template<typename MAT_T>
    bool _GENERIC_MAT_TRANSPOSE_<MAT>::memory_conflict(Matrix<MAT_T> const& target) const
    {return (address_overlap(this->mat_, target));}





   /***********************************************************************************************
    * index_t block_size()
    *
    * Computes a cache block size based on current thread's chunk size
    * 
    * Template Parameters:
    *     MAT (implicit)    Matrix or atomic matrix expression type
    * Parameters:
    * Returns:
    *     index_t             Number of elements in a block
    ***********************************************************************************************/ 
    template<typename MAT>
    index_t _GENERIC_MAT_TRANSPOSE_<MAT>::block_size() const
    {
        index_t type_size = sizeof(_GENERIC_MAT_TRANSPOSE_<MAT>::dtype);
        return chunk_size()/type_size;
    }


}

