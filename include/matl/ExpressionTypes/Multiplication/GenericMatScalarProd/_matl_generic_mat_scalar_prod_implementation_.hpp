
namespace matl
{
 
   /***********************************************************************************************
    * _GENERIC_MAT_SCALAR_PROD_(AE_T const&, T) 
    *
    * Constructor for generic matrix scalar product
    * 
    * Template Parameters:
    *     AE     (implicit) Atomic expression or matrix type
    *     T      (implicit) Scalar data type
    * 
    * Parameters:
    *     mat     Matrix or atomic expression
    *     val     Scalar multiplier
    * 
    * Returns: 
    ***********************************************************************************************/ 
    template<typename AE, class T>
    _GENERIC_MAT_SCALAR_PROD_<AE,T>::_GENERIC_MAT_SCALAR_PROD_(AE const& mat, T val) : 
    mat_(mat), val_(val), rows_(mat.rows()), cols_(mat.cols()) {}






   /***********************************************************************************************
    * dtype _GENERIC_MAT_SCALAR_PROD_::operator() (index_t i, index_t j) const
    *
    * Retrieve evaluated (double) indexed value 
    * 
    * Template Parameters:
    *     AE     (implicit) (_ATOMIC_EXPR_ or Matrix) Generic atomic matrix expression
    *     T      (implicit) Scalar data type
    * 
    * Parameters:
    *     i       (index_t) Row index
    *     j       (index_t) Column index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix  
    ***********************************************************************************************/ 
    template<typename AE, class T>
    typename _GENERIC_MAT_SCALAR_PROD_<AE,T>::dtype 
        _GENERIC_MAT_SCALAR_PROD_<AE,T>::operator() (index_t i, index_t j) const
    {return this->mat_(i,j) * this->val_;} 




   /**********************************************************************************************
    * void _GENERIC_MAT_SUM_::eval(Matrix&) const
    * 
    * Evaluates the expression to target Matrix
    * 
    * Template Parameters:
    *     AE     (implicit) (_ATOMIC_EXPR_ or Matrix) Generic atomic matrix expression
    *     T      (implicit) Scalar data type
    *     MAT_T     (implicit) (Matrix) Target matrix type 
    * 
    * Parameters:
    *     target    (Matrix&) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename AE, class T>
    template<typename MAT_T>
    void _GENERIC_MAT_SCALAR_PROD_<AE,T>::eval(Matrix<MAT_T>& target) const
    {
        const index_t block = this->block_size();
        const index_t M = this->rows();
        const index_t N = this->cols();

        #ifdef MATL_DIM_CHECK // Check that target matrix is correctly sized
        if  ((this->rows() != target.rows()) || (this->cols() != target.cols()))
        {
            std::cerr << "Error in 'void _GENERIC_MAT_SCALAR_PROD_<AE,T>::eval(Matrix<MAT_T>& target)"
                    << "\n    target Matrix with dimensions " << target.rows() << "x" 
                    << target.cols() << " incompatible with _GENERIC_MAT_SCALAR_PROD_ with " 
                    << "dimensions " << this->rows() << "x" << this->cols() 
                    << std::endl;
            abort();          
        }    
        #endif

        // Cache blocking to mitigate miss-rate for intermingled row and column major ordering
        for (index_t ii=1; ii<(M+1); ii+=block)
            for (index_t jj=1; jj<(N+1); jj+=block)
            {
                const index_t istop = branchless_min(ii+block, M+1);
                const index_t jstop = branchless_min(jj+block, N+1);
                for (index_t i=ii; i!=istop; ++i)
                    for (index_t j=jj; j!=jstop; ++j)
                        target(i,j) = (*this)(i,j);
            }
  
        return;
    }
    


   /**********************************************************************************************
    * index_t rows() const
    * 
    * Gets number of rows in expression
    * 
    * Template Parameters:
    *     AE     (implicit) Matrix or atomic expression type
    *     T      (implicit) Scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t      Number of rows in expression
    **********************************************************************************************/ 
    template<typename AE, class T>
    index_t _GENERIC_MAT_SCALAR_PROD_<AE,T>::rows() const
    {return mat_.rows();}


   /**********************************************************************************************
    * index_t cols() const
    * 
    * Gets number of columns in expression
    * 
    * Template Parameters:
    *     AE     (implicit) Matrix or atomic expression type
    *     T      (implicit) Scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t      Number of rows in expression
    **********************************************************************************************/ 
    template<typename AE, class T>
    index_t _GENERIC_MAT_SCALAR_PROD_<AE,T>::cols() const
    {return mat_.cols();}     




   /**********************************************************************************************
    * index_t size() const
    * 
    * Gets size of resultant matrix
    * 
    * Template Parameters:
    *     AE     (implicit) Matrix or atomic expression type
    *     T      (implicit) Scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t      size of resultant matrix
    **********************************************************************************************/ 
    template<typename AE, class T>
    index_t _GENERIC_MAT_SCALAR_PROD_<AE,T>::size() const
    {return mat_.size();}      
       





    
   /***********************************************************************************************
    * index_t block_size()
    *
    * Computes a cache block size based on current thread's chunk size
    ***********************************************************************************************/ 
    template<typename AE, class T>
    index_t _GENERIC_MAT_SCALAR_PROD_<AE,T>::block_size() const
    {
        index_t type_size = sizeof(_GENERIC_MAT_SCALAR_PROD_<AE,T>::dtype);
        return chunk_size()/type_size;
    }  


}
