
namespace matl
{
   /***********************************************************************************************
    * _GENERIC_MAT_SCALAR_SUB_(EXPR const&, T)
    *
    * Constructor for matrix-scalar subtraction expression template
    * 
    * Template Parameters:
    *     EXPR      (implicit) Atomic expression template or Matrix type
    *     T         (implicit) scalar data type
    * 
    * Parameters:
    *     mat        Atomic expression or Matrix
    *     val        Scalar value
    * 
    * Returns:
    ***********************************************************************************************/ 
     template<typename EXPR, class T>
    _GENERIC_MAT_SCALAR_SUB_<EXPR,T>::_GENERIC_MAT_SCALAR_SUB_(EXPR const& mat, T val) : 
        mat_(mat), val_(val) {}

   
   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    *
    * Retrieve evaluated (double) indexed value 
    * 
    * Template Parameters:
    *     EXPR      (implicit) Atomic expression template or Matrix type
    *     T         (implicit) scalar data type
    * 
    * Parameters:
    *     i       (index_t) Row index
    *     j       (index_t) Column index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix  
    ***********************************************************************************************/ 
    template<typename EXPR, class T>
    typename _GENERIC_MAT_SCALAR_SUB_<EXPR,T>::dtype 
        _GENERIC_MAT_SCALAR_SUB_<EXPR,T>::operator() (index_t i, index_t j) const
    {return this->mat_(i,j) - this->val_;} 




   /**********************************************************************************************
    * void eval(Matrix& target) const
    * 
    * Evaluates the expression to target Matrix
    * 
    * Template Parameters:
    *     EXPR      (implicit) Atomic expression template or Matrix type
    *     T         (implicit) scalar data type
    *     MAT_T     (implicit) (Matrix) Matrix type 
    * 
    * Parameters:
    *     target    (Matrix&) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename EXPR, class T>
    template<typename MAT_T>
    void _GENERIC_MAT_SCALAR_SUB_<EXPR,T>::eval(Matrix<MAT_T>& target) const
    {
        const index_t block = this->block_size();
        const index_t M = this->rows();
        const index_t N = this->cols();

        #ifdef MATL_DIM_CHECK // Check that target matrix is correctly sized
        if (  ((this->rows_ != target->rows()) || (this->rows_ != 0))
            ||((this->cols_ != target->cols()) || (this->cols_ != 0)) )
        {
            std::cerr << "Error in 'void _GENERIC_MAT_SCALAR_SUB_<EXPR,T>::eval(Matrix<MAT_T>& target)"
                    << "\n    target Matrix with dimensions " << target->rows() << "x" 
                    << target->cols() << " incompatible with _GENERIC_MAT_SUB_ with " 
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
    * Number of rows in matrix-scalar subtraction expression
    * 
    * Template Parameters:
    *     EXPR      (implicit) Atomic expression template or Matrix type
    *     T         (implicit) scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t    Number of rows in matrix subtraction expression
    **********************************************************************************************/ 
    template<typename EXPR, class T>
    index_t _GENERIC_MAT_SCALAR_SUB_<EXPR,T>::rows() const {return this->mat_.rows();}




   /**********************************************************************************************
    * index_t cols() const
    * 
    * Number of columns in matrix-scalar subtraction expression
    * 
    * Template Parameters:
    *     EXPR      (implicit) Atomic expression template or Matrix type
    *     T         (implicit) scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t    Number of columns in matrix subtraction expression
    **********************************************************************************************/ 
    template<typename EXPR, class T>
    index_t _GENERIC_MAT_SCALAR_SUB_<EXPR,T>::cols() const {return this->mat_.cols();} 





   /**********************************************************************************************
    * index_t size() const
    * 
    * Size in matrix-scalar subtraction expression
    * 
    * Template Parameters:
    *     EXPR      (implicit) Atomic expression template or Matrix type
    *     T         (implicit) scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t    Number of rows in matrix subtraction expression
    **********************************************************************************************/ 
    template<typename EXPR, class T>
    index_t _GENERIC_MAT_SCALAR_SUB_<EXPR,T>::size() const {return this->mat_.size();} 



    
   /***********************************************************************************************
    * index_t block_size()
    *
    * Computes a cache block size based on current thread's chunk size
    ***********************************************************************************************/ 
    template<typename EXPR, class T>
    index_t _GENERIC_MAT_SCALAR_SUB_<EXPR,T>::block_size()
    {
        index_t type_size = sizeof(_GENERIC_MAT_SCALAR_SUB_<EXPR,T>::dtype);
        return chunk_size()/type_size;
    }  

}
