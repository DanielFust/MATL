//#ifndef _MATL_GENERIC_MAT_SUB__HPP
//#define _MATL_GENERIC_MAT_SUB__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 **************************************************************************************************/
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"

namespace matl
{
   /***********************************************************************************************
    * _GENERIC_MAT_SUB_(AE_T1 const&, AE_T2 const&)
    *
    * Constructor for matrix subtraction expression template
    * 
    * Template Parameters:
    *     AE_T1     (implicit) Atomic expression template or Matrix type
    *     AE_T2     (implicit) Atomic expression template or Matrix type
    * 
    * Parameters:
    *     mat1       Atomic expression or Matrix
    *     mat2       Atomic expression or Matrix
    * 
    * Returns:
    ***********************************************************************************************/ 
     template<typename AE_T1, typename AE_T2>
    _GENERIC_MAT_SUB_<AE_T1,AE_T2>::_GENERIC_MAT_SUB_(AE_T1 const& mat1, AE_T2 const& mat2) : 
        mat1_(mat1), mat2_(mat2) {}

   
   /***********************************************************************************************
    * dtype _GENERIC_MAT_SUB_::operator() (index_t i, index_t j) const
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
    typename _GENERIC_MAT_SUB_<AE_T1,AE_T2>::dtype 
        _GENERIC_MAT_SUB_<AE_T1,AE_T2>::operator() (index_t i, index_t j) const
    {return this->mat1_(i,j) - this->mat2_(i,j);} 




   /**********************************************************************************************
    * void _GENERIC_MAT_SUB_::eval(Matrix& target) const
    * 
    * Evaluates the expression to target Matrix
    * 
    * Template Parameters:
    *     AE_T1     (implicit) (_ATOMIC_EXPR_) Generic atomic matrix expression
    *     AE_T2     (implicit) (_ATOMIC_EXPR_) Generic atomic matrix expression
    *     MAT_T     (implicit) (Matrix) Matrix type 
    * 
    * Parameters:
    *     target    (Matrix&) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename AE1, typename AE2>
    template<typename MAT_T>
    void _GENERIC_MAT_SUB_<AE1,AE2>::eval(Matrix<MAT_T>& target) const
    {
        const index_t block = this->block_size();
        const index_t M = this->rows();
        const index_t N = this->cols();

        #ifdef MATL_DIM_CHECK // Check that target matrix is correctly sized
        if (  ((this->rows() != target.rows()) || (this->cols() != target.cols())) && target.size() != 0)
            
        {
            std::cerr << "Error in 'void _GENERIC_MAT_SUB_<E1,E2>::eval(Matrix<MAT_T>* target)"
                    << "\n    target Matrix with dimensions " << target.rows() << "x" 
                    << target.cols() << " incompatible with _GENERIC_MAT_SUB_ with " 
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
    * Number of rows in matrix subtraction expression
    * 
    * Template Parameters:
    *     AE_T1     (implicit) (_ATOMIC_EXPR_) Generic atomic matrix expression
    *     AE_T2     (implicit) (_ATOMIC_EXPR_) Generic atomic matrix expression
    * 
    * Parameters:
    * Returns:
    *     index_t    Number of rows in matrix subtraction expression
    **********************************************************************************************/ 
    template<typename AE1, typename AE2>
    index_t _GENERIC_MAT_SUB_<AE1,AE2>::rows() const {return mat1_.rows();}


   /**********************************************************************************************
    * index_t cols() const
    * 
    * Number of columns in matrix subtraction expression
    * 
    * Template Parameters:
    *     AE_T1     (implicit) (_ATOMIC_EXPR_) Generic atomic matrix expression
    *     AE_T2     (implicit) (_ATOMIC_EXPR_) Generic atomic matrix expression
    * 
    * Parameters:
    * Returns:
    *     index_t    Number of columns in matrix subtraction expression
    **********************************************************************************************/ 
    template<typename AE1, typename AE2>
    index_t _GENERIC_MAT_SUB_<AE1,AE2>::cols() const {return mat1_.cols();}


   /**********************************************************************************************
    * index_t size() const
    * 
    * Number of elements in matrix subtraction expression
    * 
    * Template Parameters:
    *     AE_T1     (implicit) (_ATOMIC_EXPR_) Generic atomic matrix expression
    *     AE_T2     (implicit) (_ATOMIC_EXPR_) Generic atomic matrix expression
    * 
    * Parameters:
    * Returns:
    *     index_t    Number of elements in matrix subtraction expression
    **********************************************************************************************/ 
    template<typename AE1, typename AE2>
    index_t _GENERIC_MAT_SUB_<AE1,AE2>::size() const {return mat1_.size();}  

    


    
   /***********************************************************************************************
    * index_t block_size()
    *
    * Computes a cache block size based on current thread's chunk size
    ***********************************************************************************************/ 
    template<typename AE_T1, typename AE_T2>
    index_t _GENERIC_MAT_SUB_<AE_T1,AE_T2>::block_size() const
    {
        index_t type_size = sizeof(_GENERIC_MAT_SUB_<AE_T1,AE_T2>::dtype);
        return chunk_size()/type_size;
    }  


}
//#endif