#ifndef _MATL_GENERIC_MAT_SCALAR_SUM_IMPLEMENTATION__HPP
#define _MATL_GENERIC_MAT_SCALAR_SUM_IMPLEMENTATION__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method delarations and implementations of Row Matrix Scalar 
 * multiplications
 **************************************************************************************************/
//#include "_matl_forward_declr_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"

namespace matl
{
    
   /***********************************************************************************************
    * _GENERIC_MAT_SCALAR_SUM_(AE_T const& mat, T val)
    * 
    * Constructor
    * 
    * Template Parameters:
    *     AE_T     (implicit) (_ATOMIC_EXPR_) Matrix or atomic matrix expression
    *     T        (implicit) Scalar data type
    * 
    * Parameters:
    *     mat       Matrix or atomic matrix expression
    *     val       Scalar     
    * 
    * Returns:
    ***********************************************************************************************/ 
    template<typename AE_T, class T>
    _GENERIC_MAT_SCALAR_SUM_<AE_T,T>::_GENERIC_MAT_SCALAR_SUM_(AE_T const& mat, T val) : 
        mat_(mat), val_(val) {}

     


   /***********************************************************************************************
    * dtype _GENERIC_MAT_SCALAR_SUM_::operator() (index_t i) const
    * 
    * Retrieve evaluated (single) indexed value
    * 
    * Template Parameters:
    *     AE_T     (implicit) (_ATOMIC_EXPR_) Matrix or atomic matrix expression
    *     T        (implicit) Scalar data type
    * 
    * Parameters:
    *     i       (index_t) Element index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix 
    ***********************************************************************************************/ 
    template<typename AE_T, class T>
    typename _GENERIC_MAT_SCALAR_SUM_<AE_T,T>::dtype 
        _GENERIC_MAT_SCALAR_SUM_<AE_T,T>::operator() (index_t i) const
    {return this->mat_(i) + this->val_;}




   /***********************************************************************************************
    * dtype _GENERIC_MAT_SCALAR_SUM_::operator() (index_t i, index_t j) const
    *
    * Retrieve evaluated (double) indexed value 
    * 
    * Template Parameters:
    *     AE_T     (implicit) (_ATOMIC_EXPR_) Matrix or atomic matrix expression
    *     T        (implicit) Scalar data type
    * 
    * Parameters:
    *     i       (index_t) Row index
    *     j       (index_t) Column index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix  
    ***********************************************************************************************/ 
    template<typename AE_T, class T>
    typename _GENERIC_MAT_SCALAR_SUM_<AE_T,T>::dtype 
        _GENERIC_MAT_SCALAR_SUM_<AE_T,T>::operator() (index_t i, index_t j) const
    {return this->mat_(i,j) + this->val_;} 





   /**********************************************************************************************
    * void _GENERIC_MAT_SCALAR_SUM_::eval(Matrix& target) const
    * 
    * Evaluates the expression to target Matrix of any kind.
    * 
    * Template Parameters:
    *     AE_T     (implicit) (_ATOMIC_EXPR_) Matrix or atomic matrix expression
    *     T        (implicit) Scalar data type
    *     MAT_T     (implicit) (Matrix) Matrix type 
    * 
    * Parameters:
    *     target    (Matrix*) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename AE_T, class T>
    template<typename MAT_T>
    void _GENERIC_MAT_SCALAR_SUM_<AE_T,T>::eval(Matrix<MAT_T>& target) const
    {
        const index_t block = this->block_size();
        const index_t M = this->rows();
        const index_t N = this->cols();
        
        #ifdef MATL_DIM_CHECK // Check that target matrix is correctly sized
        if (  (this->rows() != target.rows()) || (this->cols() != target.cols()) 
            &&(this->size()!=0) )
        {
            std::cerr << "Error in 'void _GENERIC_MAT_SCALAR_PROD_<AE_T,T>::eval(Matrix<MAT_T>* target)"
                      << "\n    target Matrix with dimensions " << target.rows() << "x" 
                      << target.cols() << " incompatible with _GENERIC_MAT_SCALAR_SUM_ with dimensions " 
                      << this->rows() << "x" << this->cols() 
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
    *     AE_T      (implicit) Matrix or matrix expression type
    *     T         (implicit) Scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t   Number of rows in expression
    **********************************************************************************************/ 
    template<typename AE_T, class T>
    index_t _GENERIC_MAT_SCALAR_SUM_<AE_T,T>::rows() const 
    {return this->mat_.rows();}
        

   /**********************************************************************************************
    * index_t cols() const
    * 
    * Gets number of columns in expression
    * 
    * Template Parameters:
    *     AE_T      (implicit) Matrix or matrix expression type
    *     T         (implicit) Scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t   Number of columns in expression
    **********************************************************************************************/ 
    template<typename AE_T, class T>
    index_t _GENERIC_MAT_SCALAR_SUM_<AE_T,T>::cols() const 
    {return this->mat_.cols();}



   /**********************************************************************************************
    * index_t size() const
    * 
    * Gets size of expression
    * 
    * Template Parameters:
    *     AE_T      (implicit) Matrix or matrix expression type
    *     T         (implicit) Scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t   Number of elements in expression
    **********************************************************************************************/ 
    template<typename AE_T, class T>
    index_t _GENERIC_MAT_SCALAR_SUM_<AE_T,T>::size() const 
    {return this->mat_.size();}






   /***********************************************************************************************
    * index_t block_size()
    *
    * Computes a cache block size based on current thread's chunk size
    ***********************************************************************************************/ 
    template<typename AE_T, class T>
    index_t _GENERIC_MAT_SCALAR_SUM_<AE_T,T>::block_size() const
    {
        index_t type_size = sizeof(_GENERIC_MAT_SCALAR_SUM_<AE_T,T>::dtype);
        return chunk_size()/type_size;
    }  
    

}


#endif