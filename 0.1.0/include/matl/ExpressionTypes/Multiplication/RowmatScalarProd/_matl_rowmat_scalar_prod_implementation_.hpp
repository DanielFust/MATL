//#ifndef _MATL_ROWMAT_SCALAR_PROD_IMPLEMENTATION__HPP
//#define _MATL_ROWMAT_SCALAR_PROD_IMPLEMENTATION__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method delarations and implementations of Row Matrix Scalar 
 * multiplications
 **************************************************************************************************/
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"

namespace matl
{
   
   /***********************************************************************************************
    * _ROWMAT_SCALAR_PROD_(RMAE const&, T)
    * 
    * Constructor for RowMatrix - scalar multiplication expression template
    * 
    * Template Parameters:
    *     RMAE     (implicit) RowMatrix or row major atomic expression
    *     T        (implicit) Scalar data type
    * 
    * Parameters:
    *     mat       RowMatrix or _ROWMAT_ATOMIC_EXPR_ 
    * 
    * Returns:
    ***********************************************************************************************/  
    template<typename RMAE, class T>
    _ROWMAT_SCALAR_PROD_<RMAE,T>::_ROWMAT_SCALAR_PROD_(RMAE const& mat, T val) : mat_(mat), 
                                                                                 val_(val) {}




   /***********************************************************************************************
    * dtype _ROWMAT_SCALAR_PROD_::operator() (index_t i) const
    * 
    * Retrieve evaluated (single) indexed value
    * 
    * Template Parameters:
    *     RMAE     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     T        (implicit) Scalar data type
    * 
    * Parameters:
    *     i       (index_t) Element index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix 
    ***********************************************************************************************/ 
    template<typename RMAE, class T>
    typename _ROWMAT_SCALAR_PROD_<RMAE,T>::dtype 
        _ROWMAT_SCALAR_PROD_<RMAE,T>::operator() (index_t i) const
    {return this->mat_(i) * this->val_;}




   /***********************************************************************************************
    * dtype _ROWMAT_SCALAR_PROD_::operator() (index_t i, index_t j) const
    *
    * Retrieve evaluated (double) indexed value 
    * 
    * Template Parameters:
    *     RMAE     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     T        (implicit) Scalar data type
    * 
    * Parameters:
    *     i       (index_t) Row index
    *     j       (index_t) Column index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix  
    ***********************************************************************************************/ 
    template<typename RMAE, class T>
    typename _ROWMAT_SCALAR_PROD_<RMAE,T>::dtype 
        _ROWMAT_SCALAR_PROD_<RMAE,T>::operator() (index_t i, index_t j) const
    {return this->mat_(i,j) * this->val_;} 





   /**********************************************************************************************
    * void _ROWMAT_SCALAR_PROD_::eval(Matrix& target) const
    * 
    * Evaluates the expression to target Matrix of any kind.
    * 
    * Template Parameters:
    *     RMAE     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     T        (implicit) Scalar data type
    *     MAT_T     (implicit) (Matrix) Matrix type 
    * 
    * Parameters:
    *     target    (Matrix*) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    template<typename MAT_T>
    void _ROWMAT_SCALAR_PROD_<RMAE1,RMAE2>::eval(Matrix<MAT_T>& target) const
    {
        std::cout << "'void _ROWMAT_SUM_<RMAE1,RMAE2>::eval(Matrix<MAT_T>& target) const' called" << std::endl;
        const index_t block = this->block_size();
        const index_t M = this->rows();
        const index_t N = this->cols();
        
        #ifdef DIM_CHECK // Check that target matrix is correctly sized
        if (  (this->rows() != target.rows()) || (this->cols() != target.cols()) 
            &&(this->size()!=0) )
        {
            std::cerr << "Error in 'void _ROWMAT_SCALAR_PROD_<RMAE,T>::eval(Matrix<MAT_T>* target)"
                      << "\n    target Matrix with dimensions " << target->rows() << "x" 
                      << target->cols() << " incompatible with _ROWMAT_SUM_ with dimensions " 
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
    * void _ROWMAT_SCALAR_PROD_::eval(RowMatrix&) const
    * 
    * Evaluates the expression to target RowMatrix.
    * 
    * Template Parameters:
    *     RMAE     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     T1       (implicit) Scalar data type
    *     T2       (implicit) (Matrix) Target RowMatrix data type 
    * 
    * Parameters:
    *     target    (RowMatrix) Target RowMatrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename RMAE, typename T1>
    template<class T2>
    void _ROWMAT_SCALAR_PROD_<RMAE,T1>::eval(RowMatrix<T2>& target) const
    {
        const index_t istop = this->rows()+1;
        const index_t jstop = this->cols()+1;
        for (index_t i=1; i!=istop; ++i)
            for (index_t j=1; j!=jstop; ++j)
                target(i,j) = (*this)(i,j);
        return;
    }




   /**********************************************************************************************
    * index_t rows() const
    * 
    * Gets number of rows in the expression 
    * 
    * Template Parameters:
    *     RMAE       RowMatrix or _ROWMAT_ATOMIC_EXPR_ 
    *     T          Scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t    Number of rows in expression
    **********************************************************************************************/
    template<typename RMAE, typename T>
    index_t _ROWMAT_SCALAR_PROD_<RMAE,T>::rows() const
    {return this->mat_.rows();}





   /**********************************************************************************************
    * index_t cols() const
    * 
    * Gets number of columns in the expression 
    * 
    * Template Parameters:
    *     RMAE       RowMatrix or _ROWMAT_ATOMIC_EXPR_ 
    *     T          Scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t    Number of columns in expression
    **********************************************************************************************/
    template<typename RMAE, typename T>
    index_t _ROWMAT_SCALAR_PROD_<RMAE,T>::cols() const 
    {return this->mat_.cols();} 




   /**********************************************************************************************
    * index_t size() const
    * 
    * Gets number of rows in the expression 
    * 
    * Template Parameters:
    *     RMAE       RowMatrix or _ROWMAT_ATOMIC_EXPR_ 
    *     T          Scalar data type
    * 
    * Parameters:
    * Returns:
    *     index_t    Size of expression
    **********************************************************************************************/
    template<typename RMAE, typename T>
    index_t _ROWMAT_SCALAR_PROD_<RMAE,T>::size() const
    {return this->mat_.size();}






   /***********************************************************************************************
    * index_t block_size()
    *
    * Computes a cache block size based on current thread's chunk size
    ***********************************************************************************************/ 
    template<typename RMAE, class T>
    index_t _ROWMAT_SCALAR_PROD_<RMAE,T>::block_size()
    {
        index_t type_size = sizeof(_ROWMAT_SUM_<RMAE,T>::dtype);
        return chunk_size()/type_size;
    }  
    

}


//#endif