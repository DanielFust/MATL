//#ifndef _MATL_ROWMAT_SUM_IMPLEMENTATION__HPP
//#define _MATL_ROWMAT_SUM_IMPLEMENTATION__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method delarations and implementations of Row Matrix Summations as well
 * as operators to construct Row Matrix Summation objects
 **************************************************************************************************/
//#include "_matl_forward_declr_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"

namespace matl
{


   /***********************************************************************************************
    * _ROWMAT_SUM_(RMAE1 const& mat1, RMAE2 const& mat2)
    * 
    * Constructor
    * 
    * Template Parameters:
    *     RMAE1     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     RMAE2     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    * 
    * Parameters:
    *     mat1      RowMatrix or row matrix/major atomic expression
    *     mat2      RowMatrix or row matrix/major atomic expression
    * 
    * Returns:
    ***********************************************************************************************/
    template<typename RMAE1, typename RMAE2>
    _ROWMAT_SUM_<RMAE1,RMAE2>::_ROWMAT_SUM_(RMAE1 const& mat1, RMAE2 const& mat2) : 
        mat1_(mat1), mat2_(mat2) {}






   /***********************************************************************************************
    * dtype _ROWMAT_SUM_::operator() (index_t i) const
    * 
    * Retrieve evaluated (single) indexed value
    * 
    * Template Parameters:
    *     RMAE1     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     RMAE2     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    * 
    * Parameters:
    *     i       (index_t) Element index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix 
    ***********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    typename _ROWMAT_SUM_<RMAE1,RMAE2>::dtype _ROWMAT_SUM_<RMAE1,RMAE2>::operator() (index_t i) const
    {return this->mat1_(i) + this->mat2_(i);}


   /***********************************************************************************************
    * dtype _ROWMAT_SUM_::operator() (index_t i, index_t j) const
    *
    * Retrieve evaluated (double) indexed value 
    * 
    * Template Parameters:
    *     RMAE1     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     RMAE2     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    * 
    * Parameters:
    *     i       (index_t) Row index
    *     j       (index_t) Column index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix  
    ***********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    typename _ROWMAT_SUM_<RMAE1,RMAE2>::dtype 
        _ROWMAT_SUM_<RMAE1,RMAE2>::operator() (index_t i, index_t j) const
    {return this->mat1_(i,j) + this->mat2_(i,j);} 



    // DEPRECIATED
   /***********************************************************************************************
    * RowMatrix<dtype> _ROWMAT_SUM_::eval() const
    *
    * Evaluate the summation, creating a new RowMatrix.
    * 
    * Template Parameters:
    *     RMAE1     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     RMAE2     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    * 
    * Parameters:
    * 
    * Returns:
    *     RowMatrix    New RowMatrix object resulting from the evaluation of all elements
    ***********************************************************************************************/ 
    /*
    template<typename RMAE1, typename RMAE2>
    RowMatrix<typename _ROWMAT_SUM_<RMAE1,RMAE2>::dtype> _ROWMAT_SUM_<RMAE1,RMAE2>::eval() const
    {
        RowMatrix<dtype> row_matrix(this->rows(), this->cols());
        const index_t N = this->size();
        
        for (index_t i = 1; i != (N+1); ++i) row_matrix(i) = (*this)(i);
        return row_matrix;
    }
    */




   /**********************************************************************************************
    * void _ROWMAT_SUM_::eval(Matrix& target) const
    * 
    * Evaluates the expression to target Matrix of any kind.
    * 
    * Template Parameters:
    *     RMAE1     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     RMAE2     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     MAT_T     (implicit) (Matrix) Matrix type 
    * 
    * Parameters:
    *     target    (Matrix*) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    template<typename MAT_T>
    void _ROWMAT_SUM_<RMAE1,RMAE2>::eval(Matrix<MAT_T>& target) const
    {
        std::cout << "'void _ROWMAT_SUM_<RMAE1,RMAE2>::eval(Matrix<MAT_T>& target) const' called" << std::endl;
        const index_t block = this->block_size();
        const index_t M = this->rows();
        const index_t N = this->cols();
        
        #ifdef MATL_DIM_CHECK // Check that target matrix is correctly sized
        if (  (this->rows() != target.rows()) || (this->cols() != target.cols()) 
            &&(this->size()!=0) )
        {
            std::cerr << "Error in 'void _ROWMAT_SUM_<E1,E2>::eval(Matrix<MAT_T>* target)"
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
    * void _ROWMAT_SUM_::eval(RowMatrix&) const
    * 
    * Evaluates the expression to target RowMatrix.
    * 
    * Template Parameters:
    *     RMAE1     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     RMAE2     (implicit) (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     T2        (implicit) (Matrix) Target RowMatrix data type 
    * 
    * Parameters:
    *     target    (RowMatrix) Target RowMatrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    template<class T>
    void _ROWMAT_SUM_<RMAE1,RMAE2>::eval(RowMatrix<T>& target) const
    {
        const index_t istop = this->rows()+1;
        const index_t jstop = this->cols()+1;
        for (index_t i=1; i!=istop; ++i)
            for (index_t j=1; j!=jstop; ++j)
                target(i,j) = (*this)(i,j);
        return;
    }







   /***********************************************************************************************
    * index_t block_size()
    *
    * Computes a cache block size based on current thread's chunk size
    ***********************************************************************************************/ 
    template<typename RMAE_T1, typename RMAE_T2>
    index_t _ROWMAT_SUM_<RMAE_T1,RMAE_T2>::block_size()
    {
        index_t type_size = sizeof(_ROWMAT_SUM_<RMAE_T1,RMAE_T2>::dtype);
        return chunk_size()/type_size;
    }  
    

}


//#endif