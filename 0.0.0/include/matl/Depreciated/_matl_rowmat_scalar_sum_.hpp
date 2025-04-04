#ifndef _BLAS_ROWMAT_SCALAR_SUM__HPP
#define _BLAS_ROWMAT_SCALAR_SUM__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method delarations and implementations of Row Matrix Scalar 
 * multiplications
 **************************************************************************************************/
#include "_blas_forward_declr_.hpp"
#include "_blas_type_traits_.hpp"
#include "_blas_utils_.hpp"

namespace blas
{
   /***********************************************************************************************
    * class _ROWMAT_SCALAR_SUM_
    * 
    * Addition of RowMatrix objects or expressions with scalars 
    * 
    * Template Arguments:
    *     RMAE     (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     T        Scalar data type
    ***********************************************************************************************/ 
    template<typename RMAE, class T>
    class _ROWMAT_SCALAR_SUM_ : public _ROWMAT_ATOMIC_EXPR_< _ROWMAT_SCALAR_SUM_<RMAE,T> >
    {
      private:
        RMAE const& mat_;
        T    const  val_;

        index_t block_size();

      public:
        // Return Type
        typedef typename promote<typename TypeTraits<RMAE>::dtype, T>::type dtype;

        // Constructor
        _ROWMAT_SCALAR_SUM_(RMAE const& mat, T val) : mat_(mat), val_(val) {}

        // Compute and return (single) indexed value
        dtype operator() (index_t) const;

        // Compute and return (double) indexed value
        dtype operator() (index_t, index_t) const;

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const;

        // Evaluate expression to target RowMatrix
        template<typename T2>
        void eval(RowMatrix<T2>& target) const;

        // Number of rows in resultant matrix
        index_t rows() const {return this->mat_.rows();}
        
        // Number of columns in resultant matrix
        index_t cols() const {return this->mat_.cols();}

        // Total size of resultant matrix
        index_t size() const {return this->mat_.size();}

    };





   /***********************************************************************************************
    * dtype _ROWMAT_SCALAR_SUM_::operator() (index_t i) const
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
    typename _ROWMAT_SCALAR_SUM_<RMAE,T>::dtype 
        _ROWMAT_SCALAR_SUM_<RMAE,T>::operator() (index_t i) const
    {return this->mat_(i) + this->val_;}




   /***********************************************************************************************
    * dtype _ROWMAT_SCALAR_SUM_::operator() (index_t i, index_t j) const
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
    typename _ROWMAT_SCALAR_SUM_<RMAE,T>::dtype 
        _ROWMAT_SCALAR_SUM_<RMAE,T>::operator() (index_t i, index_t j) const
    {return this->mat_(i,j) + this->val_;} 





   /**********************************************************************************************
    * void _ROWMAT_SCALAR_SUM_::eval(Matrix& target) const
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
    void _ROWMAT_SCALAR_SUM_<RMAE1,RMAE2>::eval(Matrix<MAT_T>& target) const
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
    * void _ROWMAT_SCALAR_SUM_::eval(RowMatrix&) const
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
    void _ROWMAT_SCALAR_SUM_<RMAE,T1>::eval(RowMatrix<T2>& target) const
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
    template<typename RMAE, class T>
    index_t _ROWMAT_SCALAR_SUM_<RMAE,T>::block_size()
    {
        index_t type_size = sizeof(_ROWMAT_SUM_<RMAE,T>::dtype);
        return chunk_size()/type_size;
    }  
    

}


#endif