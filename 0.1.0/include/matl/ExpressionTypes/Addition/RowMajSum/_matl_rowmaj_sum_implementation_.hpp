#ifndef _MATL_ROWMAJ_SUM_IMPLEMENTATION__HPP
#define _MATL_ROWMAJ_SUM_IMPLEMENTATION__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method delarations and implementations of Row Major Summations as well
 * as operators to construct Row Major Summation objects
 **************************************************************************************************/
//#include "_matl_forward_declr_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"

namespace matl
{

   /***********************************************************************************************
    * _ROWMAJ_SUM_(RME1 const& mat1, RME2 const& mat2)
    * 
    * Constructor
    * 
    * Template Parameters:
    *     RME1     (implicit)  RowMajorMatrix or row major expression
    *     RME2     (implicit)  RowMajorMatrix or row major expression
    * 
    * Parameters:
    *     mat1      RowMajorMatrix or row matrix/major expression
    *     mat2      RowMajorMatrix or row matrix/major expression
    * 
    * Returns:
    ***********************************************************************************************/
    template<typename RME1, typename RME2>
    _ROWMAJ_SUM_<RME1,RME2>::_ROWMAJ_SUM_(RME1 const& mat1, RME2 const& mat2) : 
        mat1_(mat1), mat2_(mat2) {}






   /***********************************************************************************************
    * dtype _ROWMAJ_SUM_::operator() (index_t i) const
    * 
    * Retrieve evaluated (single) indexed value
    * 
    * Template Parameters:
    *     RME1     (implicit)  RowMajorMatrix or row major expression
    *     RME2     (implicit)  RowMajorMatrix or row major expression
    * 
    * Parameters:
    *     i       (index_t) Element index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix 
    ***********************************************************************************************/ 
    template<typename RME1, typename RME2>
    typename _ROWMAJ_SUM_<RME1,RME2>::dtype _ROWMAJ_SUM_<RME1,RME2>::operator() (index_t i) const
    {return this->mat1_(i) + this->mat2_(i);}


   /***********************************************************************************************
    * dtype _ROWMAJ_SUM_::operator() (index_t i, index_t j) const
    *
    * Retrieve evaluated (double) indexed value 
    * 
    * Template Parameters:
    *     RME1     (implicit)  RowMajorMatrix or row major expression
    *     RME2     (implicit)  RowMajorMatrix or row major expression
    * 
    * Parameters:
    *     i       (index_t) Row index
    *     j       (index_t) Column index
    * 
    * Returns:
    *     dtype   Data type of resultant matrix  
    ***********************************************************************************************/ 
    template<typename RME1, typename RME2>
    typename _ROWMAJ_SUM_<RME1,RME2>::dtype 
        _ROWMAJ_SUM_<RME1,RME2>::operator() (index_t i, index_t j) const
    {return this->mat1_(i,j) + this->mat2_(i,j);} 



    // DEPRECIATED
   /***********************************************************************************************
    * RowMajorMatrix<dtype> _ROWMAJ_SUM_::eval() const
    *
    * Evaluate the summation, creating a new RowMajorMatrix.
    * 
    * Template Parameters:
    *     RME1     (implicit)  RowMajorMatrix or row major expression
    *     RME2     (implicit)  RowMajorMatrix or row major expression
    * 
    * Parameters:
    * 
    * Returns:
    *     RowMajorMatrix    New RowMajorMatrix object resulting from the evaluation of all elements
    ***********************************************************************************************/ 
    /*
    template<typename RME1, typename RME2>
    RowMajorMatrix<typename _ROWMAJ_SUM_<RME1,RME2>::dtype> _ROWMAJ_SUM_<RME1,RME2>::eval() const
    {
        RowMajorMatrix<dtype> row_matrix(this->rows(), this->cols());
        const index_t N = this->size();
        
        for (index_t i = 1; i != (N+1); ++i) row_matrix(i) = (*this)(i);
        return row_matrix;
    }
    */




   /**********************************************************************************************
    * void _ROWMAJ_SUM_::eval(Matrix& target) const
    * 
    * Evaluates the expression to target Matrix of any kind.
    * 
    * Template Parameters:
    *     RME1     (implicit)  RowMajorMatrix or row major expression
    *     RME2     (implicit)  RowMajorMatrix or row major expression
    *     MAT_T    (implicit) (Matrix) Matrix type 
    * 
    * Parameters:
    *     target    (Matrix*) Target matrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename RME1, typename RME2>
    template<typename MAT_T>
    void _ROWMAJ_SUM_<RME1,RME2>::eval(Matrix<MAT_T>& target) const
    {
        std::cout << "'void _ROWMAJ_SUM_<RME1,RME2>::eval(Matrix<MAT_T>& target) const' called" << std::endl;
        const index_t block = this->block_size();
        const index_t M = this->rows();
        const index_t N = this->cols();
        
        #ifdef MATL_DIM_CHECK // Check that target matrix is correctly sized
        if (  (this->rows() != target.rows()) || (this->cols() != target.cols()) 
            &&(this->size()!=0) )
        {
            std::cerr << "Error in 'void _ROWMAJ_SUM_<E1,E2>::eval(Matrix<MAT_T>* target)"
                      << "\n    target Matrix with dimensions " << target->rows() << "x" 
                      << target->cols() << " incompatible with _ROWMAJ_SUM_ with dimensions " 
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
    * void _ROWMAJ_SUM_::eval(RowMajorMatrix&) const
    * 
    * Evaluates the expression to target RowMajorMatrix.
    * 
    * Template Parameters:
    *     RME1     (implicit)  RowMajorMatrix or row major expression
    *     RME2     (implicit)  RowMajorMatrix or row major expression
    *     T2        (implicit) (Matrix) Target RowMajorMatrix data type 
    * 
    * Parameters:
    *     target    (RowMajorMatrix) Target RowMajorMatrix for expression evaluation
    * 
    * Returns:
    **********************************************************************************************/ 
    template<typename RME1, typename RME2>
    template<class T>
    void _ROWMAJ_SUM_<RME1,RME2>::eval(RowMajorMatrix<T>& target) const
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
    template<typename RME_T1, typename RME_T2>
    index_t _ROWMAJ_SUM_<RME_T1,RME_T2>::block_size()
    {
        index_t type_size = sizeof(_ROWMAJ_SUM_<RME_T1,RME_T2>::dtype);
        return chunk_size()/type_size;
    }  
    

}


#endif