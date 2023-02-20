#ifndef _MATL_GENERIC_MAT_SCALAR_PROD__HPP
#define _MATL_GENERIC_MAT_SCALAR_PROD__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 **************************************************************************************************/
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"
#include "_matl_utils_.hpp"

namespace matl
{
   /***********************************************************************************************
    * class _GENERIC_MAT_SCALAR_PROD_
    * 
    * Template expression for the scalar multiplication of generic atomic matrix expressions. This
    * expression template is the default representation for matrix scalar multiplication. 
    * 
    * Template Arguments:
    *     AE        (_ATOMIC_EXPR_ or Matrix) Atomic matrix expression or matrix
    *     T         Scalar data type
    ***********************************************************************************************/ 
    template<typename AE, class T>
    class _GENERIC_MAT_SCALAR_PROD_ : public _ATOMIC_EXPR_< _GENERIC_MAT_SCALAR_PROD_<AE,T> >
    {
      private:
        AE const& mat_;
        T         val_;

        index_t block_size();

      public:
        // Return Type
        typedef typename promote<typename TypeTraits<AE>::dtype, T>::type dtype;

        // Constructor
        _GENERIC_MAT_SCALAR_PROD_(AE const& mat, T val) : mat_(mat), val_(val) {}

        // Compute and return (double) indexed value
        dtype operator() (index_t, index_t) const;

        // Evaluate expression 
        //Matrix<typename TypeTraits<AE_T1>::mat_type> eval() const;

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const;

        // Number of rows in resultant matrix
        index_t rows() const {return mat_.rows();}
        
        // Number of columns in resultant matrix
        index_t cols() const {return mat_.cols();}

        // Total size of resultant matrix
        index_t size() const {return mat_.size();}

    };







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

        #ifdef DIM_CHECK // Check that target matrix is correctly sized
        if (  ((this->rows() != target->rows()) || (this->rows_ != 0))
            ||((this->cols() != target->cols()) || (this->cols_ != 0)) )
        {
            std::cerr << "Error in 'void _GENERIC_MAT_SCALAR_PROD_<AE,T>::eval(Matrix<MAT_T>& target)"
                    << "\n    target Matrix with dimensions " << target->rows() << "x" 
                    << target->cols() << " incompatible with _GENERIC_MAT_SCALAR_PROD_ with " 
                    << "dimensions " << this->rows() << "x" << this->cols() 
                    << std::endl;
            abort();          
        }    
        #endif

        // Cache blocking to mitigate miss-rate for intermingled row and column major ordering
        for (index_t ii=1; ii<(M+1); ii+=block)
            for (index_t jj=1; jj<(N+1); jj+=block)
            {
                const index_t istop = branchless_min(ii+block, M);
                const index_t jstop = branchless_min(jj+block, N);
                for (index_t i=ii; i!=istop; ++i)
                    for (index_t j=jj; j!=jstop; ++j)
                        target(i,j) = (*this)(i,j);
            }
  
        return;
    }
    



    
   /***********************************************************************************************
    * index_t block_size()
    *
    * Computes a cache block size based on current thread's chunk size
    ***********************************************************************************************/ 
    template<typename AE, class T>
    index_t _GENERIC_MAT_SCALAR_PROD_<AE,T>::block_size()
    {
        index_t type_size = sizeof(_GENERIC_MAT_SCALAR_PROD_<AE,T>::dtype);
        return chunk_size()/type_size;
    }  


}
#endif