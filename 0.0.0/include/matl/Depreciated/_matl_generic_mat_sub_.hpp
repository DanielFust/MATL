#ifndef _MATL_GENERIC_MAT_SUB__HPP
#define _MATL_GENERIC_MAT_SUB__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 **************************************************************************************************/
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"
#include "_matl_utils_.hpp"

namespace matl
{
   /***********************************************************************************************
    * class _GENERIC_MAT_SUB_
    * 
    * Template expression for the subtraction of generic atomic matrix expressions. This expression
    * template is the default representation for matrix subtraction. This expression template 
    * supports lazy execution and some cache blocking, however more specific types are recommended
    * to take advantage of specific data storage architectures of Matrix types.
    * 
    * Template Arguments:
    *     AE_T1     (_ATOMIC_EXPR_) Atomic matrix expression
    *     AE_T2     (_ATOMIC_EXPR_) Atomic matrix expression
    ***********************************************************************************************/ 
    template<typename AE_T1, typename AE_T2>
    class _GENERIC_MAT_SUB_ : public _ATOMIC_EXPR_< _GENERIC_MAT_SUB_<AE_T1,AE_T2> >
    {
      private:
        AE_T1 const& mat1_;
        AE_T2 const& mat2_;

        index_t block_size();

      public:
        // Return Type
        typedef typename promote<typename TypeTraits<AE_T1>::dtype, 
                                 typename TypeTraits<AE_T2>::dtype>::type dtype;

        // Constructor
        _GENERIC_MAT_SUB_(AE_T1 const& mat1, AE_T2 const& mat2) : mat1_(mat1), mat2_(mat2) {}

        // Compute and return (double) indexed value
        dtype operator() (index_t, index_t) const;

        // Evaluate expression 
        //Matrix<typename TypeTraits<AE_T1>::mat_type> eval() const;

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const;

        // Number of rows in resultant matrix
        index_t rows() const {return mat1_.rows();}
        
        // Number of columns in resultant matrix
        index_t cols() const {return mat1_.cols();}

        // Total size of resultant matrix
        index_t size() const {return mat1_.size();}

    };







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

        #ifdef DIM_CHECK // Check that target matrix is correctly sized
        if (  ((this->rows_ != target->rows()) || (this->rows_ != 0))
            ||((this->cols_ != target->cols()) || (this->cols_ != 0)) )
        {
            std::cerr << "Error in 'void _GENERIC_MAT_SUB_<E1,E2>::eval(Matrix<MAT_T>* target)"
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
    template<typename AE_T1, typename AE_T2>
    index_t _GENERIC_MAT_SUB_<AE_T1,AE_T2>::block_size()
    {
        index_t type_size = sizeof(_GENERIC_MAT_SUB_<AE_T1,AE_T2>::dtype);
        return chunk_size()/type_size;
    }  


}
#endif