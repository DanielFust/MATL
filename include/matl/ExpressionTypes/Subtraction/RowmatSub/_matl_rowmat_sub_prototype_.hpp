//#ifndef _BLAS_ROWMAT_SUB_PROTOTYPE__HPP
//#define _BLAS_ROWMAT_SUB_PROTOTYPE__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method delarations of Row Matrix Subtractions
 **************************************************************************************************/
//#include "_blas_forward_declarations_.hpp"
//#include "_blas_type_traits_.hpp"
//#include "_blas_utils_.hpp"

namespace matl
{
   /***********************************************************************************************
    * class _ROWMAT_SUB_
    * 
    * Subtraction of atomic RowMatrix objects or expressions  
    * 
    * Template Arguments:
    *     RMAE1     (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     RMAE2     (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    ***********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    class _ROWMAT_SUB_ : public _ROWMAT_ATOMIC_EXPR_< _ROWMAT_SUB_<RMAE1,RMAE2> >
    {
      private:
        RMAE1 const& mat1_;
        RMAE2 const& mat2_;

        index_t block_size();

      public:
        // Return Type
        typedef typename promote<typename TypeTraits<RMAE1>::dtype, 
                                 typename TypeTraits<RMAE2>::dtype>::type dtype;

        // Constructor
        _ROWMAT_SUB_(RMAE1 const& mat1, RMAE2 const& mat2);

        // Compute and return (single) indexed value
        dtype operator() (index_t) const;

        // Compute and return (double) indexed value
        dtype operator() (index_t, index_t) const;

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const;

        // Evaluate expression to target RowMatrix
        template<typename T>
        void eval(RowMatrix<T>& target) const;

        // Number of rows in resultant matrix
        index_t rows() const;
        
        // Number of columns in resultant matrix
        index_t cols() const;

        // Total size of resultant matrix
        index_t size() const;

    };

}


//#endif