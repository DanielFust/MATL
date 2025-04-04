//#ifndef _MATL_ROWMAT_SCALAR_SUM_PROTOTYPE__HPP
//#define _MATL_ROWMAT_SCALAR_SUM_PROTOTYPE__HPP
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
    * class _ROWMAT_SCALAR_SUM_
    * 
    * Addition of RowMatrix objects or expressions with scalars 
    * 
    * Template Arguments:
    *     RMAE     (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     T        Scalar data type
    ***********************************************************************************************/ 
    template<typename RMAE, class T>
    class _ROWMAT_SCALAR_SUM_ : public _MATEXPRESSION_< _ROWMAT_SCALAR_SUM_<RMAE,T> >
    {
      private:
        RMAE const& mat_;
        T    const  val_;

        index_t block_size();

      public:
        // Return Type
        typedef typename promote<typename TypeTraits<RMAE>::dtype, T>::type dtype;

        // Constructor
        _ROWMAT_SCALAR_SUM_(RMAE const& mat, T val);

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

}


//#endif