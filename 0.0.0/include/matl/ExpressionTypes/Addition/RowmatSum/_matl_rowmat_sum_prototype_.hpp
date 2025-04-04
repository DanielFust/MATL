//#ifndef _MATL_ROWMAT_SUM_PROTOTYPE__HPP
//#define _MATL_ROWMAT_SUM_PROTOTYPE__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method delarations and implementations of Row Matrix Summations as well
 * as operators to construct Row Matrix Summation objects
 **************************************************************************************************/
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"

namespace matl
{
   /***********************************************************************************************
    * class _ROWMAT_SUM_
    * 
    * Summation of atomic RowMatrix objects or expressions  
    * 
    * Template Arguments:
    *     RMAE1     (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    *     RMAE2     (_ROWMAT_ATOMIC_EXPR_) RowMatrix or row major atomic expression
    ***********************************************************************************************/ 
    template<typename RMAE1, typename RMAE2>
    class _ROWMAT_SUM_ : public _MATEXPRESSION_< _ROWMAT_SUM_<RMAE1,RMAE2> >
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
        _ROWMAT_SUM_(RMAE1 const& mat1, RMAE2 const& mat2);

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
        index_t rows() const {return this->mat1_.rows();}
        
        // Number of columns in resultant matrix
        index_t cols() const {return this->mat1_.cols();}

        // Total size of resultant matrix
        index_t size() const {return this->mat1_.size();}

    };

}


//#endif