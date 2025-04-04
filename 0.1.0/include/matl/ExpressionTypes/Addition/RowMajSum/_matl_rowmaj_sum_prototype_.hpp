#ifndef _MATL_ROWMAJ_SUM_PROTOTYPE__HPP
#define _MATL_ROWMAJ_SUM_PROTOTYPE__HPP
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
    * class _ROWMAJ_SUM_
    * 
    * Summation of RowMajorMatrix objects or expressions  
    * 
    * Template Arguments:
    *     RME1      RowMajorMatrix or row major expression
    *     RME2      RowMajorMatrix or row major expression
    ***********************************************************************************************/ 
    template<typename RMA1, typename RMA2>
    class _ROWMAJ_SUM_ : public _MATEXPRESSION_< _ROWMAJ_SUM_<RME1,RME2> >
    {
      private:
        RME1 const& mat1_;
        RME2 const& mat2_;

        index_t block_size();

      public:
        // Return Type
        typedef typename promote<typename TypeTraits<RME1>::dtype, 
                                 typename TypeTraits<RME2>::dtype>::type dtype;

        // Constructor
        _ROWMAJ_SUM_(RME1 const& mat1, RME2 const& mat2);

        // Compute and return (single) indexed value
        dtype operator() (index_t) const;

        // Compute and return (double) indexed value
        dtype operator() (index_t, index_t) const;

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(MATRIX<MAT_T>& target) const;

        // Evaluate expression to target RowMajorMatrix
        template<typename T>
        void eval(RowMajorMatrix<T>& target) const;

        // Number of rows in resultant matrix
        index_t rows() const {return this->mat1_.rows();}
        
        // Number of columns in resultant matrix
        index_t cols() const {return this->mat1_.cols();}

        // Total size of resultant matrix
        index_t size() const {return this->mat1_.size();}

    };

}


#endif