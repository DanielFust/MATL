//#ifndef _MATL_GENERIC_MAT_SCALAR_SUM_PROTOTYPE__HPP
//#define _MATL_GENERIC_MAT_SCALAR_SUM_PROTOTYPE__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method delarations and implementations of Row Matrix Scalar 
 * multiplications
 **************************************************************************************************/
#include <iostream>
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"

namespace matl
{
   /***********************************************************************************************
    * class _GENERIC_MAT_SCALAR_SUM_
    * 
    * Addition of Matrix objects or expressions with scalars 
    * 
    * Template Arguments:
    *     AE_T     (_ATOMIC_EXPR_) Matrix or atomic matrix expression
    *     T        Scalar data type
    ***********************************************************************************************/ 
    template<typename AE_T, class T>
    class _GENERIC_MAT_SCALAR_SUM_ : public _MATEXPRESSION_< _GENERIC_MAT_SCALAR_SUM_<AE_T,T> >
    {
      private:
        AE_T const& mat_;
        T    const  val_;

        index_t block_size() const;

      public:
        // Return Type
        typedef typename TypeTraits<AE_T>::dtype dtype;

        // Constructor
        _GENERIC_MAT_SCALAR_SUM_(AE_T const& mat, T val);
        ~_GENERIC_MAT_SCALAR_SUM_(){std::cout<<"_GENERIC_MAT_SCALAR_SUM_ Destructor called"<<std::endl;}

        // Compute and return (single) indexed value
        dtype operator() (index_t) const;

        // Compute and return (double) indexed value
        dtype operator() (index_t, index_t) const;

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const;

        // Number of rows in resultant matrix
        index_t rows() const;
        
        // Number of columns in resultant matrix
        index_t cols() const;

        // Total size of resultant matrix
        index_t size() const;

        // Check if there is memory conflict with target matrix
        template<typename MAT_T>
        bool memory_conflict(Matrix<MAT_T> const&) const;

        // Is the expression evaluatable element-by-element
        constexpr bool elementwise() {return true;}

    };    

}


//#endif