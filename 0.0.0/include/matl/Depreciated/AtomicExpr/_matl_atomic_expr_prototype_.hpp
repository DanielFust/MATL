//#ifndef _MATL_ATOMIC_EXPR_PROTOTYPE__HPP
//#define _MATL_ATOMIC_EXPR_PROTOTYPE__HPP

//#include <iostream>
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"

namespace matl
{
   
   /**********************************************************************************************
    * class _ATOMIC_EXPR_ Prototype
    * 
    * Base class for atomic matrix expression templates. Atomic expressions are operations that
    * can be performed independently on all elements of a matrix concurrently.
    * 
    * Template Parameters:
    *     EXPR_T         (_ATOMIC_EXPR_) Expression type
    **********************************************************************************************/
    template<typename EXPR_T>
    class _ATOMIC_EXPR_ : public _MATEXPRESSION_< _ATOMIC_EXPR_<EXPR_T> >
    {
      public:

        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename TypeTraits<EXPR_T>::dtype dtype;

        // Resultant matrix type (e.g. row matrix, compressed sparse column, etc.)
        //typedef typename TypeTraits<EXPR_T>::mat_type mat_type;

        // Single index data getter
        dtype operator() (index_t i) const;

        // Double index data getter
        dtype operator() (index_t i, index_t j) const;

        // Evaluate expression
        //mat_type eval() const
        //{return static_cast<EXPR_T const&>(*this).eval();}

        // Evaluate Expression to target matrix
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const;

        // Get number of rows of resultant matrix
        index_t rows() const;

        // Get number of columns of resultant matrix
        index_t cols() const;

        // Get number of elements in resultant matrix
        index_t size() const;

        // Check if there is memory conflict with target matrix
        template<typename MAT_T>
        bool memory_conflict(Matrix<MAT_T> const&) const;
              
    };
}

//#endif