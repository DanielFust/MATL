#ifndef _MATL_ROWMAJEXPRESSION_PROTOTYPE__HPP
#define _MATL_ROWMAJEXPRESSION_PROTOTYPE__HPP

//#include <iostream>
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"

namespace matl
{
   /**********************************************************************************************
    * class _ROWMAJEXPRESSION_ (Prototype)
    * 
    * Intermediate derived type to identify row major matrix expressions
    * 
    * Template Parameters:
    *     RME_T    (_ROWMAJEXPRESSION_) Expression type
    **********************************************************************************************/
    template<typename RME_T>
    class _ROWMAJEXPRESSION_ : public _MATEXPRESSION_<_ROWMAJEXPRESSION_<RME_T> >
    {
      public:

        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename TypeTraits<RME_T>::dtype dtype;

        // Preferred matrix return type
        typedef typename TypeTraits<RME_T>::matrix_type matrix_type;

        // Evaluate Expression to target matrix
        template<typename MAT_T> void eval(MATRIX<MAT_T>& target) const;

        // Evaluate Expression to new matrix
        //MATRIX<matrix_type> eval() const;

        // Get indexed result of expression
        dtype operator() (index_t, index_t) const;

        // Get number of rows of resultant matrix
        index_t rows() const;

        // Get number of columns of resultant matrix
        index_t cols() const;

        // Get number of elements in resultant matrix
        index_t size() const;

        // Check if there is memory conflict with target matrix
        //template<typename MAT_T>
        //bool memory_conflict(MATRIX<MAT_T> const&) const;

        // Is the expression evaluatable element-by-element
        constexpr bool elementwise();
              
    };
}

#endif