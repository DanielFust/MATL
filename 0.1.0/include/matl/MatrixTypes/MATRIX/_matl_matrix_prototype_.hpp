#ifndef _MATL_MATRIX__HPP
#define _MATL_MATRIX__HPP

//#include <iostream>
//#include "../../_matl_forward_declarations_.hpp"
//#include "../../Utilities/_matl_type_traits_.hpp"


namespace matl 
{
   /**********************************************************************************************
    * class MATRIX
    * 
    * Base class for 'matl' matrix types.
    * 
    * Template Parameters:
    *     MAT_T       (Matrix) Matrix type (e.g. RowMatrix, ColMatrix, CSRMatrix, etc.)
    **********************************************************************************************/ 
    template<typename MAT_T>
    class MATRIX
    {
      public:  

        // Matrix datatype
        typedef typename TypeTraits<MAT_T>::dtype dtype;

        /****************** Sized Constructor ******************/
        // I'm not sure if this can exist, but I need something that can create both statically-sized and dynamically-sized matrices
        //template<index_t M, index_t N>
        //MATRIX<MAT_T> MATRIX()

        /**************** Assignment & Copying *****************/
        
        // Matrix copy
        template<typename MAT_T2>
        void operator= (MATRIX<MAT_T2>& mat);

        // Expression evaluation trigger
        template<typename EXPR_T>
        void operator= (_MATEXPRESSION_<EXPR_T>& mat);

        // Creates a blank copy of the matrix with same 'size'   <-- hopefully get rid of this?
        //MATRIX<MAT_T> blank_copy() const;

        /********************* Data Access *********************/

        // Single index getter
        dtype operator() (index_t i) const;

        // Double index getter
        dtype operator() (index_t i, index_t j) const;

        // Single index accessor
        dtype& operator() (index_t i);

        // Double index accessor
        dtype& operator() (index_t i, index_t j);

        /********************* Utilities *********************/

        // Scalar data writing
        template<class T, enable_if_scalar<T> = true>
        void operator= (T val);

        // Row Getter
        index_t rows() const;

        // Coloumn Getter
        index_t cols() const;

        // Size Getter
        index_t size() const;

        // Buffer size
        index_t buffer() const;

        // Shaped Resize
        //void resize(index_t, index_t);  <-- get rid of this? (should still be implemented in dynamically-sized derived types)

        void println() const;

    };


}



#endif