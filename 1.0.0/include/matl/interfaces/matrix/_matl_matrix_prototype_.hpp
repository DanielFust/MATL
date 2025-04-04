#ifndef _MATL_MATRIX_PROTOTYPE__HPP
#define _MATL_MATRIX_PROTOTYPE__HPP
namespace matl 
{
   /**********************************************************************************************
    * class _MATRIX_
    * 
    * Base class for 'matl' matrix types.
    * 
    * Template Parameters:
    *     MAT       (Matrix) Matrix type (e.g. RowMatrix, ColMatrix, CSRMatrix, etc.)
    **********************************************************************************************/ 
    template<typename MAT>
    class _MATRIX_ : public _MATRIX_ALGEBRA_OBJECT_< _MATRIX_<MAT> >,
                     public _ELEMENTWISE_< _MATRIX_<MAT> >
    {
      public:  

        // Matrix datatype
        typedef typename type_traits<MAT>::dtype dtype;

        /**************** Assignment & Copying *****************/
        
        // Matrix copy
        template<typename MAT2>
        MAT& operator= (_MATRIX_<MAT2>& mat);

        // Expression evaluation trigger
        template<typename EXPR>
        MAT& operator=(_EXPRESSION_<EXPR> const& expr);

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

        // Scalar data writing (moved to derived classes. Not necessary here)
        //template<class T, enable_if_scalar<T> = true>
        //void operator= (T val);

        // Row Getter
        index_t rows() const;

        // Coloumn Getter
        index_t cols() const;

        // Size Getter
        index_t size() const;
    };
}
#endif
