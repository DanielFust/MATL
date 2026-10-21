/**************************************************************************************************
 * File: _matl_scalar_minus_matrix_prototype_.hpp
 * 
 * Description: Prototype for the scalar minus matrix expression template.
 **************************************************************************************************/

namespace matl
{
   /***********************************************************************************************
    * class _SCALAR_MINUS_MATRIX_
    * 
    * Generic scalar-matrix subtraction expression template
    ***********************************************************************************************/
    template<class T, typename MAO>
    class _SCALAR_MINUS_MATRIX_ : public _EXPRESSION_<_SCALAR_MINUS_MATRIX_<T, MAO> >
    {
      public:
        // Typedefs
        typedef typename promote<T, typename type_traits<MAO>::dtype>::type dtype;

      private:
        // Expression or matrix
        typename member_type<MAO>::elementwise_expression mat_;

        // Scalar
        dtype scalar_;  
      
      public:  

        // Constructor (elementwise object)
        template<typename MAO_               = MAO,
                 enable_if_elementwise<MAO_> = true,
                 enable_if_scalar<T>         = true>
        _SCALAR_MINUS_MATRIX_(T, MAO_ const&);

        // Constructor (non-elementwise object)
        template<typename MAO_                = MAO,
                 disable_if_elementwise<MAO_> = true,
                 enable_if_scalar<T>          = true>
        _SCALAR_MINUS_MATRIX_(T, MAO_ const&);

        // Check if evaluation target is safe
        template<typename MAT>
        constexpr bool unsafe_evaluation_target(_MATRIX_<MAT> const&) const;

        // Check if expression references address of target matrix
        template<typename MAT>
        constexpr bool references_address(_MATRIX_<MAT> const*) const;

        // Evaluate the expression to the target matrix
        template<typename MAT_T>
        void eval(_MATRIX_<MAT_T>&) const;

        // Get the number of rows
        inline index_t rows() const;

        // Get the number of columns
        inline index_t cols() const;

        // Get the number of elements
        inline index_t size() const;

        // Get the element at the i-th row and j-th column
        dtype operator() (index_t i, index_t j) const;

        // Get the element at the i-th index
        dtype operator() (index_t i) const;

    }; 
}