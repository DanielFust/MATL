/**************************************************************************************************
 * File: _matl_matrix_plus_scalar_prototype_.hpp
 * 
 * Description: Prototypes for the matrix plus scalar expression.
 **************************************************************************************************/
namespace matl
{
   /***********************************************************************************************
    * class _MATRIX_PLUS_SCALAR_ 
    * 
    * Expression template class for the matrix plus scalar operation.
    * Template parameters:
    *     MAO: The expression or matrix type.
    *     T: The scalar type.
    ***********************************************************************************************/
    template<typename MAO, class T>
    class _MATRIX_PLUS_SCALAR_ : public _EXPRESSION_<_MATRIX_PLUS_SCALAR_<MAO, T> >
    {
      public:
        // Typedefs
        typedef typename promote<typename type_traits<MAO>::dtype, T>::type dtype;

      private:
        // Expression or matrix
        typename member_type<MAO>::elementwise_expression mat_;
        
        // Scalar
        dtype scalar_;  
        
      public:

        // Constructor (elementwise + scalar)
        template<typename MAO_               = MAO, 
                 class    T_                 = T,
                 enable_if_elementwise<MAO_> = true, 
                 enable_if_scalar<T_>        = true>
        _MATRIX_PLUS_SCALAR_(MAO_ const&, T_);

        // Constructor (non-elementwise + scalar)
        template<typename MAO_                = MAO, 
                 class    T_                  = T,
                 disable_if_elementwise<MAO_> = true, 
                 enable_if_scalar<T_>         = true>
        _MATRIX_PLUS_SCALAR_(MAO_ const&, T_);

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
        dtype operator() (size_t i, size_t j) const;

        // Get the element at the i-th index
        dtype operator() (size_t i) const;

    };
}