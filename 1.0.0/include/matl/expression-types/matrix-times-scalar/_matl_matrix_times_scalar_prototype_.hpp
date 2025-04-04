#ifndef _MATL_MATRIX_TIMES_SCALAR_PROTOTYPE__HPP
#define _MATL_MATRIX_TIMES_SCALAR_PROTOTYPE__HPP

namespace matl 
{
   /***********************************************************************************************
    * class _MATRIX_TIMES_SCALAR_ 
    * 
    * Expression template class for the matrix times scalar operation.
    ***********************************************************************************************/
    template<typename MAO, class T> 
    class _MATRIX_TIMES_SCALAR_ : public _EXPRESSION_<_MATRIX_TIMES_SCALAR_<MAO, T> >
    {
      public:
        // Typedefs
        typedef typename promote<typename type_traits<MAO>::dtype, T>::type dtype;

      private:
        // Expression or matrix
        member_type<MAO>::elementwise_expression mat_;

        // Scalar
        dtype scalar_; 

      public:

        // Constructor (elementwise * scalar)
        template<typename MAO_               = MAO, 
                 enable_if_elementwise<MAO_> = true, 
                 enable_if_scalar<T>         = true>
        _MATRIX_TIMES_SCALAR_(MAO_ const&, T);

        // Constructor (non-elementwise * scalar)
        template<typename MAO_                = MAO, 
                 disable_if_elementwise<MAO_> = true, 
                 enable_if_scalar<T>          = true>
        _MATRIX_TIMES_SCALAR_(MAO_ const&, T);

        // Check if evaluation target is safe
        template<typename MAT>
        constexpr bool unsafe_evaluation_target(_MATRIX_<MAT> const&) const;

        // Check if expression references address of target matrix
        template<typename MAT>
        constexpr bool references_address(_MATRIX_<MAT> const*) const;

        // Evaluate the expression to the target matrix
        template<typename MAT>
        void eval(_MATRIX_<MAT>&) const;

        // Get the number of rows
        inline index_t rows() const;

        // Get the number of columns
        inline index_t cols() const;

        // Get the number of elements
        inline index_t size() const;

        // Get the element at the i-th row and j-th column
        dtype operator() (index_t, index_t) const;

        // Get the element at the i-th index
        dtype operator() (index_t) const;
 
    };
}

#endif