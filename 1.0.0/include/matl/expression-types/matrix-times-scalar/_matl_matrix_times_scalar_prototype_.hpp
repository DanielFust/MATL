#ifndef _MATL_MATRIX_TIMES_SCALAR_PROTOTYPE__HPP
#define _MATL_MATRIX_TIMES_SCALAR_PROTOTYPE__HPP

namespace matl 
{
   /***********************************************************************************************
    * class _MATRIX_TIMES_SCALAR_
    * 
    * Expression template for the product of a matrix or expression and a scalar.
    * 
    * The expression is a formula: it stores its operands and computes nothing until it is
    * assigned to a matrix, which evaluates it through evaluator<_MATRIX_TIMES_SCALAR_>.
    * 
    * Template Parameters:
    *     MAO        Operand type (matrix or expression)
    *     T          Scalar type
    ***********************************************************************************************/
    template<typename MAO, class T>
    class _MATRIX_TIMES_SCALAR_ : public _EXPRESSION_<_MATRIX_TIMES_SCALAR_<MAO,T> >
    {
      public:
        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename promote<typename type_traits<MAO>::dtype, T>::type dtype;

        // Operand type as stored (without references)
        typedef std::remove_cvref_t<typename member_type<MAO>::elementwise_expression> operand_type;

      private:
        // Operand (expression by value, matrix by reference)
        typename member_type<MAO>::elementwise_expression mat_;

        // Scalar
        dtype scalar_;

      public:

        // Constructor
        template<typename OP>
        _MATRIX_TIMES_SCALAR_(OP const&, T);

        // Matrix or expression operand
        operand_type const& operand() const;

        // Scalar operand
        dtype scalar() const;

        // Get the number of rows
        index_t rows() const;

        // Get the number of columns
        index_t cols() const;

        // Get the number of elements
        index_t size() const;

        // Check if evaluation target is safe
        template<typename MAT>
        constexpr bool unsafe_evaluation_target(_MATRIX_<MAT> const&) const;

        // Check if expression references address of target matrix
        template<typename MAT>
        constexpr bool references_address(_MATRIX_<MAT> const*) const;

        // Evaluate the expression to the target matrix
        template<typename MAT>
        void eval(_MATRIX_<MAT>&) const;
    };
}

#endif
