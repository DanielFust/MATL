/**************************************************************************************************
 * File: _matl_matrix_times_matrix_prototype_.hpp
 *
 * Description: Prototype for the matrix times matrix expression template.
 **************************************************************************************************/
#ifndef _MATL_MATRIX_TIMES_MATRIX_PROTOTYPE__HPP
#define _MATL_MATRIX_TIMES_MATRIX_PROTOTYPE__HPP
namespace matl
{
   /***********************************************************************************************
    * class _MATRIX_TIMES_MATRIX_
    *
    * Expression template for the product of two matrices or expressions.
    *
    * The expression is a formula: it stores its operands and computes nothing until it is
    * assigned to a matrix. When assigned directly (C = A*B), the product is computed into the
    * target. When nested inside a larger expression (C = A*B + D), evaluator<_MATRIX_TIMES_MATRIX_>
    * computes the product into a temporary matrix at the start of the assignment.
    *
    * Template Parameters:
    *     MAO1       Left operand type (matrix or expression)
    *     MAO2       Right operand type (matrix or expression)
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    class _MATRIX_TIMES_MATRIX_ : public _EXPRESSION_<_MATRIX_TIMES_MATRIX_<MAO1,MAO2> >
    {
      public:
        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename promote<typename type_traits<MAO1>::dtype,
                                 typename type_traits<MAO2>::dtype>::type dtype;

        // Operand types as stored (without references)
        typedef std::remove_cvref_t<typename member_type<MAO1>::elementwise_expression> lhs_type;
        typedef std::remove_cvref_t<typename member_type<MAO2>::elementwise_expression> rhs_type;

      private:
        // Operands (expressions by value, matrices by reference)
        typename member_type<MAO1>::elementwise_expression mat1_;
        typename member_type<MAO2>::elementwise_expression mat2_;

      public:

        // Constructor
        template<typename OP1, typename OP2>
        _MATRIX_TIMES_MATRIX_(OP1 const&, OP2 const&);

        // Left operand
        lhs_type const& lhs() const;

        // Right operand
        rhs_type const& rhs() const;

        // Get number of rows of resultant matrix
        index_t rows() const;

        // Get number of columns of resultant matrix
        index_t cols() const;

        // Get number of elements in resultant matrix
        index_t size() const;

        // Check if evaluation target is safe
        template<typename MAT>
        constexpr bool unsafe_evaluation_target(_MATRIX_<MAT> const&) const;

        // Check if expression references address of target matrix
        template<typename MAT>
        constexpr bool references_address(_MATRIX_<MAT> const*) const;

        // Evaluate matrix multiplication to target matrix
        template<typename MAT>
        void eval(_MATRIX_<MAT>& ) const;
    };
}
#endif
