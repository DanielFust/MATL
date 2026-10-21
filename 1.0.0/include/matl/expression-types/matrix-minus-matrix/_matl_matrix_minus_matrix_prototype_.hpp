/**************************************************************************************************
 * File: _matl_matrix_minus_matrix_prototype_.hpp
 * 
 * Prototype for the matrix minus matrix expression template.
 **************************************************************************************************/
#ifndef _MATL_MATRIX_MINUS_MATRIX_PROTOTYPE__HPP
#define _MATL_MATRIX_MINUS_MATRIX_PROTOTYPE__HPP

namespace matl 
{
   /***********************************************************************************************
    * class _MATRIX_MINUS_MATRIX_
    * 
    * Expression template for the elementwise difference of two matrices or expressions.
    * 
    * The expression is a formula: it stores its operands and computes nothing until it is
    * assigned to a matrix, which evaluates it through evaluator<_MATRIX_MINUS_MATRIX_>.
    * 
    * Template Parameters:
    *     MAO1       Left operand type (matrix or expression)
    *     MAO2       Right operand type (matrix or expression)
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    class _MATRIX_MINUS_MATRIX_ : public _EXPRESSION_<_MATRIX_MINUS_MATRIX_<MAO1,MAO2> >
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
        _MATRIX_MINUS_MATRIX_(OP1 const&, OP2 const&);

        // Left operand
        lhs_type const& lhs() const;

        // Right operand
        rhs_type const& rhs() const;

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
