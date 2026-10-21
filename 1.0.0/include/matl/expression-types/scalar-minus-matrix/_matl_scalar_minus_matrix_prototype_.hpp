/**************************************************************************************************
 * File: _matl_scalar_minus_matrix_prototype_.hpp
 * 
 * Description: Prototype for the scalar minus matrix expression template.
 **************************************************************************************************/
#ifndef _MATL_SCALAR_MINUS_MATRIX_PROTOTYPE__HPP
#define _MATL_SCALAR_MINUS_MATRIX_PROTOTYPE__HPP

namespace matl 
{
   /***********************************************************************************************
    * class _SCALAR_MINUS_MATRIX_
    * 
    * Expression template for the difference of a scalar and a matrix or expression (every element subtracted from the scalar).
    * 
    * The expression is a formula: it stores its operands and computes nothing until it is
    * assigned to a matrix, which evaluates it through evaluator<_SCALAR_MINUS_MATRIX_>.
    * 
    * Template Parameters:
    *     MAO        Operand type (matrix or expression)
    *     T          Scalar type
    ***********************************************************************************************/
    template<class T, typename MAO>
    class _SCALAR_MINUS_MATRIX_ : public _EXPRESSION_<_SCALAR_MINUS_MATRIX_<T,MAO> >
    {
      public:
        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename promote<T, typename type_traits<MAO>::dtype>::type dtype;

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
        _SCALAR_MINUS_MATRIX_(T, OP const&);

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
