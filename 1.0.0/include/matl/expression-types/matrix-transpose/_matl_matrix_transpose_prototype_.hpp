#ifndef _MATL_MATRIX_TRANSPOSE_PROTOTYPE__HPP
#define _MATL_MATRIX_TRANSPOSE_PROTOTYPE__HPP

namespace matl
{
   /***********************************************************************************************
    * class _MATRIX_TRANSPOSE_
    *
    * Expression template for the transpose of a matrix or expression.
    *
    * The expression is a formula: it stores its operand and computes nothing until it is assigned
    * to a matrix, which evaluates it through evaluator<_MATRIX_TRANSPOSE_>.
    *
    * Template Parameters:
    *     MAO        Operand type (matrix or expression)
    *************************************************************************************************/
    template<typename MAO>
    class _MATRIX_TRANSPOSE_ : public _EXPRESSION_<_MATRIX_TRANSPOSE_<MAO> >
    {
      public:
        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename type_traits<MAO>::dtype dtype;

        // Operand type as stored (without references)
        typedef std::remove_cvref_t<typename member_type<MAO>::elementwise_expression> operand_type;

      private:
        // Operand (expression by value, matrix by reference)
        typename member_type<MAO>::elementwise_expression mat_;

      public:

        // Constructor
        explicit _MATRIX_TRANSPOSE_(MAO const&);

        // Operand being transposed
        operand_type const& operand() const;

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

        // Evaluate the expression to the target matrix
        template<typename MAT>
        void eval(_MATRIX_<MAT>&) const;
    };
}
#endif
