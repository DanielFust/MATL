/**************************************************************************************************
 * File: _matl_transpose_prototype_.hpp
 *
 * Description: Prototype for the transpose expression template.
 **************************************************************************************************/
#ifndef MATL_TRANSPOSE_PROTOTYPE_HPP
#define MATL_TRANSPOSE_PROTOTYPE_HPP

namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * traits<expr::Transpose>
        *
        * The transpose of an M x N operand is N x M, with the opposite storage order.
        *******************************************************************************************/
        template<typename MAO>
        struct traits< expr::Transpose<MAO> >
        {
            typedef dtype_t<MAO> dtype;
            static constexpr bool         fixed_size    = traits<MAO>::fixed_size;
            static constexpr index_t      rows          = traits<MAO>::cols;
            static constexpr index_t      cols          = traits<MAO>::rows;
            static constexpr StorageOrder storage_order = transposed(traits<MAO>::storage_order);
        };

        template<typename MAO>
        inline constexpr bool is_transpose_v< expr::Transpose<MAO> > = true;

        // Single-index access (see linear_index_v) if the operand has it
        template<typename MAO>
        inline constexpr bool linear_index_v< expr::Transpose<MAO> > = linear_index_v<MAO>;
    }

    namespace expr
    {
       /*******************************************************************************************
        * class Transpose<MAO>
        *
        * Expression template for the transpose of a matrix or expression.
        *
        * The expression is a formula: it stores its operand and computes nothing until it is
        * assigned to a matrix, which evaluates it through detail::evaluator<Transpose>.
        *
        * Template Parameters:
        *     MAO        Operand type (matrix or expression)
        *******************************************************************************************/
        template<typename MAO>
        class Transpose : public ExpressionBase< Transpose<MAO> >
        {
          public:
            // Resultant matrix data type (e.g. double, int, float, etc.)
            typedef typename detail::traits<Transpose>::dtype dtype;

            // Operand type
            typedef MAO operand_type;

          private:
            // Operand (expression by value, matrix by reference)
            detail::stored_t<MAO> mat_;

          public:

            // Constructor
            explicit Transpose(MAO const&);

            // Operand being transposed
            MAO const& operand() const;

            // Get number of rows of resultant matrix
            index_t rows() const;

            // Get number of columns of resultant matrix
            index_t cols() const;

            // Get number of elements in resultant matrix
            index_t size() const;

            // Check the dimensions of the operand (recursively)
            void check_dimensions() const;

            // Check if evaluation target is safe
            template<WritableMatrix MAT>
            bool unsafe_evaluation_target(MAT const&) const;

            // Check if the expression reads memory that may hold elements of a target (see detail::Footprint)
            bool reads(detail::Footprint const&) const;

            // Evaluate the expression to the target matrix
            template<WritableMatrix MAT>
            void eval(MAT&) const;
        };
    }
}
#endif
