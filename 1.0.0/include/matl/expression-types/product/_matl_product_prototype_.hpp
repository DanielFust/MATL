/**************************************************************************************************
 * File: _matl_product_prototype_.hpp
 *
 * Description: Prototype for the matrix product expression template.
 **************************************************************************************************/
#ifndef MATL_PRODUCT_PROTOTYPE_HPP
#define MATL_PRODUCT_PROTOTYPE_HPP
namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * traits<expr::Product>
        *
        * The product of an M x K and a K x N operand is M x N. Products are computed into a
        * default (row-major) matrix.
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        struct traits< expr::Product<LHS, RHS> >
        {
            typedef promote_t< dtype_t<LHS>, dtype_t<RHS> > dtype;
            static constexpr bool         fixed_size    = traits<LHS>::fixed_size && traits<RHS>::fixed_size;
            static constexpr index_t      rows          = fixed_size ? traits<LHS>::rows : 0;
            static constexpr index_t      cols          = fixed_size ? traits<RHS>::cols : 0;
            static constexpr StorageOrder storage_order = StorageOrder::RowMajor;
        };

        template<typename LHS, typename RHS>
        inline constexpr bool is_product_v< expr::Product<LHS, RHS> > = true;
    }

    namespace expr
    {
       /*******************************************************************************************
        * class Product<LHS, RHS>
        *
        * Expression template for the product of two matrices or expressions.
        *
        * The expression is a formula: it stores its operands and computes nothing until it is
        * assigned to a matrix. When assigned directly (C = A*B), the product is computed into the
        * target. When nested inside a larger expression (C = A*B + D), detail::evaluator<Product>
        * computes the product into a temporary matrix at the start of the assignment.
        *
        * Template Parameters:
        *     LHS        Left operand type (matrix or expression)
        *     RHS        Right operand type (matrix or expression)
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        class Product : public ExpressionBase< Product<LHS, RHS> >
        {
          public:
            // Resultant matrix data type (e.g. double, int, float, etc.)
            typedef typename detail::traits<Product>::dtype dtype;

            // Operand types
            typedef LHS lhs_type;
            typedef RHS rhs_type;

          private:
            // Operands (matrices by reference, expressions by value)
            detail::stored_t<LHS> lhs_;
            detail::stored_t<RHS> rhs_;

          public:

            // Constructor
            Product(LHS const&, RHS const&);

            // Left operand
            LHS const& lhs() const;

            // Right operand
            RHS const& rhs() const;

            // Get number of rows of resultant matrix
            index_t rows() const;

            // Get number of columns of resultant matrix
            index_t cols() const;

            // Get number of elements in resultant matrix
            index_t size() const;

            // Check if evaluation target is safe
            template<MatrixType MAT>
            bool unsafe_evaluation_target(MAT const&) const;

            // Check if expression references the matrix at the given address
            bool references_address(void const*) const;

            // Evaluate matrix multiplication to target matrix
            template<MatrixType MAT>
            void eval(MAT&) const;
        };
    }
}
#endif
