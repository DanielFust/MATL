/**************************************************************************************************
 * File: _matl_elementwise_prototype_.hpp
 *
 * Description: Prototype for the elementwise expression template, which represents every
 *              elementwise operation (+, -, * and / by a scalar) on matrices, expressions and
 *              scalars.
 **************************************************************************************************/
#ifndef MATL_ELEMENTWISE_PROTOTYPE_HPP
#define MATL_ELEMENTWISE_PROTOTYPE_HPP
namespace matl
{
    namespace expr
    {
       /*******************************************************************************************
        * Elementwise operations
        *
        * Each operation combines one element of each operand.
        *******************************************************************************************/
        struct Add
        {
            template<typename A, typename B>
            static constexpr auto apply(A a, B b) {return a + b;}
        };

        struct Subtract
        {
            template<typename A, typename B>
            static constexpr auto apply(A a, B b) {return a - b;}
        };

        struct Multiply
        {
            template<typename A, typename B>
            static constexpr auto apply(A a, B b) {return a * b;}
        };

        struct Divide
        {
            template<typename A, typename B>
            static constexpr auto apply(A a, B b) {return a / b;}
        };

        struct Power
        {
            template<typename A, typename B>
            static auto apply(A a, B b) {return std::pow(a, b);}
        };
    }

    namespace detail
    {
       /*******************************************************************************************
        * traits<expr::Elementwise>
        *
        * The result has the dimensions of the matrix operands (the left operand unless it is a
        * scalar). It is fixed-size if every matrix operand is fixed-size.
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        struct traits< expr::Elementwise<OP, LHS, RHS> >
        {
            // Operand that determines the dimensions
            typedef std::conditional_t<Scalar<LHS>, RHS, LHS> shape_type;

            typedef decltype(OP::apply(std::declval<dtype_t<LHS>>(), std::declval<dtype_t<RHS>>())) dtype;
            static constexpr bool         fixed_size    = traits<LHS>::fixed_size && traits<RHS>::fixed_size;
            static constexpr index_t      rows          = fixed_size ? traits<shape_type>::rows : 0;
            static constexpr index_t      cols          = fixed_size ? traits<shape_type>::cols : 0;
            static constexpr StorageOrder storage_order = combine(traits<LHS>::storage_order,
                                                                  traits<RHS>::storage_order);
        };

        // Single-index access (see linear_index_v) if both operands have it
        template<typename OP, typename LHS, typename RHS>
        inline constexpr bool linear_index_v< expr::Elementwise<OP, LHS, RHS> > = linear_index_v<LHS> && linear_index_v<RHS>;
    }

    namespace expr
    {
       /*******************************************************************************************
        * class Elementwise<OP, LHS, RHS>
        *
        * Expression template for an elementwise operation: element (i,j) of the result is
        * OP::apply(lhs(i,j), rhs(i,j)). Either operand (but not both) may be a scalar, which
        * stands for a matrix with every element equal to the scalar. For example:
        *
        *     A + B       Elementwise<Add, A, B>
        *     A - 2.0     Elementwise<Subtract, A, double>
        *     2.0 - A     Elementwise<Subtract, double, A>
        *     A * 2.0     Elementwise<Multiply, A, double>
        *     A / 2       Elementwise<Divide, A, int>  (division by a floating-point scalar is
        *                                               multiplication by its reciprocal)
        *
        * The expression is a formula: it stores its operands and computes nothing until it is
        * assigned to a matrix, which evaluates it through detail::evaluator<Elementwise>.
        *
        * Template Parameters:
        *     OP         Operation (Add, Subtract, Multiply, Divide or Power)
        *     LHS        Left operand type (matrix, expression or scalar)
        *     RHS        Right operand type (matrix, expression or scalar)
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        class Elementwise : public ExpressionBase< Elementwise<OP, LHS, RHS> >
        {
            static_assert(!(Scalar<LHS> && Scalar<RHS>), "Elementwise: at least one operand must be a matrix or expression");

          public:
            // Resultant matrix data type (e.g. double, int, float, etc.)
            typedef typename detail::traits<Elementwise>::dtype dtype;

            // Operation and operand types
            typedef OP  operation;
            typedef LHS lhs_type;
            typedef RHS rhs_type;

          private:
            // Operands (matrices by reference, expressions and scalars by value)
            detail::stored_t<LHS> lhs_;
            detail::stored_t<RHS> rhs_;

          public:

            // Constructor
            Elementwise(LHS const&, RHS const&);

            // Left operand
            LHS const& lhs() const;

            // Right operand
            RHS const& rhs() const;

            // Get the number of rows
            index_t rows() const;

            // Get the number of columns
            index_t cols() const;

            // Get the number of elements
            index_t size() const;

            // Check that the operands have the same dimensions (recursively)
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
