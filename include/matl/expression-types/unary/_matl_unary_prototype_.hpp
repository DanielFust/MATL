/**************************************************************************************************
 * File: _matl_unary_prototype_.hpp
 *
 * Description: Prototype for the unary elementwise expression template, which applies a function
 *              to every element of a matrix or expression (-A, abs(A), sqrt(A), exp(A), ...).
 **************************************************************************************************/
#ifndef MATL_UNARY_PROTOTYPE_HPP
#define MATL_UNARY_PROTOTYPE_HPP
namespace matl
{
    namespace expr
    {
       /*******************************************************************************************
        * Unary elementwise operations
        *
        * Each operation maps one element to one element. The standard library functions are
        * called qualified (std::), so they are never confused with MATL's functions of the same
        * names, which take matrices.
        *******************************************************************************************/
        struct Negate {template<typename A> static constexpr auto apply(A a) {return -a;}};
        struct Abs    {template<typename A> static auto apply(A a) {return std::abs(a);}};
        struct Sqrt   {template<typename A> static auto apply(A a) {return std::sqrt(a);}};
        struct Cbrt   {template<typename A> static auto apply(A a) {return std::cbrt(a);}};
        struct Exp    {template<typename A> static auto apply(A a) {return std::exp(a);}};
        struct Log    {template<typename A> static auto apply(A a) {return std::log(a);}};
        struct Log2   {template<typename A> static auto apply(A a) {return std::log2(a);}};
        struct Log10  {template<typename A> static auto apply(A a) {return std::log10(a);}};
        struct Sin    {template<typename A> static auto apply(A a) {return std::sin(a);}};
        struct Cos    {template<typename A> static auto apply(A a) {return std::cos(a);}};
        struct Tan    {template<typename A> static auto apply(A a) {return std::tan(a);}};
        struct Asin   {template<typename A> static auto apply(A a) {return std::asin(a);}};
        struct Acos   {template<typename A> static auto apply(A a) {return std::acos(a);}};
        struct Atan   {template<typename A> static auto apply(A a) {return std::atan(a);}};
        struct Sinh   {template<typename A> static auto apply(A a) {return std::sinh(a);}};
        struct Cosh   {template<typename A> static auto apply(A a) {return std::cosh(a);}};
        struct Tanh   {template<typename A> static auto apply(A a) {return std::tanh(a);}};
        struct Floor  {template<typename A> static auto apply(A a) {return std::floor(a);}};
        struct Ceil   {template<typename A> static auto apply(A a) {return std::ceil(a);}};
        struct Round  {template<typename A> static auto apply(A a) {return std::round(a);}};
    }

    namespace detail
    {
       /*******************************************************************************************
        * traits<expr::Unary>
        *
        * The result has the operand's dimensions and storage order. Its element type is the type
        * the operation returns (e.g. sqrt of an int matrix gives double elements).
        *******************************************************************************************/
        template<typename OP, typename MAO>
        struct traits< expr::Unary<OP, MAO> >
        {
            typedef decltype(OP::apply(std::declval<dtype_t<MAO>>())) dtype;
            static constexpr bool         fixed_size    = traits<MAO>::fixed_size;
            static constexpr index_t      rows          = traits<MAO>::rows;
            static constexpr index_t      cols          = traits<MAO>::cols;
            static constexpr StorageOrder storage_order = traits<MAO>::storage_order;
        };

        // Single-index access (see linear_index_v) if the operand has it
        template<typename OP, typename MAO>
        inline constexpr bool linear_index_v< expr::Unary<OP, MAO> > = linear_index_v<MAO>;
    }

    namespace expr
    {
       /*******************************************************************************************
        * class Unary<OP, MAO>
        *
        * Expression template for a unary elementwise operation: element (i,j) of the result is
        * OP::apply(operand(i,j)). For example -A is Unary<Negate, A>, and sqrt(A) is
        * Unary<Sqrt, A>.
        *
        * The expression is a formula: it stores its operand and computes nothing until it is
        * assigned to a matrix, which evaluates it through detail::evaluator<Unary>.
        *
        * Template Parameters:
        *     OP         Operation (Negate, Abs, Sqrt, ...)
        *     MAO        Operand type (matrix or expression)
        *******************************************************************************************/
        template<typename OP, typename MAO>
        class Unary : public ExpressionBase< Unary<OP, MAO> >
        {
          public:
            // Resultant matrix data type (e.g. double, int, float, etc.)
            typedef typename detail::traits<Unary>::dtype dtype;

            // Operation and operand types
            typedef OP  operation;
            typedef MAO operand_type;

          private:
            // Operand (matrix by reference, expression by value)
            detail::stored_t<MAO> mat_;

          public:

            // Constructor
            explicit Unary(MAO const&);

            // Operand
            MAO const& operand() const;

            // Get the number of rows
            index_t rows() const;

            // Get the number of columns
            index_t cols() const;

            // Get the number of elements
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
