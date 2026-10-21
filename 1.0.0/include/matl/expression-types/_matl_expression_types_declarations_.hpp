#ifndef MATL_EXPRESSION_TYPES_DECLARATIONS_HPP
#define MATL_EXPRESSION_TYPES_DECLARATIONS_HPP
namespace matl
{
    // Expression types. These are the types of expression objects, e.g. 'auto e = A + B;'.
    namespace expr
    {
        // Elementwise operations
        struct Add;
        struct Subtract;
        struct Multiply;
        struct Divide;

        // Elementwise operation on two operands (matrices, expressions, or one scalar)
        template<typename OP, typename LHS, typename RHS> class Elementwise;

        // Matrix product
        template<typename LHS, typename RHS> class Product;

        // Transpose
        template<typename MAO> class Transpose;

        // Expression that owns a (formerly temporary) matrix
        template<typename MAT> class Owned;
    }

    namespace detail
    {
        // Evaluators (see expression-types/evaluation/_matl_evaluation_.hpp)
        template<typename T> struct evaluator;
    }
}
#endif
