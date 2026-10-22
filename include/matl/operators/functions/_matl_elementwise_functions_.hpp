/**************************************************************************************************
 * File: _matl_elementwise_functions_.hpp
 *
 * Description: Elementwise functions of matrices and expressions:
 *
 *     -A                                      negation
 *     abs, sqrt, cbrt, exp, log, log2, log10  elementwise math functions (as in <cmath>)
 *     sin, cos, tan, asin, acos, atan,
 *     sinh, cosh, tanh, floor, ceil, round
 *     times(A, B)                             elementwise product    (MATLAB A .* B)
 *     rdivide(A, B)                           elementwise quotient   (MATLAB A ./ B)
 *     pow(A, B)                               elementwise power      (MATLAB A .^ B)
 *
 * In times, rdivide and pow either operand (but not both) may be a scalar, e.g. pow(A, 2) or
 * rdivide(1.0, A). Like the arithmetic operators, these build expressions and compute nothing
 * until the expression is assigned, and a temporary matrix operand is moved into the expression.
 *
 * Note: inside namespace matl, these names hide the <cmath> functions of the same names, so MATL's
 * own code calls the standard functions qualified (std::sqrt).
 **************************************************************************************************/
#ifndef MATL_ELEMENTWISE_FUNCTIONS_HPP
#define MATL_ELEMENTWISE_FUNCTIONS_HPP
namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * capture(x)
        *
        * How an expression built by a function stores an operand passed as a forwarding
        * reference: a temporary matrix is moved into an expr::Owned expression, anything else
        * is passed on as it is (and stored as usual, see stored_t).
        *******************************************************************************************/
        template<typename X>
        decltype(auto) capture(X&& x)
        {
            if constexpr (TemporaryMatrix<X>) return own_temporary(std::forward<X>(x));
            else                              return (x);
        }

        template<typename X>
        using captured_t = std::remove_cvref_t< decltype(capture(std::declval<X>())) >;

        // Unary elementwise expression of x
        template<typename OP, typename X>
        expr::Unary<OP, captured_t<X>> unary(X&& x)
        {return expr::Unary<OP, captured_t<X>>(capture(std::forward<X>(x)));}

        // Binary elementwise expression of lhs and rhs
        template<typename OP, typename L, typename R>
        expr::Elementwise<OP, captured_t<L>, captured_t<R>> binary(L&& lhs, R&& rhs)
        {
            return expr::Elementwise<OP, captured_t<L>, captured_t<R>>(capture(std::forward<L>(lhs)),
                                                                       capture(std::forward<R>(rhs)));
        }
    }

    // Elementwise operands: a matrix or expression and another operand that fits it (an operand of
    // the same compile-time dimensions, or a scalar)
    template<typename L, typename R>
    concept ElementwiseOperands = Operand<L> && Operand<R> &&
                                  (MatrixAlgebraObject<L> || MatrixAlgebraObject<R>) && SizesMatch<L, R>;


    /****************************************** Negation ******************************************/

    template<MatrixAlgebraObject X>
    auto operator- (X&& x) {return detail::unary<expr::Negate>(std::forward<X>(x));}


    /******************************* Elementwise math functions *********************************/

    template<MatrixAlgebraObject X> auto abs  (X&& x) {return detail::unary<expr::Abs>  (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto sqrt (X&& x) {return detail::unary<expr::Sqrt> (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto cbrt (X&& x) {return detail::unary<expr::Cbrt> (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto exp  (X&& x) {return detail::unary<expr::Exp>  (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto log  (X&& x) {return detail::unary<expr::Log>  (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto log2 (X&& x) {return detail::unary<expr::Log2> (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto log10(X&& x) {return detail::unary<expr::Log10>(std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto sin  (X&& x) {return detail::unary<expr::Sin>  (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto cos  (X&& x) {return detail::unary<expr::Cos>  (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto tan  (X&& x) {return detail::unary<expr::Tan>  (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto asin (X&& x) {return detail::unary<expr::Asin> (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto acos (X&& x) {return detail::unary<expr::Acos> (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto atan (X&& x) {return detail::unary<expr::Atan> (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto sinh (X&& x) {return detail::unary<expr::Sinh> (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto cosh (X&& x) {return detail::unary<expr::Cosh> (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto tanh (X&& x) {return detail::unary<expr::Tanh> (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto floor(X&& x) {return detail::unary<expr::Floor>(std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto ceil (X&& x) {return detail::unary<expr::Ceil> (std::forward<X>(x));}
    template<MatrixAlgebraObject X> auto round(X&& x) {return detail::unary<expr::Round>(std::forward<X>(x));}


    /************************************ Binary elementwise ************************************/

    // Elementwise product (MATLAB A .* B); A*B is the matrix product
    template<typename L, typename R> requires ElementwiseOperands<L, R>
    auto times(L&& lhs, R&& rhs)
    {return detail::binary<expr::Multiply>(std::forward<L>(lhs), std::forward<R>(rhs));}

    // Elementwise quotient (MATLAB A ./ B)
    template<typename L, typename R> requires ElementwiseOperands<L, R>
    auto rdivide(L&& lhs, R&& rhs)
    {return detail::binary<expr::Divide>(std::forward<L>(lhs), std::forward<R>(rhs));}

    // Elementwise power (MATLAB A .^ B), computed with std::pow
    template<typename L, typename R> requires ElementwiseOperands<L, R>
    auto pow(L&& lhs, R&& rhs)
    {return detail::binary<expr::Power>(std::forward<L>(lhs), std::forward<R>(rhs));}
}
#endif
