/**************************************************************************************************
 * File: _matl_arithmetic_operators_.hpp
 *
 * Description: Arithmetic operators (+, -, *, /) for constructing expression objects.
 *
 * Each operator builds an expression (a formula) from its operands and computes nothing; the
 * expression is evaluated when it is assigned to a matrix. The operand types are deduced
 * exactly, so the expression type records them, e.g.
 *
 *     A + B * 2.0    expr::Elementwise<expr::Add, A, expr::Elementwise<expr::Multiply, B, double>>
 *
 * The operators are constrained to MATL matrices, expressions (MatrixAlgebraObject) and
 * arithmetic scalars, so they never apply to unrelated types. Operands whose dimensions are known
 * at compile time must fit (SizesMatch, InnerSizesMatch), e.g. adding a 2x3 and a 3x2 fixed-size
 * matrix does not compile.
 **************************************************************************************************/
#ifndef MATL_ARITHMETIC_OPERATORS_HPP
#define MATL_ARITHMETIC_OPERATORS_HPP
namespace matl
{
    /******************************************* Addition *****************************************/

    // Matrix or expression + matrix or expression
    template<MatrixAlgebraObject LHS, MatrixAlgebraObject RHS> requires SizesMatch<LHS, RHS>
    expr::Elementwise<expr::Add, LHS, RHS> operator+ (LHS const& lhs, RHS const& rhs)
    {return expr::Elementwise<expr::Add, LHS, RHS>(lhs, rhs);}

    // Matrix or expression + scalar
    template<MatrixAlgebraObject LHS, Scalar S>
    expr::Elementwise<expr::Add, LHS, S> operator+ (LHS const& lhs, S rhs)
    {return expr::Elementwise<expr::Add, LHS, S>(lhs, rhs);}

    // Scalar + matrix or expression
    template<Scalar S, MatrixAlgebraObject RHS>
    expr::Elementwise<expr::Add, S, RHS> operator+ (S lhs, RHS const& rhs)
    {return expr::Elementwise<expr::Add, S, RHS>(lhs, rhs);}


    /***************************************** Subtraction ****************************************/

    // Matrix or expression - matrix or expression
    template<MatrixAlgebraObject LHS, MatrixAlgebraObject RHS> requires SizesMatch<LHS, RHS>
    expr::Elementwise<expr::Subtract, LHS, RHS> operator- (LHS const& lhs, RHS const& rhs)
    {return expr::Elementwise<expr::Subtract, LHS, RHS>(lhs, rhs);}

    // Matrix or expression - scalar
    template<MatrixAlgebraObject LHS, Scalar S>
    expr::Elementwise<expr::Subtract, LHS, S> operator- (LHS const& lhs, S rhs)
    {return expr::Elementwise<expr::Subtract, LHS, S>(lhs, rhs);}

    // Scalar - matrix or expression
    template<Scalar S, MatrixAlgebraObject RHS>
    expr::Elementwise<expr::Subtract, S, RHS> operator- (S lhs, RHS const& rhs)
    {return expr::Elementwise<expr::Subtract, S, RHS>(lhs, rhs);}


    /*************************************** Multiplication ***************************************/

    // Matrix or expression * matrix or expression (matrix product)
    template<MatrixAlgebraObject LHS, MatrixAlgebraObject RHS> requires InnerSizesMatch<LHS, RHS>
    expr::Product<LHS, RHS> operator* (LHS const& lhs, RHS const& rhs)
    {return expr::Product<LHS, RHS>(lhs, rhs);}

    // Matrix or expression * scalar
    template<MatrixAlgebraObject LHS, Scalar S>
    expr::Elementwise<expr::Multiply, LHS, S> operator* (LHS const& lhs, S rhs)
    {return expr::Elementwise<expr::Multiply, LHS, S>(lhs, rhs);}

    // Scalar * matrix or expression
    template<Scalar S, MatrixAlgebraObject RHS>
    expr::Elementwise<expr::Multiply, S, RHS> operator* (S lhs, RHS const& rhs)
    {return expr::Elementwise<expr::Multiply, S, RHS>(lhs, rhs);}


    /****************************************** Division ******************************************/

    // Matrix or expression / scalar. Division by a floating-point scalar is implemented as
    // multiplication by its reciprocal (one division instead of one per element). Integer
    // division cannot be replaced by a reciprocal, so it divides every element.
    template<MatrixAlgebraObject LHS, Scalar S>
    auto operator/ (LHS const& lhs, S rhs)
    {
        if constexpr (std::is_integral_v<S>)
            return expr::Elementwise<expr::Divide, LHS, S>(lhs, rhs);
        else
        {
            typedef detail::promote_t<detail::dtype_t<LHS>, S> dtype;
            return expr::Elementwise<expr::Multiply, LHS, dtype>(lhs, dtype(1) / rhs);
        }
    }
}
#endif
