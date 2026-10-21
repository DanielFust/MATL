/**************************************************************************************************
 * File: _matl_temporary_operands_.hpp
 *
 * Description: Operator overloads for temporary (rvalue) matrix operands.
 *
 * Expressions hold named matrices by reference, so a temporary matrix operand, e.g.
 *
 *     auto e = A + eye<double>(3);
 *
 * would be destroyed at the end of the statement while the expression still refers to it. These
 * overloads move temporary matrices into an expr::Owned expression, which the resulting
 * expression stores by value, and then forward to the regular operators.
 *
 * Overload resolution prefers these overloads for rvalue matrix arguments (binding an rvalue to
 * 'M&&' is a better match than binding it to the 'const&' parameters of the regular operators).
 * For lvalue matrix arguments the TemporaryMatrix constraint fails, so the regular operators are
 * used.
 **************************************************************************************************/
#ifndef MATL_TEMPORARY_OPERANDS_HPP
#define MATL_TEMPORARY_OPERANDS_HPP
namespace matl
{
   /***********************************************************************************************
    * expr::Owned<M> own_temporary(M&& mat)
    *
    * Moves a temporary matrix into an expr::Owned expression.
    ***********************************************************************************************/
    template<TemporaryMatrix M>
    expr::Owned<std::remove_cv_t<M> > own_temporary(M&& mat)
    {return expr::Owned<std::remove_cv_t<M> >(std::forward<M>(mat));}


    /******************************************* Addition *****************************************/

    // Temporary matrix + operand
    template<TemporaryMatrix M, Operand R>
    auto operator+ (M&& mat, R&& rhs)
    {return own_temporary(std::forward<M>(mat)) + std::forward<R>(rhs);}

    // Operand + temporary matrix
    template<Operand L, TemporaryMatrix M> requires (!TemporaryMatrix<L>)
    auto operator+ (L&& lhs, M&& mat)
    {return std::forward<L>(lhs) + own_temporary(std::forward<M>(mat));}


    /***************************************** Subtraction ****************************************/

    // Temporary matrix - operand
    template<TemporaryMatrix M, Operand R>
    auto operator- (M&& mat, R&& rhs)
    {return own_temporary(std::forward<M>(mat)) - std::forward<R>(rhs);}

    // Operand - temporary matrix
    template<Operand L, TemporaryMatrix M> requires (!TemporaryMatrix<L>)
    auto operator- (L&& lhs, M&& mat)
    {return std::forward<L>(lhs) - own_temporary(std::forward<M>(mat));}


    /*************************************** Multiplication ***************************************/

    // Temporary matrix * operand
    template<TemporaryMatrix M, Operand R>
    auto operator* (M&& mat, R&& rhs)
    {return own_temporary(std::forward<M>(mat)) * std::forward<R>(rhs);}

    // Operand * temporary matrix
    template<Operand L, TemporaryMatrix M> requires (!TemporaryMatrix<L>)
    auto operator* (L&& lhs, M&& mat)
    {return std::forward<L>(lhs) * own_temporary(std::forward<M>(mat));}


    /****************************************** Division ******************************************/

    // Temporary matrix / scalar
    template<TemporaryMatrix M, Scalar T>
    auto operator/ (M&& mat, T scalar)
    {return own_temporary(std::forward<M>(mat)) / scalar;}


    /****************************************** Transpose *****************************************/

    // transpose(temporary matrix)
    template<TemporaryMatrix M>
    auto transpose(M&& mat)
    {return transpose(own_temporary(std::forward<M>(mat)));}
}
#endif
