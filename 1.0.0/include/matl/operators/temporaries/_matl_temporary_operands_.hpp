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
 * overloads move temporary matrices into an _OWNED_MATRIX_ expression, which the resulting
 * expression stores by value, and then forward to the regular operators.
 *
 * Overload resolution prefers these overloads for rvalue matrix arguments (binding an rvalue to
 * 'M&&' is an exact match, whereas the regular operators bind '_MATRIX_<MAT> const&' through a
 * derived-to-base conversion). For lvalue matrix arguments the TemporaryMatrix constraint fails,
 * so the regular operators are used.
 **************************************************************************************************/
#ifndef _MATL_TEMPORARY_OPERANDS__HPP
#define _MATL_TEMPORARY_OPERANDS__HPP
namespace matl
{
   /***********************************************************************************************
    * _OWNED_MATRIX_ own_temporary(M&& mat)
    *
    * Moves a temporary matrix into an _OWNED_MATRIX_ expression.
    ***********************************************************************************************/
    template<TemporaryMatrix M>
    _OWNED_MATRIX_<std::remove_cv_t<M> > own_temporary(M&& mat)
    {return _OWNED_MATRIX_<std::remove_cv_t<M> >(std::forward<M>(mat));}


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
