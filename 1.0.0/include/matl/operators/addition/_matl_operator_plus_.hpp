/**************************************************************************************************
 * File: _matl_operator_plus_.hpp
 * 
 * Description: Implementations for the (+) operator for constructing _EXPRESSION_ objects
 **************************************************************************************************/
#ifndef _MATL_OPERATOR_PLUS__HPP
#define _MATL_OPERATOR_PLUS__HPP

namespace matl
{
    // Matrix + Matrix
    template<typename MAT1, typename MAT2>
    _MATRIX_PLUS_MATRIX_<MAT1, MAT2> operator+ (_MATRIX_<MAT1> const& mat1, _MATRIX_<MAT2> const& mat2)
    {return _MATRIX_PLUS_MATRIX_<MAT1, MAT2>(mat1, mat2);}

    // Matrix + Expression
    template<typename MAT, typename EXPR>
    _MATRIX_PLUS_MATRIX_<MAT, _EXPRESSION_<EXPR> > operator+ (_MATRIX_<MAT> const& mat, _EXPRESSION_<EXPR> const& expr)
    {return _MATRIX_PLUS_MATRIX_<MAT, _EXPRESSION_<EXPR> >(mat, expr);}

    // Expression + Matrix
    template<typename EXPR, typename MAT>
    _MATRIX_PLUS_MATRIX_<_EXPRESSION_<EXPR>, MAT> operator+ (_EXPRESSION_<EXPR> const& expr, _MATRIX_<MAT> const& mat)
    {return _MATRIX_PLUS_MATRIX_<_EXPRESSION_<EXPR>, MAT>(expr, mat);}

    // Expression + Expression
    template<typename EXPR1, typename EXPR2>
    _MATRIX_PLUS_MATRIX_<_EXPRESSION_<EXPR1>, _EXPRESSION_<EXPR2> > operator+ (_EXPRESSION_<EXPR1> const& expr1, _EXPRESSION_<EXPR2> const& expr2)
    {return _MATRIX_PLUS_MATRIX_<_EXPRESSION_<EXPR1>, _EXPRESSION_<EXPR2> >(expr1, expr2);}

    // Matrix + Scalar
    template<typename MAT, class T, enable_if_scalar<T> = true>
    _MATRIX_PLUS_SCALAR_<MAT, T> operator+ (_MATRIX_<MAT> const& mat, T scalar)
    {return _MATRIX_PLUS_SCALAR_<MAT, T>(mat, scalar);}

    // Expression + Scalar
    template<typename EXPR, class T, enable_if_scalar<T> = true>
    _MATRIX_PLUS_SCALAR_<_EXPRESSION_<EXPR>, T> operator+ (_EXPRESSION_<EXPR> const& expr, T scalar)
    {return _MATRIX_PLUS_SCALAR_<_EXPRESSION_<EXPR>, T>(expr, scalar);}

    // Scalar + Matrix
    template<class T, typename MAT, enable_if_scalar<T> = true>
    _MATRIX_PLUS_SCALAR_<MAT, T> operator+ (T scalar, _MATRIX_<MAT> const& mat)
    {return _MATRIX_PLUS_SCALAR_<MAT, T>(mat, scalar);}

    // Scalar + Expression
    template<class T, typename EXPR, enable_if_scalar<T> = true>
    _MATRIX_PLUS_SCALAR_<_EXPRESSION_<EXPR>, T> operator+ (T scalar, _EXPRESSION_<EXPR> const& expr)
    {return _MATRIX_PLUS_SCALAR_<_EXPRESSION_<EXPR>, T>(expr, scalar);}
}

#endif