#ifndef _MATL_OPERATOR_DIVISION__HPP
#define _MATL_OPERATOR_DIVISION__HPP

namespace matl
{
    // Matrix divided by non-integer scalar
    template<typename MAT, class T, enable_if_scalar<T> = true, disable_if_integral<T> = true>
    _MATRIX_TIMES_SCALAR_<MAT, typename type_traits< _MATRIX_TIMES_SCALAR_<MAT, T> >::dtype> 
    operator/ (_MATRIX_<MAT> const& mat, T scalar)
    {
        typename type_traits< _MATRIX_TIMES_SCALAR_<MAT, T> >::dtype reciprocal = 1.0 / scalar;
        return _MATRIX_TIMES_SCALAR_<MAT, T>(mat, reciprocal);
    }

    // Matrix divided by integer scalar
    template<typename MAT, class T, enable_if_scalar<T> = true, enable_if_integral<T> = true>
    _MATRIX_DIVIDED_BY_SCALAR_<MAT, T> operator/ (_MATRIX_<MAT> const& mat, T scalar)
    {return _MATRIX_DIVIDED_BY_SCALAR_<MAT, T>(mat, scalar);}

    // Expression divided by non-integer scalar
    template<typename EXPR, class T, enable_if_scalar<T> = true, disable_if_integral<T> = true>
    _MATRIX_TIMES_SCALAR_<EXPR, typename type_traits< _MATRIX_TIMES_SCALAR_<EXPR, T> >::dtype> 
    operator/ (_EXPRESSION_<EXPR> const& expr, T scalar)
    {
        typename type_traits< _MATRIX_TIMES_SCALAR_<EXPR, T> >::dtype reciprocal = 1.0 / scalar;
        return _MATRIX_TIMES_SCALAR_<EXPR, T>(expr, reciprocal);
    }

    // Expression divided by integer scalar
    template<typename EXPR, class T, enable_if_scalar<T> = true, enable_if_integral<T> = true>
    _MATRIX_DIVIDED_BY_SCALAR_<EXPR, T> operator/ (_EXPRESSION_<EXPR> const& expr, T scalar)
    {return _MATRIX_DIVIDED_BY_SCALAR_<EXPR, T>(expr, scalar);}

}


#endif