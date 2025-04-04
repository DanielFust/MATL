#ifndef _MATL_EXPRESSION_TYPES_DECLARATIONS__HPP
#define _MATL_EXPRESSION_TYPES_DECLARATIONS__HPP
namespace matl
{
    // General Matrix + Matrix summations
    template<typename MAO1, typename MAO2> class _MATRIX_PLUS_MATRIX_;

    // General Matrix and scalar summations
    template<typename MAO, class T> class _MATRIX_PLUS_SCALAR_;

    // General Matrix-Matrix subtraction
    template<typename MAO1, typename MAO2> class _MATRIX_MINUS_MATRIX_;

    // General Matrix and scalar subtraction
    template<typename MAO, class T> class _MATRIX_MINUS_SCALAR_;
    template<class T, typename MAO> class _SCALAR_MINUS_MATRIX_;

    // General Matrix * scalar multiplications
    template<typename MAO, class T> class _MATRIX_TIMES_SCALAR_;

    // General Matrix * Matrix multiplications
    template<typename MAT1, typename MAT2> class _MATRIX_TIMES_MATRIX_;

    // General Matrix / scalar divisions
    template<typename MAO, class T> class _MATRIX_DIVIDED_BY_SCALAR_;

    // Matrix transpose
    template<typename MAO> class _MATRIX_TRANSPOSE_;

}
#endif
