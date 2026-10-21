
/***********************************************************************************************
 * struct is_elementwise specializations
 * 
 * The is_elementwise struct is used to determine whether a matrix algebra object (MAO) implements
 * the _ELEMENTWISE_ interface, which is used for flow control at compile time.
 * 
 * _ELEMENTWISE_ interface:
 *     Determines if a given object implements the _ELEMENTWISE_ interface for
 *     memory-safe, double-indexed data retrieval.
 ***********************************************************************************************/ 
#ifndef _MATL_IS_ELEMENTWISE_SPECIALIZATIONS__HPP
#define _MATL_IS_ELEMENTWISE_SPECIALIZATIONS__HPP
namespace matl
{
    // Default for no specialization
    template<typename T>
    struct is_elementwise<T>
    {static constexpr bool value = false;};

    // Matrix Algebra Object
    template<typename MAO>
    struct is_elementwise< _MATRIX_ALGEBRA_OBJECT_<MAO> >
    {static constexpr bool value = is_elementwise<MAO>::value;};

    // Elementwise interface
    template<typename EW>
    struct is_elementwise< _ELEMENTWISE_<EW> >
    {static constexpr bool value = true;};

    // Matrix base class
    template<typename MAT>
    struct is_elementwise< _MATRIX_<MAT> >
    {static constexpr bool value = true;};

    // Expression base class
    template<typename EXPR>
    struct is_elementwise< _EXPRESSION_<EXPR> >
    {static constexpr bool value = is_elementwise<EXPR>::value;};

    // Fixed-size interface
    template<typename FS>
    struct is_elementwise< _FIXED_SIZE_<FS> >
    {static constexpr bool value = is_elementwise<FS>::value;};

    // Fixed-size matrix interface
    template<typename FSM>
    struct is_elementwise< _FIXED_SIZE_MATRIX_<FSM> >
    {static constexpr bool value = true;};

    // RowMajorMatrix class
    template<class T>
    struct is_elementwise< RowMajorMatrix<T> >
    {static constexpr bool value = true;};

    // FixedSizeRowMajorMatrix class
    template<class T, index_t M, index_t N>
    struct is_elementwise< FixedSizeRowMajorMatrix<T, M, N> >
    {static constexpr bool value = true;};

    // Matrix-Matrix Sum class
    template<typename MAO1,typename MAO2>
    struct is_elementwise< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >
    {static constexpr bool value = true;};

    // Matrix-Scalar Sum class
    template<typename MAO, class T>
    struct is_elementwise< _MATRIX_PLUS_SCALAR_<MAO, T> >
    {static constexpr bool value = true;};

    // Matrix-Matrix Difference class
    template<typename MAO1,typename MAO2>
    struct is_elementwise< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >
    {static constexpr bool value = true;};

    // Matrix-Scalar Difference class
    template<typename MAO, class T>
    struct is_elementwise< _MATRIX_MINUS_SCALAR_<MAO, T> >
    {static constexpr bool value = true;};

    // Scalar-Matrix Difference class
    template<class T, typename MAO>
    struct is_elementwise< _SCALAR_MINUS_MATRIX_<T, MAO> >
    {static constexpr bool value = true;};

    // Matrix-Matrix Product class
    template<typename MAO1,typename MAO2>
    struct is_elementwise<_MATRIX_TIMES_MATRIX_<MAO1, MAO2> >
    {static constexpr bool value = false;};

    // Matrix-Scalar Product class
    template<typename MAO, class T>
    struct is_elementwise< _MATRIX_TIMES_SCALAR_<MAO, T> >
    {static constexpr bool value = true;};

    // Matrix-Scalar Division class
    template<typename MAO, class T>
    struct is_elementwise< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >
    {static constexpr bool value = true;};

    // Matrix Transpose class
    template<typename MAO>
    struct is_elementwise< _MATRIX_TRANSPOSE_<MAO> >
    {static constexpr bool value = true;};

    // Owned matrix class
    template<typename MAT>
    struct is_elementwise< _OWNED_MATRIX_<MAT> >
    {static constexpr bool value = true;};
}
#endif