/**************************************************************************************************
 * file: _matl_is_matrix_algebra_object_specializations_.hpp
 * 
 * The is_matrix_algebra_object struct provides compile-time logic to determine whether a matrix 
 * algebra object (MAO) from the MATL library or not. This is used to clearly indicate what
 * types of data are allowed in constructors and operators, and to reduce duplication of code
 * when specializing operators while disallowing unintended types.
 **************************************************************************************************/ 
#ifndef _MATL_IS_MATRIX_ALGEBRA_OBJECT_SPECIALIZATIONS__HPP
#define _MATL_IS_MATRIX_ALGEBRA_OBJECT_SPECIALIZATIONS__HPP

namespace matl
{
    // Default for no specialization
    template<typename T>
    struct is_matrix_algebra_object<T>
    {static constexpr bool value = false;};

    // Matrix Algebra Object
    template<typename MAO>
    struct is_matrix_algebra_object< _MATRIX_ALGEBRA_OBJECT_<MAO> >
    {static constexpr bool value = true;};

    // Elementwise interface
    template<typename EW>
    struct is_matrix_algebra_object< _ELEMENTWISE_<EW> >
    {static constexpr bool value = true;};

    // Matrix base class
    template<typename MAT>
    struct is_matrix_algebra_object< _MATRIX_<MAT> >
    {static constexpr bool value = true;};

    // Expression base class
    template<typename EXPR>
    struct is_matrix_algebra_object< _EXPRESSION_<EXPR> >
    {static constexpr bool value = is_matrix_algebra_object<EXPR>::value;};

    // Fixed-size interface
    template<typename FS>
    struct is_matrix_algebra_object< _FIXED_SIZE_<FS> >
    {static constexpr bool value = is_matrix_algebra_object<FS>::value;};

    // Fixed-size matrix interface
    template<typename FSM>
    struct is_matrix_algebra_object< _FIXED_SIZE_MATRIX_<FSM> >
    {static constexpr bool value = true;};

    // RowMajorMatrix class
    template<class T>
    struct is_matrix_algebra_object< RowMajorMatrix<T> >
    {static constexpr bool value = true;};
    
    // FixedSizeRowMajorMatrix class
    template<class T, index_t M, index_t N>
    struct is_matrix_algebra_object< FixedSizeRowMajorMatrix<T, M, N> >
    {static constexpr bool value = true;};

    // Matrix-Matrix Sum class
    template<typename MAO1,typename MAO2>
    struct is_matrix_algebra_object< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >
    {static constexpr bool value = true;};

    // Matrix-Scalar Sum class
    template<typename MAO, class T>
    struct is_matrix_algebra_object< _MATRIX_PLUS_SCALAR_<MAO, T> >
    {static constexpr bool value = true;};

    // Matrix-Matrix Difference class
    template<typename MAO1,typename MAO2>
    struct is_matrix_algebra_object< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >
    {static constexpr bool value = true;};

    // Matrix-Scalar Difference class
    template<typename MAO, class T>
    struct is_matrix_algebra_object< _MATRIX_MINUS_SCALAR_<MAO, T> >
    {static constexpr bool value = true;};

    // Scalar-Matrix Difference class
    template<class T, typename MAO>
    struct is_matrix_algebra_object< _SCALAR_MINUS_MATRIX_<T, MAO> >
    {static constexpr bool value = true;};

    // Matrix-Matrix Product class
    template<typename MAO1,typename MAO2>
    struct is_matrix_algebra_object< _MATRIX_TIMES_MATRIX_<MAO1, MAO2> >
    {static constexpr bool value = true;};
   
    // Matrix-Scalar Product class
    template<typename MAO, class T>
    struct is_matrix_algebra_object< _MATRIX_TIMES_SCALAR_<MAO, T> >
    {static constexpr bool value = true;};

    // Matrix-Scalar Division class
    template<typename MAO, class T>
    struct is_matrix_algebra_object< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >
    {static constexpr bool value = true;};

    // Matrix Transpose class
    template<typename MAO>
    struct is_matrix_algebra_object< _MATRIX_TRANSPOSE_<MAO> >
    {static constexpr bool value = true;};

    // Owned matrix class
    template<typename MAT>
    struct is_matrix_algebra_object< _OWNED_MATRIX_<MAT> >
    {static constexpr bool value = true;};

}

#endif