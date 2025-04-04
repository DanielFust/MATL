#ifndef _MATL_IS_COLMAJOR_STORAGE_SPECIALIZTIONS__HPP
#define _MATL_IS_COLMAJOR_STORAGE_SPECIALIZATIONS__HPP

namespace matl 
{
   /***********************************************************************************************
    * class is_colmajor_storage
    * 
    * Type trait to determine if a matrix or expression is stored in row-major order.
    **********************************************************************************************/
    
    // Recursion
    template<typename FIRST_TYPE, typename ... OTHER_TYPES>
    struct is_colmajor_storage<FIRST_TYPE, OTHER_TYPES...>
    {
        static constexpr bool value = is_colmajor_storage<FIRST_TYPE>::value && 
                                       is_colmajor_storage<OTHER_TYPES...>::value;
    };
    

    // Default (false)
    template<typename T>
    struct is_colmajor_storage<T>
    {static constexpr bool value = false;};
    
    /****************** Interface specializations ******************/

    // Matrix Algebra Object
    template<typename MAO> 
    struct is_colmajor_storage< _MATRIX_ALGEBRA_OBJECT_<MAO> > 
    {static constexpr bool value = is_colmajor_storage<MAO>::value;}; 
    
    // Matrix
    template<typename MAT>
    struct is_colmajor_storage< _MATRIX_<MAT> >
    {static constexpr bool value = is_colmajor_storage<MAT>::value;};

    // Expression
    template<typename EXPR>
    struct is_colmajor_storage< _EXPRESSION_<EXPR> >
    {static constexpr bool value = is_colmajor_storage<EXPR>::value;};

    // Fixed-size interface
    template<typename FS>
    struct is_colmajor_storage< _FIXED_SIZE_<FS> >
    {static constexpr bool value = is_colmajor_storage<FS>::value;};

    // Fixed-size matrix interface
    template<typename FSM>
    struct is_colmajor_storage< _FIXED_SIZE_MATRIX_<FSM> >
    {static constexpr bool value = is_colmajor_storage<FSM>::value;};



    /********************** Matrix Specializations ***********************/

    // RowMajorMatrix
    template<class T>
    struct is_colmajor_storage< RowMajorMatrix<T> >
    {static constexpr bool value = true;};

    // FixedSizeRowMajorMatrix
    template<class T, index_t M, index_t N>
    struct is_colmajor_storage< FixedSizeRowMajorMatrix<T, M, N> >
    {static constexpr bool value = true;};


    /********************** Expression Specializations ***********************/

    // Matrix-Matrix Sum
    template<typename MAO1, typename MAO2>
    struct is_colmajor_storage< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >
    {static constexpr bool value = is_colmajor_storage<MAO1>::value && is_colmajor_storage<MAO2>::value;};

    // Matrix-Scalar Sum
    template<typename MAO, class T>
    struct is_colmajor_storage< _MATRIX_PLUS_SCALAR_<MAO, T> >
    {static constexpr bool value = is_colmajor_storage<MAO>::value;};

    // Matrix-Matrix Difference
    template<typename MAO1, typename MAO2>
    struct is_colmajor_storage< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >
    {static constexpr bool value = is_colmajor_storage<MAO1>::value && is_colmajor_storage<MAO2>::value;};

    // Matrix-Scalar Difference
    template<typename MAO, class T>
    struct is_colmajor_storage< _MATRIX_MINUS_SCALAR_<MAO, T> >
    {static constexpr bool value = is_colmajor_storage<MAO>::value;};

    // Scalar-Matrix Difference
    template<class T, typename MAO>
    struct is_colmajor_storage< _SCALAR_MINUS_MATRIX_<T, MAO> >
    {static constexpr bool value = is_colmajor_storage<MAO>::value;};

    // Matrix-Matrix Product
    template<typename MAO1, typename MAO2>
    struct is_colmajor_storage< _MATRIX_TIMES_MATRIX_<MAO1, MAO2> >
    {static constexpr bool value = is_colmajor_storage<MAO1>::value && is_colmajor_storage<MAO2>::value;};

    // Matrix-Scalar Product
    template<typename MAO, class T>
    struct is_colmajor_storage< _MATRIX_TIMES_SCALAR_<MAO, T> >
    {static constexpr bool value = is_colmajor_storage<MAO>::value;};

    // Matrix-Scalar Division
    template<typename MAO, class T>
    struct is_colmajor_storage< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >
    {static constexpr bool value = is_colmajor_storage<MAO>::value;};

    // Matrix transpose
    template<typename MAO>
    struct is_colmajor_storage< _MATRIX_TRANSPOSE_<MAO> >
    {static constexpr bool value = is_rowmajor_storage<MAO>::value;};

}

#endif