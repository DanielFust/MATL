#ifndef _MATL_IS_TRANSPOSED_SPECIALIZATIONS__HPP
#define _MATL_IS_TRANSPOSED_SPECIALIZATIONS__HPP

namespace matl 
{
   /***********************************************************************************************
    * struct is_transposed
    * 
    * The 'is_transposed' struct provides compile-time logic to identify matrix and expression 
    * transposes. This is used to determine if an expression may be safely evaluated to a target
    * matrix without memory corruption, or if a temporary matrix is required.
    * 
    * Note: this struc is trivially simple, but greatly improves readability of the compile-time
    *       branching logic in the expression classes.
    *************************************************************************************************/

    // Default
    template<typename T>
    struct is_transposed<T>
    {static constexpr bool value = false;};

    // Matrix transpose
    template<typename MAO>
    struct is_transposed<_MATRIX_TRANSPOSE_<MAO> >
    {static constexpr bool value = true;};
    

    /********* Interfaces (downcasted) *********/

    // Matrix algebra object
    template<typename MAO>
    struct is_transposed<_MATRIX_ALGEBRA_OBJECT_<MAO> >
    {static constexpr bool value = is_transposed<MAO>::value;};

    // Expression
    template<typename EXPR>
    struct is_transposed<_EXPRESSION_<EXPR> >
    {static constexpr bool value = is_transposed<EXPR>::value;};

    // Elementwise expression
    template<typename EW>
    struct is_transposed<_ELEMENTWISE_<EW> >
    {static constexpr bool value = is_transposed<EW>::value;};

    // Fixed-size interface
    template<typename FS>
    struct is_transposed<_FIXED_SIZE_<FS> >
    {static constexpr bool value = is_transposed<FS>::value;};

}

#endif