/***********************************************************************************************
 * struct member_type specializations
 * 
 * The member_type struct is used to specify the types of member variables in general classes
 * of matrix expressions used to store reference to or result of a matrix algebra object (MAO).
 * 
 * member types:
 *   elementwise_expression:
 *     - If the MAO implements the _ELEMENTWISE_ interface, the member type is a 'const&'
 *       reference to the MAO.
 *    
 *     - If the MAO does not implement the _ELEMENTWISE_ interface (e.g. matrix 
 *       multiplication), the member type is a _MATRIX_ type without the 'const&' qualifiers.
 *       (i.e. it is a container for the result of the expression owned by the constructing 
 *       expression object rather than a reference to a MAO).
 * 
 *   matmul_expression:
 *     - If the MAO is a _MATRIX_ type, the member type is a 'const&' reference to the MAO.
 *     
 *     - If the MAO is an expression, the member type is a _MATRIX_ type without the 'const&'.
 *       For cache efficiency and minimizing floating point ops in matrix multiplications and 
 *       similar operations, execution of preceding expressions are triggered and stored in 
 *       temporary matrices owned by the matrix multiplication expression object.
 *   
 ***********************************************************************************************/ 

#ifndef _MATL_MEMBER_TYPE_SPECIALIZATIONS__HPP
#define _MATL_MEMBER_TYPE_SPECIALIZATIONS__HPP

namespace matl 
{
    // Matrix Algebra Object
    template<typename MAO>
    struct member_type< _MATRIX_ALGEBRA_OBJECT_<MAO> >
    {
        typedef typename member_type<MAO>::elementwise_expression elementwise_expression;
        typedef typename member_type<MAO>::matmul_expression      matmul_expression;
    };

    // Elementwise interface
    template<typename EW>
    struct member_type< _ELEMENTWISE_<EW> >
    {
        typedef typename member_type<EW>::elementwise_expression elementwise_expression;
        typedef typename member_type<EW>::matmul_expression matmul_expression;
    };

    // Matrix base class
    template<typename MAT>
    struct member_type< _MATRIX_<MAT> >
    {
        typedef MAT const& elementwise_expression;
        typedef MAT const& matmul_expression;
    };

    // Expression base class
    template<typename EXPR>
    struct member_type< _EXPRESSION_<EXPR> >
    {
        typedef typename member_type<EXPR>::elementwise_expression elementwise_expression;
        typedef typename member_type<EXPR>::matmul_expression      matmul_expression;
    };

    // Fixed-size interface
    template<typename FS>
    struct member_type< _FIXED_SIZE_<FS> >
    {
        typedef typename member_type<FS>::elementwise_expression elementwise_expression;
        typedef typename member_type<FS>::matmul_expression      matmul_expression;
    };

    // Fixed-size matrix interface
    template<typename FSM>
    struct member_type< _FIXED_SIZE_MATRIX_<FSM> >
    {
        typedef typename member_type<FSM>::elementwise_expression elementwise_expression;
        typedef typename member_type<FSM>::matmul_expression      matmul_expression;
    };

    // RowMajorMatrix class
    template<class T>
    struct member_type< RowMajorMatrix<T> >
    {
        typedef RowMajorMatrix<T> const& elementwise_expression;
        typedef RowMajorMatrix<T> const& matmul_expression;
    };

    // FixedSizeRowMajorMatrix class
    template<class T, index_t M, index_t N>
    struct member_type< FixedSizeRowMajorMatrix<T, M, N> >
    {
        typedef FixedSizeRowMajorMatrix<T, M, N> const& elementwise_expression;
        typedef FixedSizeRowMajorMatrix<T, M, N> const& matmul_expression;
    };

    // Matrix-Matrix Sum class
    template<typename MAO1,typename MAO2>
    struct member_type< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >
    {
        typedef typename _MATRIX_PLUS_MATRIX_<MAO1, MAO2>::dtype dtype;
        typedef          _MATRIX_PLUS_MATRIX_<MAO1, MAO2>                             elementwise_expression;
        typedef typename is_fixed_size< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >::matrix_type matmul_expression;
    };

    // Matrix-Scalar Sum class
    template<typename MAO, class T>
    struct member_type< _MATRIX_PLUS_SCALAR_<MAO, T> >
    {
        typedef typename _MATRIX_PLUS_SCALAR_<MAO, T>::dtype dtype;
        typedef          _MATRIX_PLUS_SCALAR_<MAO, T>                              elementwise_expression;
        typedef typename is_fixed_size< _MATRIX_PLUS_SCALAR_<MAO, T> >::matrix_type matmul_expression;
    };

    // Matrix-Matrix Difference class
    template<typename MAO1,typename MAO2>
    struct member_type< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >
    {
        typedef typename _MATRIX_MINUS_MATRIX_<MAO1, MAO2>::dtype dtype;
        typedef          _MATRIX_MINUS_MATRIX_<MAO1, MAO2>                              elementwise_expression;
        typedef typename is_fixed_size< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >::matrix_type matmul_expression;
    };

    // Matrix-Scalar Difference class
    template<typename MAO, class T>
    struct member_type< _MATRIX_MINUS_SCALAR_<MAO, T> >
    {
        typedef typename _MATRIX_MINUS_SCALAR_<MAO, T>::dtype dtype;
        typedef          _MATRIX_MINUS_SCALAR_<MAO, T>                              elementwise_expression;
        typedef typename is_fixed_size< _MATRIX_MINUS_SCALAR_<MAO, T> >::matrix_type matmul_expression;
    };

    // Scalar-Matrix Difference class
    template<class T, typename MAO>
    struct member_type< _SCALAR_MINUS_MATRIX_<T, MAO> >
    {
        typedef typename _SCALAR_MINUS_MATRIX_<T, MAO>::dtype dtype;
        typedef          _SCALAR_MINUS_MATRIX_<T, MAO>                              elementwise_expression;
        typedef typename is_fixed_size< _SCALAR_MINUS_MATRIX_<T, MAO> >::matrix_type matmul_expression;
    };

    // Matrix-Matrix Product class
    template<typename MAO1,typename MAO2>
    struct member_type<_MATRIX_TIMES_MATRIX_<MAO1, MAO2> >
    {
        typedef typename _MATRIX_TIMES_MATRIX_<MAO1, MAO2>::dtype dtype;
        typedef typename is_fixed_size< _MATRIX_TIMES_MATRIX_<MAO1, MAO2> >::matrix_type elementwise_expression;
        typedef typename is_fixed_size< _MATRIX_TIMES_MATRIX_<MAO1, MAO2> >::matrix_type matmul_expression;
    };

    // Matrix-Scalar Product class
    template<typename MAO, class T>
    struct member_type< _MATRIX_TIMES_SCALAR_<MAO, T> >
    {
        typedef typename _MATRIX_TIMES_SCALAR_<MAO, T>::dtype       dtype;
        typedef          _MATRIX_TIMES_SCALAR_<MAO, T>                              elementwise_expression;
        typedef typename is_fixed_size< _MATRIX_TIMES_SCALAR_<MAO, T> >::matrix_type matmul_expression;
    };

    // Matrix-Scalar Division class
    template<typename MAO, class T>
    struct member_type< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >
    {
        typedef typename _MATRIX_DIVIDED_BY_SCALAR_<MAO, T>::dtype dtype;
        typedef          _MATRIX_DIVIDED_BY_SCALAR_<MAO, T>                              elementwise_expression;
        typedef typename is_fixed_size< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >::matrix_type matmul_expression;
    };

    // Matrix Transpose class
    template<typename MAO>
    struct member_type< _MATRIX_TRANSPOSE_<MAO> >
    {
        typedef typename _MATRIX_TRANSPOSE_<MAO>::dtype dtype;
        typedef          _MATRIX_TRANSPOSE_<MAO>                              elementwise_expression;
        typedef typename is_fixed_size< _MATRIX_TRANSPOSE_<MAO> >::matrix_type matmul_expression;
    };
    // Owned matrix class (by value: copying it copies a shared pointer, not the matrix). Products
    // evaluate the owned matrix into their own matrix operand.
    template<typename MAT>
    struct member_type< _OWNED_MATRIX_<MAT> >
    {
        typedef typename _OWNED_MATRIX_<MAT>::dtype dtype;
        typedef          _OWNED_MATRIX_<MAT>        elementwise_expression;
        typedef typename is_fixed_size< _OWNED_MATRIX_<MAT> >::matrix_type matmul_expression;
    };
}
#endif