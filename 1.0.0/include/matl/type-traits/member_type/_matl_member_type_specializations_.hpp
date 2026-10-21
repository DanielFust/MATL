/***********************************************************************************************
 * struct member_type specializations
 * 
 * The member_type struct is used to specify the types of member variables in general classes
 * of matrix expressions used to store reference to or result of a matrix algebra object (MAO).
 * 
 * member types:
 *   elementwise_expression (how an expression stores an operand of this type):
 *     - Matrices are stored as a 'const&' reference, so an expression is a live formula that
 *       reads the matrices' current values when it is evaluated.
 *    
 *     - Expressions (including matrix products) are stored by value, so an expression object
 *       remains valid after the temporaries it was built from are destroyed. Expressions are
 *       small (references, scalars and nested expressions); they hold no results. Products and
 *       other intermediate results are computed into temporaries only during evaluation (see
 *       expression-types/evaluation/_matl_evaluation_.hpp).
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
    };

    // Elementwise interface
    template<typename EW>
    struct member_type< _ELEMENTWISE_<EW> >
    {
        typedef typename member_type<EW>::elementwise_expression elementwise_expression;
    };

    // Matrix base class
    template<typename MAT>
    struct member_type< _MATRIX_<MAT> >
    {
        typedef MAT const& elementwise_expression;
    };

    // Expression base class
    template<typename EXPR>
    struct member_type< _EXPRESSION_<EXPR> >
    {
        typedef typename member_type<EXPR>::elementwise_expression elementwise_expression;
    };

    // Fixed-size interface
    template<typename FS>
    struct member_type< _FIXED_SIZE_<FS> >
    {
        typedef typename member_type<FS>::elementwise_expression elementwise_expression;
    };

    // Fixed-size matrix interface
    template<typename FSM>
    struct member_type< _FIXED_SIZE_MATRIX_<FSM> >
    {
        typedef typename member_type<FSM>::elementwise_expression elementwise_expression;
    };

    // RowMajorMatrix class
    template<class T>
    struct member_type< RowMajorMatrix<T> >
    {
        typedef RowMajorMatrix<T> const& elementwise_expression;
    };

    // FixedSizeRowMajorMatrix class
    template<class T, index_t M, index_t N>
    struct member_type< FixedSizeRowMajorMatrix<T, M, N> >
    {
        typedef FixedSizeRowMajorMatrix<T, M, N> const& elementwise_expression;
    };

    // Matrix-Matrix Sum class
    template<typename MAO1,typename MAO2>
    struct member_type< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >
    {
        typedef typename _MATRIX_PLUS_MATRIX_<MAO1, MAO2>::dtype dtype;
        typedef          _MATRIX_PLUS_MATRIX_<MAO1, MAO2>                             elementwise_expression;
    };

    // Matrix-Scalar Sum class
    template<typename MAO, class T>
    struct member_type< _MATRIX_PLUS_SCALAR_<MAO, T> >
    {
        typedef typename _MATRIX_PLUS_SCALAR_<MAO, T>::dtype dtype;
        typedef          _MATRIX_PLUS_SCALAR_<MAO, T>                              elementwise_expression;
    };

    // Matrix-Matrix Difference class
    template<typename MAO1,typename MAO2>
    struct member_type< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >
    {
        typedef typename _MATRIX_MINUS_MATRIX_<MAO1, MAO2>::dtype dtype;
        typedef          _MATRIX_MINUS_MATRIX_<MAO1, MAO2>                              elementwise_expression;
    };

    // Matrix-Scalar Difference class
    template<typename MAO, class T>
    struct member_type< _MATRIX_MINUS_SCALAR_<MAO, T> >
    {
        typedef typename _MATRIX_MINUS_SCALAR_<MAO, T>::dtype dtype;
        typedef          _MATRIX_MINUS_SCALAR_<MAO, T>                              elementwise_expression;
    };

    // Scalar-Matrix Difference class
    template<class T, typename MAO>
    struct member_type< _SCALAR_MINUS_MATRIX_<T, MAO> >
    {
        typedef typename _SCALAR_MINUS_MATRIX_<T, MAO>::dtype dtype;
        typedef          _SCALAR_MINUS_MATRIX_<T, MAO>                              elementwise_expression;
    };

    // Matrix-Matrix Product class
    template<typename MAO1,typename MAO2>
    struct member_type<_MATRIX_TIMES_MATRIX_<MAO1, MAO2> >
    {
        typedef typename _MATRIX_TIMES_MATRIX_<MAO1, MAO2>::dtype dtype;
        typedef          _MATRIX_TIMES_MATRIX_<MAO1, MAO2>                              elementwise_expression;
    };

    // Matrix-Scalar Product class
    template<typename MAO, class T>
    struct member_type< _MATRIX_TIMES_SCALAR_<MAO, T> >
    {
        typedef typename _MATRIX_TIMES_SCALAR_<MAO, T>::dtype       dtype;
        typedef          _MATRIX_TIMES_SCALAR_<MAO, T>                              elementwise_expression;
    };

    // Matrix-Scalar Division class
    template<typename MAO, class T>
    struct member_type< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >
    {
        typedef typename _MATRIX_DIVIDED_BY_SCALAR_<MAO, T>::dtype dtype;
        typedef          _MATRIX_DIVIDED_BY_SCALAR_<MAO, T>                              elementwise_expression;
    };

    // Matrix Transpose class
    template<typename MAO>
    struct member_type< _MATRIX_TRANSPOSE_<MAO> >
    {
        typedef typename _MATRIX_TRANSPOSE_<MAO>::dtype dtype;
        typedef          _MATRIX_TRANSPOSE_<MAO>                              elementwise_expression;
    };
    // Owned matrix class (by value: copying it copies a shared pointer, not the matrix)
    template<typename MAT>
    struct member_type< _OWNED_MATRIX_<MAT> >
    {
        typedef typename _OWNED_MATRIX_<MAT>::dtype dtype;
        typedef          _OWNED_MATRIX_<MAT>        elementwise_expression;
    };
}
#endif