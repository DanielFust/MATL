#ifndef _MATL_TYPE_TRAITS_SPECIALIZATIONS__HPP
#define _MATL_TYPE_TRAITS_SPECIALIZATIONS__HPP

namespace matl
{
    // Matrix Algebra Object return types
    template<typename MAO>
    struct type_traits< _MATRIX_ALGEBRA_OBJECT_<MAO> >
    {
        typedef typename type_traits<MAO>::dtype            dtype;
        typedef typename type_traits<MAO>::unqualified_type unqualified_type;
    };

    // Elementwise interface return types
    template<typename EW>
    struct type_traits< _ELEMENTWISE_<EW> >
    {
        typedef typename type_traits<EW>::dtype            dtype;
        typedef typename type_traits<EW>::unqualified_type unqualified_type;
    };


    // Matrix base class return types
    template<typename MAT>
    struct type_traits< _MATRIX_<MAT> >
    {
        typedef typename type_traits<MAT>::dtype            dtype;
        typedef typename type_traits<MAT>::unqualified_type unqualified_type;
        
    };

    // Expression base class return types
    template<typename EXPR>
    struct type_traits< _EXPRESSION_<EXPR> >
    {
        typedef typename type_traits<EXPR>::dtype            dtype;
        typedef typename type_traits<EXPR>::unqualified_type unqualified_type;
    };

    // Fixed-size interface return types
    template<typename FS>
    struct type_traits< _FIXED_SIZE_<FS> >
    {
        typedef typename type_traits<FS>::dtype            dtype;
        typedef typename type_traits<FS>::unqualified_type unqualified_type;
    };

    // Fixed-size matrix interface return types
    template<typename FSM>
    struct type_traits< _FIXED_SIZE_MATRIX_<FSM> >
    {
        typedef typename type_traits<FSM>::dtype            dtype;
        typedef typename type_traits<FSM>::unqualified_type unqualified_type;
    };

    // RowMajorMatrix class return types
    template<class T>
    struct type_traits< RowMajorMatrix<T> >
    {
        typedef T                 dtype;
        typedef RowMajorMatrix<T> unqualified_type;
    };

    // FixedSizeRowMajorMatrix class return types
    template<class T, index_t M, index_t N>
    struct type_traits< FixedSizeRowMajorMatrix<T, M, N> >
    {
        typedef T dtype;
        typedef FixedSizeRowMajorMatrix<T, M, N> unqualified_type;
    };

    // Matrix-Matrix Sum class return types
    template<typename MAO1,typename MAO2>
    struct type_traits< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >
    {
        typedef typename promote<typename type_traits<MAO1>::dtype, 
                                 typename type_traits<MAO2>::dtype>::type dtype;
        typedef          _MATRIX_PLUS_MATRIX_<MAO1,MAO2>                  unqualified_type;
    };

    // Matrix-Scalar Sum class return types
    template<typename MAO, class T>
    struct type_traits< _MATRIX_PLUS_SCALAR_<MAO, T> >
    {
        typedef typename promote<typename type_traits<MAO>::dtype, T>::type dtype;
        typedef          _MATRIX_PLUS_SCALAR_<MAO,T>                        unqualified_type;
    };

    // Matrix-Matrix Difference class return types
    template<typename MAO1,typename MAO2>
    struct type_traits< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >
    {
        typedef typename promote<typename type_traits<MAO1>::dtype, 
                                 typename type_traits<MAO2>::dtype>::type dtype;
        typedef          _MATRIX_MINUS_MATRIX_<MAO1,MAO2>                 unqualified_type;
    };

    // Matrix-Scalar Difference class return types
    template<typename MAO, class T>
    struct type_traits< _MATRIX_MINUS_SCALAR_<MAO, T> >
    {
        typedef typename promote<typename type_traits<MAO>::dtype, T>::type dtype;
        typedef          _MATRIX_MINUS_SCALAR_<MAO,T>                       unqualified_type;
    };

    // Scalar-Matrix Difference class return types
    template<class T, typename MAO>
    struct type_traits< _SCALAR_MINUS_MATRIX_<T, MAO> >
    {
        typedef typename promote<T, typename type_traits<MAO>::dtype>::type dtype;
        typedef          _SCALAR_MINUS_MATRIX_<T,MAO>                       unqualified_type;
    };

    // Matrix-Scalar Product class return types
    template<typename MAO, class T>
    struct type_traits< _MATRIX_TIMES_SCALAR_<MAO, T> >
    {
        typedef typename promote<typename type_traits<MAO>::dtype, T>::type dtype;
        typedef          _MATRIX_TIMES_SCALAR_<MAO, T>                      unqualified_type;
    };

    // Matrix-Matrix Product class return types
    template<typename MAO1,typename MAO2>
    struct type_traits<_MATRIX_TIMES_MATRIX_<MAO1, MAO2> >
    {
        typedef typename promote<typename type_traits<MAO1>::dtype,
                                 typename type_traits<MAO2>::dtype>::type dtype;
        typedef          _MATRIX_TIMES_MATRIX_<MAO1,MAO2>                 unqualified_type;
    };


    // Matrix-Scalar Division class return types
    template<typename MAO, class T>
    struct type_traits< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >
    {
        typedef typename promote<typename type_traits<MAO>::dtype, T>::type                         dtype;
        typedef          _MATRIX_DIVIDED_BY_SCALAR_<typename type_traits<MAO>::unqualified_type, T> unqualified_type;
    };


    // Matrix transpose class return types
    template<typename MAO>
    struct type_traits< _MATRIX_TRANSPOSE_<MAO> >
    {
        typedef typename type_traits<MAO>::dtype dtype;
        typedef          _MATRIX_TRANSPOSE_<MAO> unqualified_type;
    };
    
} 
#endif
