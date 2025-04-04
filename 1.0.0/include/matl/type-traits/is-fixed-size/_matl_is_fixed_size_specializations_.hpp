/**************************************************************************************************
 * struct matl::fixed_size
 * 
 * The is_fixed_size struct is used to determine whether a matrix algebra object (MAO) implements
 * the _FIXED_SIZE_ interface, i.e. whether an object's dimensions can be determined at compile 
 * time. This can be used to change the member types of expressions to be fixed size matrices or
 * matrix/expression references, thereby avoiding unnecessary heap allocations or recycling of
 * memory by the memory manager.
 * 
 * The is_fixed_size struct has three members and two type definitions for all specializations:
 *   - value:       A boolean value indicating whether the MAO implements the _FIXED_SIZE_ 
 *                  interface.
 * 
 *   - ROWS:        The number of rows in the MAO. If the MAO is not fixed-size, this will be 0.
 * 
 *   - COLS:        The number of columns in the MAO. If the MAO is not fixed-size, this will be 0.
 * 
 *   - dtype:       The datatype of the MAO. This is redundant with the type_traits struct, but
 *                  is included to help readability when utilizing the selector structs.
 * 
 *   - matrix_type: The typename of the matrix implied by the MAO. For expressions without
 *                  implied storage ordering, the matrix_type will be either the dynamic or
 *                  fixed-size variant of the MATL default Matrix type. 
 **************************************************************************************************/
#ifndef _MATL_IS_FIXED_SIZE_SPECIALIZATIONS__HPP
#define _MATL_IS_FIXED_SIZE_SPECIALIZATIONS__HPP

namespace matl 
{
    namespace is_fixed_size_selectors
    {
       /*******************************************************************************************
        * struct default_matrix_selector
        * 
        * This helper struct contains a matrix typdef that is either the dynamic or fixed-size
        * variant of the MATL default Matrix type depending on wherer the Matrix Algebra Object
        * is fixed-size or not.
        *******************************************************************************************/
        template<bool IS_FIXED_SIZE, index_t ROWS, index_t COLS, class dtype> 
        struct default_matrix_selector;
        
        // Specialization for fixed-size object
        template<index_t ROWS, index_t COLS, class dtype>
        struct default_matrix_selector<true, ROWS, COLS, dtype>
        {typedef Matrix<dtype,ROWS,COLS> matrix_type;};

        // Specialization for dynamic object
        template<index_t ROWS, index_t COLS, class dtype>
        struct default_matrix_selector<false, ROWS, COLS, dtype>
        {typedef Matrix<dtype> matrix_type;};
    }    

    // Default for no specialization
    template<typename T>
    struct is_fixed_size<T>
    {
        static constexpr bool    value = false;
        static constexpr index_t ROWS  = 0;
        static constexpr index_t COLS  = 0;
        typedef T                 dtype;
        //typedef T                 matrix_type; // Matrix type does not exist
    };

    // Matrix Algebra Object
    template<typename MAO>
    struct is_fixed_size< _MATRIX_ALGEBRA_OBJECT_<MAO> >
    {
        static constexpr bool    value = is_fixed_size<MAO>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAO>::ROWS;
        static constexpr index_t COLS  = is_fixed_size<MAO>::COLS;
        typedef typename type_traits<MAO>::dtype      dtype;
        typedef typename is_fixed_size<MAO>::matrix_type matrix_type;
    };

    // Elementwise interface
    template<typename EW>
    struct is_fixed_size< _ELEMENTWISE_<EW> >
    {
        static constexpr bool    value = is_fixed_size<EW>::value;
        static constexpr index_t ROWS  = is_fixed_size<EW>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<EW>::COLS*value;
        typedef typename type_traits<EW>::dtype      dtype;
        typedef typename is_fixed_size<EW>::matrix_type matrix_type;
    };

    // Matrix base class
    template<typename MAT>
    struct is_fixed_size< _MATRIX_<MAT> >
    {
        static constexpr bool    value = is_fixed_size<MAT>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAT>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<MAT>::COLS*value;
        typedef typename type_traits<MAT>::dtype         dtype;
        typedef typename is_fixed_size<MAT>::matrix_type matrix_type;
    };

    // Expression base class
    template<typename EXPR>
    struct is_fixed_size< _EXPRESSION_<EXPR> >
    {
        static constexpr bool    value = is_fixed_size<EXPR>::value;
        static constexpr index_t ROWS  = is_fixed_size<EXPR>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<EXPR>::COLS*value;
        typedef typename type_traits<EXPR>::dtype         dtype;
        typedef typename is_fixed_size<EXPR>::matrix_type matrix_type;
    };

    // Fixed-size interface
    template<typename FS>
    struct is_fixed_size< _FIXED_SIZE_<FS> >
    {
        static constexpr bool    value = true;
        static constexpr index_t ROWS  = FS::ROWS*value;
        static constexpr index_t COLS  = FS::COLS*value;
        typedef typename type_traits<FS>::dtype         dtype;
        typedef typename is_fixed_size<FS>::matrix_type matrix_type;
    };

    // Fixed-size matrix interface
    template<typename FSM>
    struct is_fixed_size< _FIXED_SIZE_MATRIX_<FSM> >
    {
        static constexpr bool    value = true;
        static constexpr index_t ROWS  = FSM::ROWS*value;
        static constexpr index_t COLS  = FSM::COLS*value;
        typedef typename type_traits<FSM>::dtype         dtype;
        typedef typename is_fixed_size<FSM>::matrix_type matrix_type;
    };

    // RowMajorMatrix class
    template<class T>
    struct is_fixed_size< RowMajorMatrix<T> >
    {
        static constexpr bool    value = false;
        static constexpr index_t ROWS  = 0;
        static constexpr index_t COLS  = 0;
        typedef T                 dtype;
        typedef RowMajorMatrix<T> matrix_type;
    };

    // FixedSizeRowMajorMatrix class
    template<class T, index_t M, index_t N>
    struct is_fixed_size< FixedSizeRowMajorMatrix<T, M, N> >
    {
        static constexpr bool    value = true;
        static constexpr index_t ROWS  = M;
        static constexpr index_t COLS  = N;
        typedef T                                dtype;
        typedef FixedSizeRowMajorMatrix<T, M, N> matrix_type;
    };

    // Matrix-Matrix Sum class
    template<typename MAO1,typename MAO2>
    struct is_fixed_size< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >
    {
        static constexpr bool    value = is_fixed_size<MAO1>::value && is_fixed_size<MAO2>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAO1>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<MAO2>::COLS*value;
        typedef promote<typename type_traits<MAO1>::dtype,
                        typename type_traits<MAO1>::dtype>::type                                              dtype;
        typedef typename is_fixed_size_selectors::default_matrix_selector<value,ROWS,COLS,dtype>::matrix_type matrix_type;
    };

    // Matrix-Scalar Sum class
    template<typename MAO, class T>
    struct is_fixed_size< _MATRIX_PLUS_SCALAR_<MAO, T> >
    {
        static constexpr bool    value = is_fixed_size<MAO>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAO>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<MAO>::COLS*value;
        typedef typename promote<typename type_traits<MAO>::dtype,T>::type                                    dtype;
        typedef typename is_fixed_size_selectors::default_matrix_selector<value,ROWS,COLS,dtype>::matrix_type matrix_type;
    };

    // Matrix-Matrix Difference class
    template<typename MAO1,typename MAO2>
    struct is_fixed_size< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >
    {
        static constexpr bool    value = is_fixed_size<MAO1>::value && is_fixed_size<MAO2>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAO1>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<MAO2>::COLS*value;
        typedef promote<typename type_traits<MAO1>::dtype,
                        typename type_traits<MAO1>::dtype>::type                                              dtype;
        typedef typename is_fixed_size_selectors::default_matrix_selector<value,ROWS,COLS,dtype>::matrix_type matrix_type;
    };

    // Matrix-Scalar Difference class
    template<typename MAO, class T>
    struct is_fixed_size< _MATRIX_MINUS_SCALAR_<MAO, T> >
    {
        static constexpr bool    value = is_fixed_size<MAO>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAO>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<MAO>::COLS*value;
        typedef typename promote<typename type_traits<MAO>::dtype,T>::type                                    dtype;
        typedef typename is_fixed_size_selectors::default_matrix_selector<value,ROWS,COLS,dtype>::matrix_type matrix_type;
    };

    // Scalar-Matrix Difference class
    template<class T, typename MAO>
    struct is_fixed_size< _SCALAR_MINUS_MATRIX_<T, MAO> >
    {
        static constexpr bool    value = is_fixed_size<MAO>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAO>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<MAO>::COLS*value;
        typedef typename promote<T,typename type_traits<MAO>::dtype>::type                                    dtype;
        typedef typename is_fixed_size_selectors::default_matrix_selector<value,ROWS,COLS,dtype>::matrix_type matrix_type;
    };

    // Matrix-Matrix Product class
    template<typename MAO1,typename MAO2>
    struct is_fixed_size<_MATRIX_TIMES_MATRIX_<MAO1, MAO2> >
    {
        static constexpr bool    value = is_fixed_size<MAO1>::value && is_fixed_size<MAO2>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAO1>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<MAO2>::COLS*value;
        typedef promote<typename type_traits<MAO1>::dtype,
                        typename type_traits<MAO1>::dtype>::type                                              dtype;
        typedef typename is_fixed_size_selectors::default_matrix_selector<value,ROWS,COLS,dtype>::matrix_type matrix_type;
    };

    // Matrix-Scalar Product class
    template<typename MAO, class T>
    struct is_fixed_size< _MATRIX_TIMES_SCALAR_<MAO, T> >
    {
        static constexpr bool    value = is_fixed_size<MAO>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAO>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<MAO>::COLS*value;
        typedef typename promote<typename type_traits<MAO>::dtype,T>::type                                    dtype;
        typedef typename is_fixed_size_selectors::default_matrix_selector<value,ROWS,COLS,dtype>::matrix_type matrix_type;
    };

    // Matrix-Scalar Division class
    template<typename MAO, class T>
    struct is_fixed_size< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >
    {
        static constexpr bool    value = is_fixed_size<MAO>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAO>::ROWS*value;
        static constexpr index_t COLS  = is_fixed_size<MAO>::COLS*value;
        typedef typename promote<typename type_traits<MAO>::dtype,T>::type                                    dtype;
        typedef typename is_fixed_size_selectors::default_matrix_selector<value,ROWS,COLS,dtype>::matrix_type matrix_type;
    };

    // Matrix Transpose class
    template<typename MAO>
    struct is_fixed_size< _MATRIX_TRANSPOSE_<MAO> >
    {
        static constexpr bool    value = is_fixed_size<MAO>::value;
        static constexpr index_t ROWS  = is_fixed_size<MAO>::COLS*value;
        static constexpr index_t COLS  = is_fixed_size<MAO>::ROWS*value;
        typedef typename type_traits<MAO>::dtype         dtype;
        typedef typename is_fixed_size_selectors::default_matrix_selector<value,ROWS,COLS,dtype>::matrix_type matrix_type;
    };

}

#endif