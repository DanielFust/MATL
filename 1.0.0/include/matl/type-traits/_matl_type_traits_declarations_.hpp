#ifndef _MATL_TYPE_TRAITS_DECLARATIONS__HPP
#define _MATL_TYPE_TRAITS_DECLARATIONS__HPP

namespace matl
{
    // type_traits Declarations
    template<typename... T>    struct type_traits;
    // Interfaces
    template<typename MAO>     struct type_traits< _MATRIX_ALGEBRA_OBJECT_<MAO> >;
    template<typename EW>      struct type_traits< _ELEMENTWISE_<EW> >;
    template<typename MAT>     struct type_traits< _MATRIX_<MAT> >;
    template<typename EXPR>    struct type_traits< _EXPRESSION_<EXPR> >;
    template<typename FS>      struct type_traits< _FIXED_SIZE_<FS> >;
    template<typename FSM>     struct type_traits< _FIXED_SIZE_MATRIX_<FSM> >;
    // Matrix Types
    template<class T>                        struct type_traits< RowMajorMatrix<T> >;
    template<class T, index_t M, index_t N>  struct type_traits< FixedSizeRowMajorMatrix<T, M, N> >;
    // Expression Types
    template<typename MAO1,typename MAO2>    struct type_traits< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>          struct type_traits< _MATRIX_PLUS_SCALAR_<MAO, T> >;
    template<typename MAO1,typename MAO2>    struct type_traits< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>          struct type_traits< _MATRIX_MINUS_SCALAR_<MAO, T> >;
    template<class T, typename MAO>          struct type_traits< _SCALAR_MINUS_MATRIX_<T, MAO> >;
    template<typename MAO, class T>          struct type_traits< _MATRIX_TIMES_SCALAR_<MAO, T> >;
    template<typename MAT1, typename MAT2>   struct type_traits< _MATRIX_TIMES_MATRIX_<MAT1, MAT2> >;
    template<typename MAO, class T>          struct type_traits< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >;
    template<typename MAO>                   struct type_traits< _MATRIX_TRANSPOSE_<MAO> >;
    template<typename MAT>                   struct type_traits< _OWNED_MATRIX_<MAT> >;

    // Matrix Algebra Object Declarations
    template<typename... T> struct is_matrix_algebra_object;
    // Default for no specialization
    template<typename T>  struct is_matrix_algebra_object<T>;
    // Interfaces
    template<typename MAO>  struct is_matrix_algebra_object< _MATRIX_ALGEBRA_OBJECT_<MAO> >;
    template<typename EW>   struct is_matrix_algebra_object< _ELEMENTWISE_<EW> >;
    template<typename MAT>  struct is_matrix_algebra_object< _MATRIX_<MAT> >;
    template<typename EXPR> struct is_matrix_algebra_object< _EXPRESSION_<EXPR> >;
    template<typename FS>   struct is_matrix_algebra_object< _FIXED_SIZE_<FS> >;
    template<typename FSM>  struct is_matrix_algebra_object< _FIXED_SIZE_MATRIX_<FSM> >;
    // Matrix Types
    template<class T>                       struct is_matrix_algebra_object< RowMajorMatrix<T> >;
    template<class T, index_t M, index_t N> struct is_matrix_algebra_object< FixedSizeRowMajorMatrix<T, M, N> >;
    // Expression Types
    template<typename MAO1,typename MAO2>  struct is_matrix_algebra_object< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct is_matrix_algebra_object< _MATRIX_PLUS_SCALAR_<MAO, T> >;
    template<typename MAO1,typename MAO2>  struct is_matrix_algebra_object< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct is_matrix_algebra_object< _MATRIX_MINUS_SCALAR_<MAO, T> >;
    template<class T, typename MAO>        struct is_matrix_algebra_object< _SCALAR_MINUS_MATRIX_<T, MAO> >;
    template<typename MAO, class T>        struct is_matrix_algebra_object< _MATRIX_TIMES_SCALAR_<MAO, T> >;
    template<typename MAT1, typename MAT2> struct is_matrix_algebra_object< _MATRIX_TIMES_MATRIX_<MAT1, MAT2> >;
    template<typename MAO, class T>        struct is_matrix_algebra_object< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >;
    template<typename MAO>                 struct is_matrix_algebra_object< _MATRIX_TRANSPOSE_<MAO> >;
    template<typename MAT>                 struct is_matrix_algebra_object< _OWNED_MATRIX_<MAT> >;

    // Elementwise Declarations
    template<typename... T> struct is_elementwise;
    // Interfaces
    template<typename MAO>  struct is_elementwise< _MATRIX_ALGEBRA_OBJECT_<MAO> >;
    template<typename EW>   struct is_elementwise< _ELEMENTWISE_<EW> >;
    template<typename MAT>  struct is_elementwise< _MATRIX_<MAT> >;
    template<typename EXPR> struct is_elementwise< _EXPRESSION_<EXPR> >;
    template<typename FS>   struct is_elementwise< _FIXED_SIZE_<FS> >;
    template<typename FSM>  struct is_elementwise< _FIXED_SIZE_MATRIX_<FSM> >;
    // Matrix Types
    template<class T>                       struct is_elementwise< RowMajorMatrix<T> >;
    template<class T, index_t M, index_t N> struct is_elementwise< FixedSizeRowMajorMatrix<T, M, N> >;
    // Expression Types
    template<typename MAO1,typename MAO2>  struct is_elementwise< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct is_elementwise< _MATRIX_PLUS_SCALAR_<MAO, T> >;
    template<typename MAO1,typename MAO2>  struct is_elementwise< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct is_elementwise< _MATRIX_MINUS_SCALAR_<MAO, T> >;
    template<class T, typename MAO>        struct is_elementwise< _SCALAR_MINUS_MATRIX_<T, MAO> >;
    template<typename MAO, class T>        struct is_elementwise< _MATRIX_TIMES_SCALAR_<MAO, T> >;
    template<typename MAT1, typename MAT2> struct is_elementwise< _MATRIX_TIMES_MATRIX_<MAT1, MAT2> >;
    template<typename MAO, class T>        struct is_elementwise< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >;
    template<typename MAO>                 struct is_elementwise< _MATRIX_TRANSPOSE_<MAO> >;
    template<typename MAT>                 struct is_elementwise< _OWNED_MATRIX_<MAT> >;


    // FixedSize Declarations
    template<typename... T> struct is_fixed_size;
    // Default for no specialization
    template<typename T>  struct is_fixed_size<T>;
    // Interfaces
    template<typename MAO>  struct is_fixed_size< _MATRIX_ALGEBRA_OBJECT_<MAO> >;
    template<typename EW>   struct is_fixed_size< _ELEMENTWISE_<EW> >;
    template<typename MAT>  struct is_fixed_size< _MATRIX_<MAT> >;
    template<typename EXPR> struct is_fixed_size< _EXPRESSION_<EXPR> >;
    template<typename FS>   struct is_fixed_size< _FIXED_SIZE_<FS> >;
    template<typename FSM>  struct is_fixed_size< _FIXED_SIZE_MATRIX_<FSM> >;
    // Matrix Types
    template<class T>                       struct is_fixed_size< RowMajorMatrix<T> >;
    template<class T, index_t M, index_t N> struct is_fixed_size< FixedSizeRowMajorMatrix<T, M, N> >;
    // Expression Types
    template<typename MAO1,typename MAO2>  struct is_fixed_size< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct is_fixed_size< _MATRIX_PLUS_SCALAR_<MAO, T> >;
    template<typename MAO1,typename MAO2>  struct is_fixed_size< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct is_fixed_size< _MATRIX_MINUS_SCALAR_<MAO, T> >;
    template<class T, typename MAO>        struct is_fixed_size< _SCALAR_MINUS_MATRIX_<T, MAO> >;
    template<typename MAO, class T>        struct is_fixed_size< _MATRIX_TIMES_SCALAR_<MAO, T> >;
    template<typename MAT1, typename MAT2> struct is_fixed_size< _MATRIX_TIMES_MATRIX_<MAT1, MAT2> >;
    template<typename MAO, class T>        struct is_fixed_size< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >;
    template<typename MAO>                 struct is_fixed_size< _MATRIX_TRANSPOSE_<MAO> >;
    template<typename MAT>                 struct is_fixed_size< _OWNED_MATRIX_<MAT> >;

    // is_transposed Declarations
    template<typename... T> struct is_transposed;
    // Default (false)
    template<typename T>  struct is_transposed<T>;
    // Interfaces (downcast to derived type)
    template<typename MAO>  struct is_transposed< _MATRIX_ALGEBRA_OBJECT_<MAO> >;
    template<typename EW>   struct is_transposed< _ELEMENTWISE_<EW> >;
    template<typename EXPR> struct is_transposed< _EXPRESSION_<EXPR> >;
    template<typename FS>   struct is_transposed< _FIXED_SIZE_<FS> >;
    // Only _MATRIX_TRANSPOSE_ is transposed
    template<typename MAO>  struct is_transposed< _MATRIX_TRANSPOSE_<MAO> >;
    template<typename MAT>  struct is_transposed< _OWNED_MATRIX_<MAT> >;
    

    // is_rowmajor_storage Declarations
    template<typename... T> struct is_rowmajor_storage;
    // Recusion
    template<typename FIRST_TYPE, typename ... OTHER_TYPES> struct is_rowmajor_storage<FIRST_TYPE, OTHER_TYPES...>;
    // Default (false)
    template<typename T>  struct is_rowmajor_storage<T>;
    // Interfaces (downcast to derived type)
    template<typename MAO>  struct is_rowmajor_storage< _MATRIX_ALGEBRA_OBJECT_<MAO> >;
    template<typename EW>   struct is_rowmajor_storage< _ELEMENTWISE_<EW> >;
    template<typename MAT>  struct is_rowmajor_storage< _MATRIX_<MAT> >;
    template<typename EXPR> struct is_rowmajor_storage< _EXPRESSION_<EXPR> >;
    template<typename FS>   struct is_rowmajor_storage< _FIXED_SIZE_<FS> >;
    template<typename FSM>  struct is_rowmajor_storage< _FIXED_SIZE_MATRIX_<FSM> >;
    // Matrix Types
    template<class T>                       struct is_rowmajor_storage< RowMajorMatrix<T> >;
    template<class T, index_t M, index_t N> struct is_rowmajor_storage< FixedSizeRowMajorMatrix<T, M, N> >;
    // Expression Types
    template<typename MAO1,typename MAO2>  struct is_rowmajor_storage< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct is_rowmajor_storage< _MATRIX_PLUS_SCALAR_<MAO, T> >;
    template<typename MAO1,typename MAO2>  struct is_rowmajor_storage< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct is_rowmajor_storage< _MATRIX_MINUS_SCALAR_<MAO, T> >;
    template<class T, typename MAO>        struct is_rowmajor_storage< _SCALAR_MINUS_MATRIX_<T, MAO> >;
    template<typename MAO, class T>        struct is_rowmajor_storage< _MATRIX_TIMES_SCALAR_<MAO, T> >;
    template<typename MAT1, typename MAT2> struct is_rowmajor_storage< _MATRIX_TIMES_MATRIX_<MAT1, MAT2> >;
    template<typename MAO, class T>        struct is_rowmajor_storage< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >;
    template<typename MAO>                 struct is_rowmajor_storage< _MATRIX_TRANSPOSE_<MAO> >;
    template<typename MAT>                 struct is_rowmajor_storage< _OWNED_MATRIX_<MAT> >;


    // is_colmajor_storage Declarations
    template<typename... T> struct is_colmajor_storage;
    // Recusion
    template<typename FIRST_TYPE, typename ... OTHER_TYPES> struct is_colmajor_storage<FIRST_TYPE, OTHER_TYPES...>;
    // Default (false)
    template<typename T>  struct is_colmajor_storage<T>;
    // Interfaces (downcast to derived type)
    template<typename MAO>  struct is_colmajor_storage< _MATRIX_ALGEBRA_OBJECT_<MAO> >;
    template<typename EW>   struct is_colmajor_storage< _ELEMENTWISE_<EW> >;
    template<typename MAT>  struct is_colmajor_storage< _MATRIX_<MAT> >;
    template<typename EXPR> struct is_colmajor_storage< _EXPRESSION_<EXPR> >;
    template<typename FS>   struct is_colmajor_storage< _FIXED_SIZE_<FS> >;
    template<typename FSM>  struct is_colmajor_storage< _FIXED_SIZE_MATRIX_<FSM> >;
    // Matrix Types
    template<class T>                       struct is_colmajor_storage< RowMajorMatrix<T> >;
    template<class T, index_t M, index_t N> struct is_colmajor_storage< FixedSizeRowMajorMatrix<T, M, N> >;
    // Expression Types
    template<typename MAO1,typename MAO2>  struct is_colmajor_storage< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct is_colmajor_storage< _MATRIX_PLUS_SCALAR_<MAO, T> >;
    template<typename MAO1,typename MAO2>  struct is_colmajor_storage< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct is_colmajor_storage< _MATRIX_MINUS_SCALAR_<MAO, T> >;
    template<class T, typename MAO>        struct is_colmajor_storage< _SCALAR_MINUS_MATRIX_<T, MAO> >;
    template<typename MAO, class T>        struct is_colmajor_storage< _MATRIX_TIMES_SCALAR_<MAO, T> >;
    template<typename MAT1, typename MAT2> struct is_colmajor_storage< _MATRIX_TIMES_MATRIX_<MAT1, MAT2> >;
    template<typename MAO, class T>        struct is_colmajor_storage< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >;
    template<typename MAO>                 struct is_colmajor_storage< _MATRIX_TRANSPOSE_<MAO> >;
    template<typename MAT>                 struct is_colmajor_storage< _OWNED_MATRIX_<MAT> >;

    // is_same_storage_order Declarations
    template<typename T, typename ... OTHER_TYPES> struct is_same_storage_order;
    // End of recursion
    template<typename T> struct is_same_storage_order<T>;

    // Member type declarations
    template<typename T>    struct member_type;
    // Interfaces
    template<typename MAO>  struct member_type< _MATRIX_ALGEBRA_OBJECT_<MAO> >;
    template<typename EW>   struct member_type< _ELEMENTWISE_<EW> >;
    template<typename MAT>  struct member_type< _MATRIX_<MAT> >;
    template<typename EXPR> struct member_type< _EXPRESSION_<EXPR> >;
    template<typename FS>   struct member_type< _FIXED_SIZE_<FS> >;
    template<typename FSM>  struct member_type< _FIXED_SIZE_MATRIX_<FSM> >;
    // Matrix Types
    template<class T>                       struct member_type< RowMajorMatrix<T> >;
    template<class T, index_t M, index_t N> struct member_type< FixedSizeRowMajorMatrix<T, M, N> >;
    // Expression Types
    template<typename MAO1,typename MAO2>  struct member_type< _MATRIX_PLUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct member_type< _MATRIX_PLUS_SCALAR_<MAO, T> >;
    template<typename MAO1,typename MAO2>  struct member_type< _MATRIX_MINUS_MATRIX_<MAO1, MAO2> >;
    template<typename MAO, class T>        struct member_type< _MATRIX_MINUS_SCALAR_<MAO, T> >;
    template<class T, typename MAO>        struct member_type< _SCALAR_MINUS_MATRIX_<T, MAO> >;
    template<typename MAO, class T>        struct member_type< _MATRIX_TIMES_SCALAR_<MAO, T> >;
    template<typename MAT1, typename MAT2> struct member_type< _MATRIX_TIMES_MATRIX_<MAT1, MAT2> >;
    template<typename MAO, class T>        struct member_type< _MATRIX_DIVIDED_BY_SCALAR_<MAO, T> >;
    template<typename MAO>                 struct member_type< _MATRIX_TRANSPOSE_<MAO> >;
    template<typename MAT>                 struct member_type< _OWNED_MATRIX_<MAT> >;
}

#endif