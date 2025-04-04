

//#include "./HeterogeneousLinkedList/_matl_linked_list_declarations_.hpp"


namespace matl
{
    // Type promotion struct
    template<typename... TYPES> struct promote;

    // TypeTraits Declarations
    template<typename... T> struct TypeTraits;
    
    template<typename MAT_T>                     struct TypeTraits<Matrix<MAT_T> >;
    template<typename T>                         struct TypeTraits<RowMatrix<T> >;
    template<typename EXPR_T>                    struct TypeTraits<_MATEXPRESSION_<EXPR_T> >;
    template<typename EXPR_T1, typename EXPR_T2> struct TypeTraits<_GENERIC_MAT_SUM_<EXPR_T1, EXPR_T2> >;
    template<typename RMAE_T1, typename RMAE_T2> struct TypeTraits<_ROWMAT_SUM_<RMAE_T1, RMAE_T2> >;
    template<typename EXPR_T1, typename EXPR_T2> struct TypeTraits<_GENERIC_MAT_SUB_<EXPR_T1, EXPR_T2> >;
    template<typename RMAE_T1, typename RMAE_T2> struct TypeTraits<_ROWMAT_SUB_<RMAE_T1, RMAE_T2> >;
    template<typename RMAE,    class    T>       struct TypeTraits<_ROWMAT_SCALAR_SUM_<RMAE, T> >;
    template<typename RMAE,    class    T>       struct TypeTraits<_ROWMAT_SCALAR_PROD_<RMAE, T> >;
    template<typename AE1,     typename AE2>     struct TypeTraits<_GENERIC_MAT_PROD_<AE1, AE2> >;
    template<typename AE,      class    T>       struct TypeTraits<_GENERIC_MAT_SCALAR_PROD_<AE, T> >;
    template<typename AE,      class    T>       struct TypeTraits<_GENERIC_MAT_SCALAR_SUM_<AE, T> >;
    template<typename AE,      class    T>       struct TypeTraits<_GENERIC_MAT_SCALAR_SUB_<AE, T> >;
    template<class    T,       typename AE>      struct TypeTraits<_GENERIC_SCALAR_MAT_SUB_<T, AE> >;
    template<typename MAT>                       struct TypeTraits<_GENERIC_MAT_TRANSPOSE_<MAT> >;
    //template<class T,          bool is_scalar>   struct TypeTraits<T>;               
    //template<class T,          bool>   struct TypeTraits<T>;              

}
