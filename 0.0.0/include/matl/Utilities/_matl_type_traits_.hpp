//#ifndef _MATL_TYPE_TRAITS__HPP
//#define _MATL_TYPE_TRAITS__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains type trait templates and template metaprogramming utilities required to 
 * correctly determine return types and support static polymorphism
 **************************************************************************************************/
//#include "../_matl_forward_declarations_.hpp"

namespace matl
{

   /***********************************************************************************************
    * TypeTraits
    * 
    * Template specializations provide data type information
    ***********************************************************************************************/ 
    

    
    // Matrix base class return types
    template<typename MAT_T>
    struct TypeTraits< Matrix<MAT_T> >
    {
        typedef typename TypeTraits<MAT_T>::dtype                           dtype;
        typedef typename TypeTraits<MAT_T>::matrix_type                     matrix_type;
        typedef          Matrix<typename TypeTraits<MAT_T>::downcast_type>  resolved_type;
        typedef typename TypeTraits<MAT_T>::downcast_type                   downcast_type;
    };

    // RowMatrix return types
    template<typename T>
    struct TypeTraits<RowMatrix<T> >
    {
        typedef T                        dtype;
        typedef RowMatrix<T>             matrix_type;
        typedef Matrix< RowMatrix<T> >   resolved_type;
        typedef RowMatrix<T>             downcast_type;
    };


    // Matrix Expression Base class return types
    template<typename EXPR_T>
    struct TypeTraits< _MATEXPRESSION_<EXPR_T> >
    {
        typedef typename TypeTraits<EXPR_T>::dtype                            dtype;
        typedef typename TypeTraits<EXPR_T>::matrix_type                      matrix_type;
        typedef  _MATEXPRESSION_<typename TypeTraits<EXPR_T>::downcast_type > resolved_type;
        typedef typename TypeTraits<EXPR_T>::downcast_type                    downcast_type;
    };



    // General Matrix Summation return types
    template<typename EXPR_T1, typename EXPR_T2>
    struct TypeTraits<_GENERIC_MAT_SUM_<EXPR_T1, EXPR_T2> >
    {
        typedef typename promote<typename TypeTraits<EXPR_T1>::dtype,
                                 typename TypeTraits<EXPR_T2>::dtype>::type              dtype;
        typedef typename TypeTraits<EXPR_T1>::matrix_type                                matrix_type;     // Chooses first expression matrix type by default
        typedef          _GENERIC_MAT_SUM_<typename TypeTraits<EXPR_T1>::resolved_type,
                                           typename TypeTraits<EXPR_T2>::resolved_type>  resolved_type;
        typedef          resolved_type                                                   downcast_type;
    };


    // RowMatrix-RowMatrix summation return type
    template<typename RMAE_T1, typename RMAE_T2>
    struct TypeTraits<_ROWMAT_SUM_<RMAE_T1, RMAE_T2> >
    {
        typedef typename promote<typename TypeTraits<RMAE_T1>::dtype,
                                 typename TypeTraits<RMAE_T2>::dtype>::type          dtype;
        typedef          RowMatrix<dtype>                                            matrix_type;
        typedef          _ROWMAT_SUM_<typename TypeTraits<RMAE_T1>::resolved_type,
                                      typename TypeTraits<RMAE_T2>::resolved_type>   resolved_type;
        typedef          resolved_type                                               downcast_type;
    };


    // General Matrix Subtraction return types
    template<typename EXPR_T1, typename EXPR_T2>
    struct TypeTraits<_GENERIC_MAT_SUB_<EXPR_T1, EXPR_T2> >
    {
        typedef typename promote<typename TypeTraits<EXPR_T1>::dtype,
                                 typename TypeTraits<EXPR_T2>::dtype>::type              dtype;
        typedef typename TypeTraits<EXPR_T1>::matrix_type                                matrix_type;    // Defaults to first expression matrix type
        typedef          _GENERIC_MAT_SUB_<typename TypeTraits<EXPR_T1>::resolved_type,
                                           typename TypeTraits<EXPR_T2>::resolved_type>  resolved_type;
        typedef          resolved_type                                                   downcast_type;
    };


    // RowMatrix-RowMatrix subtraction return type
    template<typename RMAE_T1, typename RMAE_T2>
    struct TypeTraits<_ROWMAT_SUB_<RMAE_T1, RMAE_T2> >
    {
        typedef typename promote<typename TypeTraits<RMAE_T1>::dtype,
                                 typename TypeTraits<RMAE_T2>::dtype>::type          dtype;
        typedef          RowMatrix<dtype>                                            matrix_type;
        typedef          _ROWMAT_SUB_<typename TypeTraits<RMAE_T1>::resolved_type,
                                      typename TypeTraits<RMAE_T2>::resolved_type>   resolved_type;
        typedef          resolved_type                                               downcast_type;
    };



    // RowMatrix Scalar Sum
    template<typename RMAE, class T>
    struct TypeTraits<_ROWMAT_SCALAR_SUM_<RMAE, T> >
    {
        typedef typename promote<typename TypeTraits<RMAE>::dtype, T>::type                dtype;
        typedef          RowMatrix<dtype>                                                  matrix_type;
        typedef          _ROWMAT_SCALAR_SUM_<typename TypeTraits<RMAE>::resolved_type,T>   resolved_type;
        typedef          resolved_type                                                     downcast_type;
    };


    // Primative Data type
    //template<typename MAT_T, class T, enable_if_scalar<T> = true>
    //operator+ (Matrix<MAT_T> const& mat, T val) 
    /*
    template<class T, enable_if_scalar<T> = true>
    */
    template<class T>
    struct TypeTraits<T>
    {
        //typedef T             dtype;
        typedef float             dtype;
        typedef RowMatrix<dtype>  matrix_type; // Shouldn't come up (but chooses RowMatrix by default I suppose)
        typedef T             resolved_type;
        typedef T             downcast_type;
    }; 
    

    // Generic Matrix Tranpose
    template<typename MAT>
    struct TypeTraits<_GENERIC_MAT_TRANSPOSE_<MAT> >
    {
        typedef typename TypeTraits<MAT>::dtype                                             dtype;
        typedef typename TypeTraits<MAT>::matrix_type                                       matrix_type; 
        typedef          _GENERIC_MAT_TRANSPOSE_<typename TypeTraits<MAT>::resolved_type>   resolved_type;
        typedef          resolved_type                                                      downcast_type;
    };
   
    
    // Generic Matrix - Matrix Multiplication
    template<typename AE1, typename AE2>
    struct TypeTraits<_GENERIC_MAT_PROD_<AE1, AE2> >
    {
        typedef typename promote<typename TypeTraits<AE1>::dtype,
                                 typename TypeTraits<AE2>::dtype>::type               dtype;         // Scalar does not change data type
        typedef typename TypeTraits<AE1>::matrix_type                                 matrix_type; 
        typedef          _GENERIC_MAT_PROD_<typename TypeTraits<AE1>::resolved_type,
                                            typename TypeTraits<AE2>::resolved_type>  resolved_type;
        typedef          resolved_type                                                downcast_type;
    };

    // RowMatrix Scalar Product return type
    template<typename RMAE, class T>
    struct TypeTraits<_ROWMAT_SCALAR_PROD_<RMAE, T> >
    {
        typedef typename TypeTraits<RMAE>::dtype                                            dtype;           // Scalar does not change data type
        typedef          RowMatrix<dtype>                                                   matrix_type;
        typedef          _ROWMAT_SCALAR_PROD_<typename TypeTraits<RMAE>::resolved_type,T>   resolved_type;
        typedef          resolved_type                                                      downcast_type;
    };



    template<typename AE, class T>
    struct TypeTraits<_GENERIC_MAT_SCALAR_PROD_<AE, T> >
    {
        typedef typename TypeTraits<AE>::dtype                                                dtype;         // Scalar does not change data type
        typedef typename TypeTraits<AE>::matrix_type                                          matrix_type; 
        typedef          _GENERIC_MAT_SCALAR_PROD_<typename TypeTraits<AE>::resolved_type,T>  resolved_type;
        typedef          resolved_type                                                        downcast_type;
    };



    // Generic Matrix - Scalar Addition
    template<typename AE, class T>
    struct TypeTraits<_GENERIC_MAT_SCALAR_SUM_<AE, T> >
    {
        typedef typename TypeTraits<AE>::dtype                                               dtype;         // Scalar does not change data type
        typedef typename TypeTraits<AE>::matrix_type                                         matrix_type; 
        typedef          _GENERIC_MAT_SCALAR_SUM_<typename TypeTraits<AE>::resolved_type,T>  resolved_type;
        typedef          resolved_type                                                       downcast_type;
    };



    // Generic Matrix - Scalar Subtraction
    template<typename AE, class T>
    struct TypeTraits<_GENERIC_MAT_SCALAR_SUB_<AE, T> >
    {
        typedef typename TypeTraits<AE>::dtype                                               dtype;         // Scalar does not change data type
        typedef typename TypeTraits<AE>::matrix_type                                         matrix_type; 
        typedef          _GENERIC_MAT_SCALAR_SUB_<typename TypeTraits<AE>::resolved_type,T>  resolved_type;
        typedef          resolved_type                                                       downcast_type;
    };



    // Generic Scalar - Matrix Subtraction
    template<class T, typename AE>
    struct TypeTraits<_GENERIC_SCALAR_MAT_SUB_<T, AE> >
    {
        typedef typename TypeTraits<AE>::dtype                                               dtype;         // Scalar does not change data type
        typedef typename TypeTraits<AE>::matrix_type                                         matrix_type; 
        typedef          _GENERIC_SCALAR_MAT_SUB_<typename TypeTraits<AE>::resolved_type,T>  resolved_type;
        typedef          resolved_type                                                       downcast_type;
    };





    // RowMatrix-RowMatrix multiplication return type
    template<typename T1, typename T2>
    struct TypeTraits<_ROWMAT_PROD_<RowMatrix<T1>, RowMatrix<T2> > >
    {
        typedef typename promote<T1, T2>::type                        dtype;
        typedef          RowMatrix<dtype>                             matrix_type;
        typedef          _ROWMAT_PROD_<RowMatrix<T1>, RowMatrix<T2> > resolved_type;
        typedef          resolved_type                                downcast_type;
    };



   /***************************************************************************
    *  downcast
    * 
    * This constexpr statically down-casts expression and matrix types to their
    * lowest applicable subclass
    ***************************************************************************/
    //template<typename... TYPES> downcast



    template<typename TYPE>
    typename TypeTraits<TYPE>::downcast_type const& downcast(TYPE const&  obj)
    {
        return static_cast<typename TypeTraits<TYPE>::downcast_type const&>(obj);
    };

    template<typename TYPE>
    typename TypeTraits<TYPE>::downcast_type& downcast(TYPE&  obj)
    {
        return static_cast<typename TypeTraits<TYPE>::downcast_type&>(obj);
    };

    
    
}


//#endif