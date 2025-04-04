#ifndef _MATL_FORWARD_DECLR__HPP
#define _MATL_FORWARD_DECLR__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains forward declarations of classes for the MATL library
 **************************************************************************************************/
namespace matl
{
    // Index type (and size) of matl matrices
    typedef unsigned long int index_t;

    // Blocking
    void set_chunk_size(index_t);
    void reset_chunk_size();
    index_t chunk_size();

    // Matrix Types
    template<typename MAT_T> class Matrix;
    template<class T>        class RowMatrix;
    //template<class T>        class ColMatrix;

    // Matrix Expression Types
    template<typename EXPR_T> class _MATEXPRESSION_;
    template<typename EXPR_T> class _ATOMIC_EXPR_;

    // Generic Atomic Expressions
    template<typename AE_T1, typename AE_T2> class _GENERIC_MAT_SUM_;
    template<typename AE_T1, typename AE_T2> class _GENERIC_MAT_SUB_;
    template<typename AE, class T>           class _GENERIC_MAT_SCALAR_PROD_;

    // RowMatrix Atomic Expressions
    template<typename RMAE_T>                class _ROWMAT_ATOMIC_EXPR_;
    template<typename RMAE1, typename RMAE2> class _ROWMAT_SUM_;
    template<typename RMAE1, typename RMAE2> class _ROWMAT_SUB_;
    template<typename RMAE,  class T>        class _ROWMAT_SCALAR_PROD_;
    template<typename RMAE,  class T>        class _ROWMAT_SCALAR_SUM_;
    template<class T>                        class _ROWMAT_ATOMIC_;

    // RowMatrix-RowMatrix Multiplication
    template<class T1, class T2>    class _ROWMAT_PROD_;  



    /********************************* Matrix Operations *****************************************/
    //template<typename AEXPR_T, typename EXPR_T1, typename EXPR_T2>
    //_ATOMIC_EXPR_<AEXPR_T> MATSUM(EXPR_T1 const&, EXPR_T2 const&); 



    /**************************************** Operators ******************************************/
    //template<typename EXPR_T, typename MATEXPR_T1, typename MATEXPR_T2>
    //_ATOMIC_EXPR_<EXPR_T> operator+ (MATEXPR_T1, MATEXPR_T2);

   

}
#endif