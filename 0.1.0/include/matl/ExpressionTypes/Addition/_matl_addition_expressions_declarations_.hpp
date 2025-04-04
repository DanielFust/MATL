#ifndef _MATL_ADDITION_EXPRESSIONS_DECLARATIONS__HPP
#define _MATL_ADDITION_EXPRESSIONS_DECLARATIONS__HPP

namespace matl
{
    // Generic Matrix Addition
    template<typename AE_T1, typename AE_T2> class _GENERIC_MAT_SUM_;

    // Generic Matrix - Scalar Addition
    template<typename AE_T, class T> class _GENERIC_MAT_SCALAR_SUM_;

    // RowMatrix - Scalar Addition
    //template<typename RMAE, class T> class _ROWMAT_SCALAR_SUM_;

    // Row Major Matrix Addition
    template<typename RME1, typename RME2> class _ROWMAJ_SUM_;
}

#endif