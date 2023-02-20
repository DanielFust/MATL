//#ifndef _MATL_UTILS__HPP
//#define _MATL_UTILS__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains utility functions intended for use by classes and subroutines of the BLAS
 * library
 **************************************************************************************************/
//#include "_matl_forward_declr_.hpp"
//#include "_matl_type_traits_.hpp"

namespace matl
{

   /***********************************************************************************************
    *                       Standard Template Library (STL) Macros
    ***********************************************************************************************/ 
    template< bool B, class T = void >
    using enable_if_t = typename std::enable_if<B,T>::type;
    
    // Enable if Scalar
    template<typename T>
    using enable_if_scalar = std::enable_if_t<std::is_scalar_v<T>, bool>;



   /***********************************************************************************************
    * Type promotion struct
    ***********************************************************************************************/  
    
    // Forward Declarations
    template<typename... TYPES> struct promote;

    // End of type promotion recursion
    template<>  
    struct promote<>
    {
        typedef bool type;
    };

    // General purpose type promotion
    template<typename FIRST_TYPE, typename... OTHER_TYPES>
    struct promote<FIRST_TYPE, OTHER_TYPES...> 
    {
        typedef decltype(FIRST_TYPE() + typename promote<OTHER_TYPES...>::type()) type;
    };


   /***********************************************************************************************
    * C-Style Casting
    * 
    * Provides clearer C++ style syntax for C-style type casts. There are pretty much no 
    * protections with C-style casts so BE CAREFUL!!!
    ***********************************************************************************************/ 
    template<class CAST_T, class T>
    CAST_T cstyle_cast(T& obj) {return (CAST_T)obj;}



   /***********************************************************************************************
    * bool address_overlap(_MATEXPRESSION_ const&, Matrix const&)
    * 
    * Checks if there is memory conflict between the target matrix and expression. If 'true', a
    * new Matrix must be allocated before assignment to target.
    * 
    * Template Parameters:
    *     EXPR        Expression type
    *     MAT         Matrix type
    * 
    * Parameters:
    *     expr        Expression template
    *     mat         Target matrix
    * 
    * Returns:
    *     'true' if target matrix could be corrupted during expression evaluation, 'false' 
    *      otherwise
    ***********************************************************************************************/
    template<typename EXPR, typename MAT>
    bool address_overlap(_MATEXPRESSION_<EXPR> const& expr, Matrix<MAT> const& mat)
    {return expr.memory_conflict(mat);}


   /***********************************************************************************************
    * bool address_overlap(_Matrix const&, Matrix const&)
    * 
    * Checks if there is memory conflict between the target matrix and source matrix (if they are 
    * the same object). If 'true', a new Matrix must be allocated before assignment to target.
    * 
    * Template Parameters:
    *     MAT1         Source matrix type
    *     MAT2         Target matrix type
    * 
    * Parameters:
    *     mat_s        Source matrix
    *     mat_t        Target matrix
    * 
    * Returns:
    *     'true' if target matrix could be corrupted during expression evaluation, 'false' 
    *      otherwise
    ***********************************************************************************************/
    template<typename MAT1, typename MAT2>
    bool address_overlap(Matrix<MAT1> const& mat_s, Matrix<MAT2> const& mat_t)
    {return ((void*)&mat_s == (void*)&mat_t);}  




   /***********************************************************************************************
    * typename promote<T1,T2>::type branchless_min(T1 const&, T2 const&)
    * 
    * Returns the minimum of two objects without branching
    * 
    * Template Arguments:
    *     val1     first value
    *     val2     second value
    ***********************************************************************************************/  
    template<class T1, class T2>
    typename promote<T1,T2>::type branchless_min(T1 const& val1, T2 const& val2)
    {
        bool gt = val1>val2;
        return (!gt)*val1 + gt*val2;
    }


   /***********************************************************************************************
    * typename promote<T1,T2>::type branchless_max(T1 const&, T2 const&)
    * 
    * Returns the maximum of two objects without branching
    * 
    * Template Arguments:
    *     val1     first value
    *     val2     second value
    ***********************************************************************************************/  
    template<class T1, class T2>
    typename promote<T1,T2>::type branchless_max(T1 const& val1, T2 const& val2)
    {
        bool lt = val1<val2;
        return (!lt)*val1 + lt*val2;
    }
}  

//#endif