#ifndef MATL_UTILITIES__HPP
#define MATL_UTILITIES__HPP
namespace matl
{

   
   /***********************************************************************************************
    * constexpr bool references_address(_MATRIX_ const*, _MATRIX_ALGEBRA_OBJECT_ const&)
    * 
    * This utility function allows syntactically identical methods of checking if a matrix
    * address is referenced by an expression or is the address of the queried matrix. 
    * 
    * Parameters:
    *     mat            matrix address to check
    *     expr/other     expression or matrix to check against
    * Returns:
    *     'true' if the address of the matrix is referenced in the expression or is the 
    *            address of the other matrix
    *     'false' otherwise
    *************************************************************************************************/
    
    // Expression specialization (maps to _EXPRESSION_::references_address method)
    template<typename MAT, typename EXPR>
    constexpr bool references_address(_MATRIX_<MAT> const* mat, _EXPRESSION_<EXPR> const& expr)
    {return (expr.references_address(mat));}

    // Matrix specialization (bottom of recursion. Checks if addresses are equal i.e. same object)
    template<typename MAT1, typename MAT2>
    constexpr bool references_address(_MATRIX_<MAT1> const* mat, _MATRIX_<MAT2> const& other)
    {return (static_cast<void const*>(mat) == static_cast<void const*>(&other));}


   /***********************************************************************************************
    * constexpr bool unsafe_evaluation_target(_MATRIX_ const&, _MATRIX_ALGEBRA_OBJECT_ const&)
    * 
    * This utility function allows a uniform syntax for checking if the evaluation of an 
    * expression to a target matrix would corrupt memory.
    * 
    * Parameters:
    *     target         target matrix
    *     expr/other     expression or matrix to check against
    * Returns:
    *     'true'  if the evaluation of the expression to target matrix would corrupt memory. 
    *     'false' otherwise or if the 'other' object is a matrix rather than expression.
    *************************************************************************************************/

    // Expression specialization (maps to _EXPRESSION_::unsafe_evaluation_target method)
    template<typename MAT, typename EXPR>
    constexpr bool unsafe_evaluation_target(_MATRIX_<MAT> const& target, _EXPRESSION_<EXPR> const& expr)
    {return (expr.unsafe_evaluation_target(target));}

    // Matrix specialization (bottom of recursion)
    template<typename MAT1, typename MAT2>
    constexpr bool unsafe_evaluation_target(_MATRIX_<MAT1> const& target, _MATRIX_<MAT2> const& other)
    {return false;}


   /***********************************************************************************************
    * constexpr bool transposed_matrix(T const&)
    * 
    * Succinct mapping to 'is_transposed<T>::value' to improve readability for instantiatable 
    * objects with long or complex qualified typenames 
    * 
    * Parameters:
    *     obj           object (likely expression or matrix) to check
    * Returns:
    *     'true' if the expression or matrix is a _MATRIX_TRANSPOSE_ expression
    *     'false' otherwise
    *************************************************************************************************/
    template<typename T>
    constexpr bool transposed_matrix(T const& obj)
    {return is_transposed<T>::value;}



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
#endif
