#ifndef _MATL_EXPRESSION_IMPLEMENTATION__HPP
#define _MATL_EXPRESSION_IMPLEMENTATION__HPP
namespace matl
{
   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const&) 
    * 
    * Check if the expression is a safe target for evaluation. Downcasts to derived type.
    * 
    * Parameters:
    *     Matrix target     Target matrix
    * Returns:
    *     'true' if direct evaluation corrupts memory 'false' otherwise
    ***********************************************************************************************/
    template<typename EXPR>
    template<typename MAT>
    constexpr bool _EXPRESSION_<EXPR>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {return static_cast<EXPR const&>(*this).unsafe_evaluation_target(target);}

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const*)
    * 
    * Check if the expression references the address of the target matrix. 
    * Downcasts to derived type.
    * 
    * Parameters:
    *     Matrix target     Target matrix
    * Returns:
    *     'true' if the recursive expression structure references the address of the target matrix
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename EXPR>
    template<typename MAT>
    constexpr bool _EXPRESSION_<EXPR>::references_address(_MATRIX_<MAT> const* target) const
    {return static_cast<EXPR const&>(*this).references_address(target);} 

   /***********************************************************************************************
    * void eval(_MATRIX_& target) const
    * 
    * Evaluate expression to target matrix. Downcasts to derived type.
    * 
    * Template Parameters:
    *     EXPR           Expression type
    *     MAT            Matrix type
    * 
    * Parameters:
    *     Matrix target     Target matrix
    * Returns:  
    ***********************************************************************************************/
    template<typename EXPR>
    template<typename MAT> 
    void _EXPRESSION_<EXPR>::eval(_MATRIX_<MAT>& target) const
    {static_cast<EXPR const&>(*this).eval(target);}





   /***********************************************************************************************
    * index_t rows() const
    * 
    * Number of rows in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR           Expression type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename EXPR>
    index_t _EXPRESSION_<EXPR>::rows() const
    {return static_cast<EXPR const&>(*this).rows();}





   /***********************************************************************************************
    * index_t cols() const
    * 
    * Number of columns in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR           Expression type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename EXPR>
    index_t _EXPRESSION_<EXPR>::cols() const
    {return static_cast<EXPR const&>(*this).cols();}




   /***********************************************************************************************
    * index_t size() const
    * 
    * Number of elements in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR           Expression type
    * 
    * Parameters:
    ***********************************************************************************************/  
    template<typename EXPR>
    index_t _EXPRESSION_<EXPR>::size() const
    {return static_cast<EXPR const&>(*this).size();}




   /***********************************************************************************************
    * bool memory_conflict(Matrix const&) const
    * 
    * Check if there is memory conflict between expression and target matrix
    * 
    * Template Parameters:
    *     EXPR           Expression type
    * 
    * Parameters:
    *     target           Target matrix
    * 
    * Returns:
    *     'true' if direct evaluation corrupts memory 'false' otherwise
    ***********************************************************************************************/  
    /*
    template<typename EXPR>
    template<typename MAT>
    bool _EXPRESSION_<EXPR>::memory_conflict(Matrix<MAT> const& target) const
    {return static_cast<EXPR const&>(*this).memory_conflict(target);}
    */

   /***********************************************************************************************
    * constexpr bool elementwise()
    * 
    * Check if the expression can be evaluated matrix entry-by-entry. 
    * 
    * Returns:
    *     'true' if altering an entry i,j in a constituent RHS matrix affects only entry i,j
    *     of the resultant LHS matrix.
    ***********************************************************************************************/  
    /*
    template<typename EXPR>
    constexpr bool _EXPRESSION_<EXPR>::elementwise()
    {return static_cast<EXPR const&>(*this).elementwise();}
    */
}
#endif
