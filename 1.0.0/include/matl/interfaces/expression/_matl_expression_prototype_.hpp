#ifndef _MATL_EXPRESSION_PROTOTYPE__HPP
#define _MATL_EXPRESSION_PROTOTYPE__HPP
namespace matl
{
   /**********************************************************************************************
    * class _EXPRESSION_ (Prototype)
    * 
    * Base class for expression templates based around the Matrix static interface in the 'matl'
    * namespace. This enables syntactically concise expressions that support nuanced behavior
    * (such as lazy execution, memory managment, cache blocking, etc.) specific to the expression
    * 
    * Template Parameters:
    *     EXPR         (_EXPRESSION_) Expression type
    **********************************************************************************************/
    template<typename EXPR>
    class _EXPRESSION_
    {
      public:

        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename type_traits<EXPR>::dtype dtype;

        // Check if target matrix is safe for evaluation
        template<typename MAT>
        constexpr bool unsafe_evaluation_target(_MATRIX_<MAT> const&) const;

        // Check if the expression references the address of the target matrix
        template<typename MAT>
        constexpr bool references_address(_MATRIX_<MAT> const*) const;

        // Evaluate Expression to target matrix
        template<typename MAT> 
        void eval(_MATRIX_<MAT>& target) const;

        // Get number of rows of resultant matrix
        index_t rows() const;

        // Get number of columns of resultant matrix
        index_t cols() const;

        // Get number of elements in resultant matrix
        index_t size() const;

    };
}
#endif
