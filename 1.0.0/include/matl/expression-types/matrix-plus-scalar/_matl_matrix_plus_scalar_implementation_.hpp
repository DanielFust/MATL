/**************************************************************************************************
 * File: _matl_matrix_plus_scalar_implementation_.hpp
 * 
 * Description: Implementations for the matrix plus scalar expression.
 **************************************************************************************************/
#ifndef _MATL_MATRIX_PLUS_SCALAR_IMPLEMENTATION__HPP
#define _MATL_MATRIX_PLUS_SCALAR_IMPLEMENTATION__HPP
namespace matl
{
   /***********************************************************************************************
    * Constructor
    * 
    * Stores the operands: an expression is copied, a matrix is referenced.
    * 
    * Parameters:
    *     mat          Matrix or expression
    *     scalar       Scalar
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename OP>
    _MATRIX_PLUS_SCALAR_<MAO,T>::_MATRIX_PLUS_SCALAR_(OP const& mat, T scalar) : 
        mat_(unqualified_cast(mat)),
        scalar_(scalar)
    {}

   /***********************************************************************************************
    * operand_type const& operand() const
    * 
    * Returns:
    *     Matrix or expression operand
    ***********************************************************************************************/
    template<typename MAO, class T>
    typename _MATRIX_PLUS_SCALAR_<MAO,T>::operand_type const& _MATRIX_PLUS_SCALAR_<MAO,T>::operand() const
    {return this->mat_;}

   /***********************************************************************************************
    * dtype scalar() const
    * 
    * Returns:
    *     Scalar operand
    ***********************************************************************************************/
    template<typename MAO, class T>
    typename _MATRIX_PLUS_SCALAR_<MAO,T>::dtype _MATRIX_PLUS_SCALAR_<MAO,T>::scalar() const
    {return this->scalar_;}

   /***********************************************************************************************
    * index_t rows() const
    * 
    * Returns:
    *      index_t      Number of rows in resultant matrix
    ***********************************************************************************************/ 
    template<typename MAO, class T>
    index_t _MATRIX_PLUS_SCALAR_<MAO,T>::rows() const
    {return this->mat_.rows();}

   /***********************************************************************************************
    * index_t cols() const
    *
    * Returns: 
    *     index_t      Number of columns in resultant matrix
    ***********************************************************************************************/
    template<typename MAO, class T>
    index_t _MATRIX_PLUS_SCALAR_<MAO,T>::cols() const
    {return this->mat_.cols();}

   /***********************************************************************************************
    * index_t size() const
    *
    * Returns: 
    *     index_t      Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename MAO, class T>
    index_t _MATRIX_PLUS_SCALAR_<MAO,T>::size() const
    {return this->mat_.size();}

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const&)
    * 
    * Check if evaluation of the expression directly to the target matrix would corrupt memory.
    * Elementwise evaluation is unsafe if the target is read through a transpose (element (i,j)
    * of the result would read an element of the target that has already been overwritten).
    * 
    * Parameters:
    *     target       Target matrix
    * Returns:
    *     'true' if evaluation directly to the target matrix would corrupt memory
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    constexpr bool _MATRIX_PLUS_SCALAR_<MAO,T>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {
        return (matl::references_address(&target, this->mat_) && transposed_matrix(this->mat_))
               || matl::unsafe_evaluation_target(target, this->mat_);
    }

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const*)
    * 
    * Check if the expression references the target matrix.
    * 
    * Parameters:
    *     target       Target matrix
    * Returns:
    *     'true' if the expression references the target matrix, 'false' otherwise
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    constexpr bool _MATRIX_PLUS_SCALAR_<MAO,T>::references_address(_MATRIX_<MAT> const* target) const
    {return matl::references_address(target, this->mat_);}

   /***********************************************************************************************
    * void eval(_MATRIX_<MAT>& target) const
    * 
    * Evaluate the expression to the target matrix.
    * 
    * Parameters:
    *     target       Target matrix (must be a safe evaluation target with matching dimensions)
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    void _MATRIX_PLUS_SCALAR_<MAO,T>::eval(_MATRIX_<MAT>& target) const
    {assign_elementwise(target, *this);}

   /***********************************************************************************************
    * struct evaluator<_MATRIX_PLUS_SCALAR_>
    * 
    * Evaluates elements of the expression from the evaluators of its operands.
    ***********************************************************************************************/
    template<typename MAO, class T>
    struct evaluator< _MATRIX_PLUS_SCALAR_<MAO,T> >
    {
        typedef _MATRIX_PLUS_SCALAR_<MAO,T> expression_type;
        typedef typename expression_type::dtype dtype;

        evaluator<typename expression_type::operand_type> operand_;
        dtype scalar_;

        explicit evaluator(expression_type const& expr) : operand_(expr.operand()), scalar_(expr.scalar()) {}

        dtype operator() (index_t i, index_t j) const {return operand_(i,j) + scalar_;}
        dtype operator() (index_t i)            const {return operand_(i) + scalar_;}
    };
}
#endif
