/**************************************************************************************************
 * File: _matl_scalar_minus_matrix_implementation_.hpp
 * 
 * Description: Implementations for the scalar minus matrix expression template.
 **************************************************************************************************/
#ifndef _MATL_SCALAR_MINUS_MATRIX_IMPLEMENTATION__HPP
#define _MATL_SCALAR_MINUS_MATRIX_IMPLEMENTATION__HPP
namespace matl
{
   /***********************************************************************************************
    * Constructor
    * 
    * Stores the operands: an expression is copied, a matrix is referenced.
    * 
    * Parameters:
    *     scalar       Scalar
    *     mat          Matrix or expression
    ***********************************************************************************************/
    template<class T, typename MAO>
    template<typename OP>
    _SCALAR_MINUS_MATRIX_<T,MAO>::_SCALAR_MINUS_MATRIX_(T scalar, OP const& mat) : 
        mat_(unqualified_cast(mat)),
        scalar_(scalar)
    {}

   /***********************************************************************************************
    * operand_type const& operand() const
    * 
    * Returns:
    *     Matrix or expression operand
    ***********************************************************************************************/
    template<class T, typename MAO>
    typename _SCALAR_MINUS_MATRIX_<T,MAO>::operand_type const& _SCALAR_MINUS_MATRIX_<T,MAO>::operand() const
    {return this->mat_;}

   /***********************************************************************************************
    * dtype scalar() const
    * 
    * Returns:
    *     Scalar operand
    ***********************************************************************************************/
    template<class T, typename MAO>
    typename _SCALAR_MINUS_MATRIX_<T,MAO>::dtype _SCALAR_MINUS_MATRIX_<T,MAO>::scalar() const
    {return this->scalar_;}

   /***********************************************************************************************
    * index_t rows() const
    * 
    * Returns:
    *      index_t      Number of rows in resultant matrix
    ***********************************************************************************************/ 
    template<class T, typename MAO>
    index_t _SCALAR_MINUS_MATRIX_<T,MAO>::rows() const
    {return this->mat_.rows();}

   /***********************************************************************************************
    * index_t cols() const
    *
    * Returns: 
    *     index_t      Number of columns in resultant matrix
    ***********************************************************************************************/
    template<class T, typename MAO>
    index_t _SCALAR_MINUS_MATRIX_<T,MAO>::cols() const
    {return this->mat_.cols();}

   /***********************************************************************************************
    * index_t size() const
    *
    * Returns: 
    *     index_t      Number of elements in resultant matrix
    ***********************************************************************************************/
    template<class T, typename MAO>
    index_t _SCALAR_MINUS_MATRIX_<T,MAO>::size() const
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
    template<class T, typename MAO>
    template<typename MAT>
    constexpr bool _SCALAR_MINUS_MATRIX_<T,MAO>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
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
    template<class T, typename MAO>
    template<typename MAT>
    constexpr bool _SCALAR_MINUS_MATRIX_<T,MAO>::references_address(_MATRIX_<MAT> const* target) const
    {return matl::references_address(target, this->mat_);}

   /***********************************************************************************************
    * void eval(_MATRIX_<MAT>& target) const
    * 
    * Evaluate the expression to the target matrix.
    * 
    * Parameters:
    *     target       Target matrix (must be a safe evaluation target with matching dimensions)
    ***********************************************************************************************/
    template<class T, typename MAO>
    template<typename MAT>
    void _SCALAR_MINUS_MATRIX_<T,MAO>::eval(_MATRIX_<MAT>& target) const
    {assign_elementwise(target, *this);}

   /***********************************************************************************************
    * struct evaluator<_SCALAR_MINUS_MATRIX_>
    * 
    * Evaluates elements of the expression from the evaluators of its operands.
    ***********************************************************************************************/
    template<class T, typename MAO>
    struct evaluator< _SCALAR_MINUS_MATRIX_<T,MAO> >
    {
        typedef _SCALAR_MINUS_MATRIX_<T,MAO> expression_type;
        typedef typename expression_type::dtype dtype;

        evaluator<typename expression_type::operand_type> operand_;
        dtype scalar_;

        explicit evaluator(expression_type const& expr) : operand_(expr.operand()), scalar_(expr.scalar()) {}

        dtype operator() (index_t i, index_t j) const {return scalar_ - operand_(i,j);}
        dtype operator() (index_t i)            const {return scalar_ - operand_(i);}
    };
}
#endif
