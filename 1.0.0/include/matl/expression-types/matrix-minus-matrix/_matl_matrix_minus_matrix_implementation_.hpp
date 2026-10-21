/**************************************************************************************************
 * File: _matl_matrix_minus_matrix_implementation_.hpp
 * 
 * Description: Implementations for the matrix minus matrix expression template.
 **************************************************************************************************/
#ifndef _MATL_MATRIX_MINUS_MATRIX_IMPLEMENTATION__HPP
#define _MATL_MATRIX_MINUS_MATRIX_IMPLEMENTATION__HPP
namespace matl
{
   /***********************************************************************************************
    * Constructor
    * 
    * Stores the operands: expressions are copied, matrices are referenced.
    * 
    * Parameters:
    *     mat1         Left operand (matrix or expression)
    *     mat2         Right operand (matrix or expression)
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename OP1, typename OP2>
    _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::_MATRIX_MINUS_MATRIX_(OP1 const& mat1, OP2 const& mat2) : 
        mat1_(unqualified_cast(mat1)), 
        mat2_(unqualified_cast(mat2))
    {}

   /***********************************************************************************************
    * lhs_type const& lhs() const
    * 
    * Returns:
    *     Left operand
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    typename _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::lhs_type const& _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::lhs() const
    {return this->mat1_;}

   /***********************************************************************************************
    * rhs_type const& rhs() const
    * 
    * Returns:
    *     Right operand
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    typename _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::rhs_type const& _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::rhs() const
    {return this->mat2_;}

   /***********************************************************************************************
    * index_t rows() const
    * 
    * Returns:
    *      index_t      Number of rows in resultant matrix
    ***********************************************************************************************/ 
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::rows() const
    {return this->mat1_.rows();}

   /***********************************************************************************************
    * index_t cols() const
    *
    * Returns: 
    *     index_t      Number of columns in resultant matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::cols() const
    {return this->mat1_.cols();}

   /***********************************************************************************************
    * index_t size() const
    *
    * Returns: 
    *     index_t      Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::size() const
    {return this->mat1_.size();}

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
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    constexpr bool _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {
        const bool branch1_unsafe = (matl::references_address(&target, this->mat1_) && transposed_matrix(this->mat1_))
                                    || matl::unsafe_evaluation_target(target, this->mat1_);
        const bool branch2_unsafe = (matl::references_address(&target, this->mat2_) && transposed_matrix(this->mat2_))
                                    || matl::unsafe_evaluation_target(target, this->mat2_);
        return branch1_unsafe || branch2_unsafe;
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
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    constexpr bool _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::references_address(_MATRIX_<MAT> const* target) const
    {return matl::references_address(target, this->mat1_) || matl::references_address(target, this->mat2_);}

   /***********************************************************************************************
    * void eval(_MATRIX_<MAT>& target) const
    * 
    * Evaluate the expression to the target matrix.
    * 
    * Parameters:
    *     target       Target matrix (must be a safe evaluation target with matching dimensions)
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    void _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::eval(_MATRIX_<MAT>& target) const
    {assign_elementwise(target, *this);}

   /***********************************************************************************************
    * struct evaluator<_MATRIX_MINUS_MATRIX_>
    * 
    * Evaluates elements of the expression from the evaluators of its operands.
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    struct evaluator< _MATRIX_MINUS_MATRIX_<MAO1,MAO2> >
    {
        typedef _MATRIX_MINUS_MATRIX_<MAO1,MAO2> expression_type;
        typedef typename expression_type::dtype dtype;

        evaluator<typename expression_type::lhs_type> lhs_;
        evaluator<typename expression_type::rhs_type> rhs_;

        explicit evaluator(expression_type const& expr) : lhs_(expr.lhs()), rhs_(expr.rhs()) {}

        dtype operator() (index_t i, index_t j) const {return lhs_(i,j) - rhs_(i,j);}
        dtype operator() (index_t i)            const {return lhs_(i) - rhs_(i);}
    };
}
#endif
