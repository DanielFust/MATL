#ifndef _MATL_MATRIX_TRANSPOSE_IMPLEMENTATION__HPP
#define _MATL_MATRIX_TRANSPOSE_IMPLEMENTATION__HPP

namespace matl
{
   /***********************************************************************************************
    * Constructor
    *
    * Stores the operand: an expression is copied, a matrix is referenced.
    *
    * Parameters:
    *     mat          Matrix or expression to transpose
    ***********************************************************************************************/
    template<typename MAO>
    _MATRIX_TRANSPOSE_<MAO>::_MATRIX_TRANSPOSE_(MAO const& mat) : mat_(mat)
    {}

   /***********************************************************************************************
    * operand_type const& operand() const
    *
    * Returns:
    *     Matrix or expression being transposed
    ***********************************************************************************************/
    template<typename MAO>
    typename _MATRIX_TRANSPOSE_<MAO>::operand_type const& _MATRIX_TRANSPOSE_<MAO>::operand() const
    {return this->mat_;}

   /***********************************************************************************************
    * index_t rows() const
    *
    * Returns:
    *     index_t      Number of rows of the transpose (columns of the operand)
    ***********************************************************************************************/
    template<typename MAO>
    index_t _MATRIX_TRANSPOSE_<MAO>::rows() const
    {return mat_.cols();}

   /***********************************************************************************************
    * index_t cols() const
    *
    * Returns:
    *     index_t      Number of columns of the transpose (rows of the operand)
    ***********************************************************************************************/
    template<typename MAO>
    index_t _MATRIX_TRANSPOSE_<MAO>::cols() const
    {return mat_.rows();}

   /***********************************************************************************************
    * index_t size() const
    *
    * Returns:
    *     index_t      Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename MAO>
    index_t _MATRIX_TRANSPOSE_<MAO>::size() const
    {return mat_.size();}

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const&)
    *
    * Check if evaluation of the expression directly to the target matrix would corrupt memory.
    *
    * Parameters:
    *     target       Target matrix
    * Returns:
    *     'true' if evaluation directly to the target matrix would corrupt memory
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO>
    template<typename MAT>
    constexpr bool _MATRIX_TRANSPOSE_<MAO>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {
        // Pure transpose expressions are safe (see eval), however, arithmetic expressions
        // involving the transpose of the target matrix are not.
        return matl::unsafe_evaluation_target(target,this->mat_);
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
    template<typename MAO>
    template<typename MAT>
    constexpr bool _MATRIX_TRANSPOSE_<MAO>::references_address(_MATRIX_<MAT> const* target) const
    {return matl::references_address(target,this->mat_);}

   /***********************************************************************************************
    * void eval(_MATRIX_<MAT>& target) const
    *
    * Evaluate the transpose to the target matrix.
    *
    * Parameters:
    *     target       Target matrix (must have the dimensions of the transpose)
    ***********************************************************************************************/
    template<typename MAO>
    template<typename MAT>
    void _MATRIX_TRANSPOSE_<MAO>::eval(_MATRIX_<MAT>& target) const
    {
        const index_t M = target.rows();
        const index_t N = target.cols();

        #ifdef MATL_BOUNDS_CHECK
        if ((M != this->mat_.cols()) || (N != this->mat_.rows()))
            throw std::out_of_range("Matrix transpose expression: target matrix dimensions do not match");
        #endif

        // Evaluator of the operand. Any products in the operand are computed here, before any
        // element of the target is written.
        const evaluator<operand_type> source(this->mat_);

        // Special case for square matrices: elements (i,j) and (j,i) are evaluated simultaneously.
        //
        // This covers cases in which the target matrix appears in the expression
        // (e.g. A = transpose(A+A)): element (i,j) of an elementwise operand depends only on
        // element (i,j) of its matrices, so reading both elements of the pair before writing
        // either is memory safe.
        if (M==N)
        {
            for (index_t i = 1; i <= M; i++)
            {
                target(i,i) = source(i,i);
                for (index_t j = i+1; j <= N; j++)
                {
                    const dtype source_ij = source(i,j);
                    const dtype source_ji = source(j,i);
                    target(i,j) = source_ji;
                    target(j,i) = source_ij;
                }
            }
        }
        // Transposing flips the storage order, so a target with the opposite storage order of the
        // operand can be filled in storage order
        else if constexpr (is_same_storage_order<MAT, _MATRIX_TRANSPOSE_<MAO> >::value)
        {
            const index_t K = this->size();
            for (index_t i = 1; i <= K; i++)
                target(i) = source(i);
        }
        else
        {
            for (index_t i = 1; i <= M; i++)
                for (index_t j = 1; j <= N; j++)
                    target(i,j) = source(j,i);
        }
    }

   /***********************************************************************************************
    * struct evaluator<_MATRIX_TRANSPOSE_>
    *
    * Evaluates elements of the transpose from the evaluator of its operand. Element (i,j) of the
    * transpose is element (j,i) of the operand. Single-index access passes through unchanged:
    * the transpose has the opposite storage order of its operand, so the i-th element in storage
    * order is the same element.
    ***********************************************************************************************/
    template<typename MAO>
    struct evaluator< _MATRIX_TRANSPOSE_<MAO> >
    {
        typedef _MATRIX_TRANSPOSE_<MAO> expression_type;
        typedef typename expression_type::dtype dtype;

        evaluator<typename expression_type::operand_type> operand_;

        explicit evaluator(expression_type const& expr) : operand_(expr.operand()) {}

        dtype operator() (index_t i, index_t j) const {return operand_(j,i);}
        dtype operator() (index_t i)            const {return operand_(i);}
    };
}
#endif
