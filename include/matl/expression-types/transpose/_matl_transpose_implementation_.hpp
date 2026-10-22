/**************************************************************************************************
 * File: _matl_transpose_implementation_.hpp
 *
 * Description: Implementation of the transpose expression template.
 **************************************************************************************************/
#ifndef MATL_TRANSPOSE_IMPLEMENTATION_HPP
#define MATL_TRANSPOSE_IMPLEMENTATION_HPP

namespace matl::expr
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
    Transpose<MAO>::Transpose(MAO const& mat) : mat_(mat)
    {}

   /***********************************************************************************************
    * MAO const& operand() const
    *
    * Returns:
    *     Matrix or expression being transposed
    ***********************************************************************************************/
    template<typename MAO>
    MAO const& Transpose<MAO>::operand() const
    {return this->mat_;}

   /***********************************************************************************************
    * index_t rows() const
    *
    * Returns:
    *     index_t      Number of rows of the transpose (columns of the operand)
    ***********************************************************************************************/
    template<typename MAO>
    index_t Transpose<MAO>::rows() const
    {return mat_.cols();}

   /***********************************************************************************************
    * index_t cols() const
    *
    * Returns:
    *     index_t      Number of columns of the transpose (rows of the operand)
    ***********************************************************************************************/
    template<typename MAO>
    index_t Transpose<MAO>::cols() const
    {return mat_.rows();}

   /***********************************************************************************************
    * index_t size() const
    *
    * Returns:
    *     index_t      Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename MAO>
    index_t Transpose<MAO>::size() const
    {return mat_.size();}

   /***********************************************************************************************
    * void check_dimensions() const
    *
    * A transpose fits any operand; checks the operand (recursively).
    ***********************************************************************************************/
    template<typename MAO>
    void Transpose<MAO>::check_dimensions() const
    {detail::check_dimensions(this->mat_);}

   /***********************************************************************************************
    * bool unsafe_evaluation_target(MAT const&)
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
    template<WritableMatrix MAT>
    bool Transpose<MAO>::unsafe_evaluation_target(MAT const& target) const
    {
        // Transposing the target itself is safe (A = transpose(A): eval swaps the elements in
        // pairs), but a matrix or view that shares elements with the target at other positions
        // is not (e.g. a view of other rows of the target), and neither are unsafe expressions.
        // (Arithmetic on the transpose of the target is checked by the elementwise expression,
        // see detail::unsafe_elementwise_operand.)
        return detail::unsafe_operand(target, this->mat_);
    }

   /***********************************************************************************************
    * bool reads(detail::Footprint const&)
    *
    * Check if the expression may read elements of the target (see detail::Footprint).
    *
    * Parameters:
    *     target       Footprint of the target (where its elements are in memory)
    * Returns:
    *     'true' if the expression references the target matrix, 'false' otherwise
    ***********************************************************************************************/
    template<typename MAO>
    bool Transpose<MAO>::reads(detail::Footprint const& target) const
    {return detail::references(target, this->mat_);}

   /***********************************************************************************************
    * void eval(MAT& target) const
    *
    * Evaluate the transpose to the target matrix.
    *
    * Parameters:
    *     target       Target matrix (must have the dimensions of the transpose)
    ***********************************************************************************************/
    template<typename MAO>
    template<WritableMatrix MAT>
    void Transpose<MAO>::eval(MAT& target) const
    {
        const index_t M = target.rows();
        const index_t N = target.cols();

        #ifdef MATL_DIMENSION_CHECK
        if ((M != this->mat_.cols()) || (N != this->mat_.rows()))
            detail::assignment_dimensions(M, N, this->mat_.cols(), this->mat_.rows());
        #endif

        // Evaluator of the operand. Any products in the operand are computed here, before any
        // element of the target is written.
        const detail::evaluator<MAO> source(this->mat_);

        constexpr index_t TILE = detail::TILE;

        // Square matrices that read the target (e.g. A = transpose(A) or A = transpose(A+A)):
        // elements (i,j) and (j,i) are evaluated simultaneously. Element (i,j) of an elementwise
        // operand depends only on element (i,j) of its matrices (a safe target, see
        // unsafe_evaluation_target, shares elements with the operand only at the same positions),
        // so reading both elements of the pair before writing either is memory safe. The pairs are visited in tiles on and above
        // the diagonal, so that the rows of one tile and the columns of its mirror tile stay in
        // cache.
        if ((M == N) && this->reads(detail::footprint(target)))
        {
            for (index_t ii = 1; ii <= M; ii += TILE)
            {
                const index_t imax = std::min(ii + TILE - 1, M);
                for (index_t jj = ii; jj <= M; jj += TILE)
                {
                    const index_t jmax = std::min(jj + TILE - 1, M);
                    for (index_t i = ii; i <= imax; i++)
                    {
                        if (jj == ii) target(i,i) = source(i,i);
                        for (index_t j = std::max(jj, i+1); j <= jmax; j++)
                        {
                            const dtype source_ij = source(i,j);
                            const dtype source_ji = source(j,i);
                            target(i,j) = source_ji;
                            target(j,i) = source_ij;
                        }
                    }
                }
            }
        }
        // Transposing flips the storage order, so a target with the opposite storage order of the
        // operand can be filled in storage order
        else if constexpr (detail::same_storage_order_v<MAT, Transpose> &&
                           detail::linear_index_v<MAT> && detail::linear_index_v<Transpose>)
        {
            const index_t K = this->size();
            for (index_t i = 1; i <= K; i++)
                target(i) = source(i);
        }
        else
        {
            // Visit the elements in tiles (see detail::TILE)
            for (index_t ii = 1; ii <= M; ii += TILE)
            {
                const index_t imax = std::min(ii + TILE - 1, M);
                for (index_t jj = 1; jj <= N; jj += TILE)
                {
                    const index_t jmax = std::min(jj + TILE - 1, N);
                    for (index_t i = ii; i <= imax; i++)
                        for (index_t j = jj; j <= jmax; j++)
                            target(i,j) = source(j,i);
                }
            }
        }
    }

}

namespace matl::detail
{
   /***********************************************************************************************
    * struct evaluator<expr::Transpose>
    *
    * Evaluates elements of the transpose from the evaluator of its operand. Element (i,j) of the
    * transpose is element (j,i) of the operand. Single-index access passes through unchanged:
    * the transpose has the opposite storage order of its operand, so the i-th element in storage
    * order is the same element.
    ***********************************************************************************************/
    template<typename MAO>
    struct evaluator< expr::Transpose<MAO> >
    {
        typedef expr::Transpose<MAO>             expression_type;
        typedef typename expression_type::dtype dtype;

        evaluator<MAO> operand_;

        explicit evaluator(expression_type const& expr) : operand_(expr.operand()) {}

        dtype operator() (index_t i, index_t j) const {return operand_(j,i);}
        dtype operator() (index_t i)            const {return operand_(i);}
    };
}
#endif
