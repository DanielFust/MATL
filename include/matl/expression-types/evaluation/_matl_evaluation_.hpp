/**************************************************************************************************
 * File: _matl_evaluation_.hpp
 *
 * Description: Evaluators and the shared routines used to evaluate expressions.
 *
 * Expression objects are pure formulas: they store their operands and compute nothing. Work
 * happens only when an expression is assigned to a matrix. The assignment builds a tree of
 * evaluator objects that mirrors the expression tree:
 *
 *   - evaluator<MATRIX> reads elements of a matrix directly, and evaluator<MatrixView> elements
 *     of a view (see matrix-types/matrixview)
 *   - evaluator<SCALAR> returns the scalar for every element
 *   - evaluators of elementwise expressions combine the elements of their child evaluators
 *   - evaluator<expr::Product> computes the product into a temporary matrix when it is
 *     constructed, and then reads elements of that temporary
 *
 * Evaluators live only for the duration of one assignment, so every temporary is released when
 * the assignment finishes, in the reverse order of allocation. Temporaries are held by
 * EvaluationTemporary objects, which take their memory from this thread's TemporaryArena (see
 * memory-management/temporary-arena) so that repeated evaluations recycle the same memory.
 * Because all products are fully computed before any element of the target is written, products
 * nested inside an expression never cause aliasing problems, even if they read the target matrix.
 *
 * The evaluator specialization for each expression type is defined alongside that expression
 * type's implementation.
 **************************************************************************************************/
#ifndef MATL_EVALUATION_HPP
#define MATL_EVALUATION_HPP
namespace matl::detail
{
   /***********************************************************************************************
    * struct evaluator<MAT> (matrices)
    *
    * Reads elements of a matrix. Holds a reference to the matrix.
    ***********************************************************************************************/
    template<MatrixType MAT>
    struct evaluator<MAT>
    {
        typedef dtype_t<MAT> dtype;

        MAT const& mat_;

        explicit evaluator(MAT const& mat) : mat_(mat) {}

        dtype operator() (index_t i, index_t j) const {return mat_(i,j);}
        dtype operator() (index_t i)            const {return mat_(i);}
    };

   /***********************************************************************************************
    * struct evaluator<S> (scalars)
    *
    * Every element of a scalar operand is the scalar itself.
    ***********************************************************************************************/
    template<Scalar S>
    struct evaluator<S>
    {
        typedef S dtype;

        S value_;

        explicit evaluator(S value) : value_(value) {}

        dtype operator() (index_t, index_t) const {return value_;}
        dtype operator() (index_t)          const {return value_;}
    };


   /***********************************************************************************************
    * bool references(Footprint const& target, X const& operand)
    *
    * Returns:
    *     'true' if the operand reads memory that may hold elements of the target (see Footprint)
    ***********************************************************************************************/
    template<typename X>
    bool references(Footprint const& target, X const& operand)
    {
        if constexpr (Scalar<X>)                          return false;
        else if constexpr (MatrixType<X> || is_view_v<X>) return overlap(target, footprint(operand));
        else                                              return operand.reads(target);
    }

   /***********************************************************************************************
    * bool unsafe_operand(Footprint const& target, X const& operand)
    *
    * Returns:
    *     'true' if evaluating the operand as part of an expression assigned directly to the target
    *     would overwrite elements of the target that the operand still has to read.
    *
    *     - Scalars are always safe, and so are products, which are computed into a temporary (see
    *       evaluator<expr::Product>) before any element of the target is written.
    *     - A matrix or view is read element by element at the position that is being written. It
    *       is safe if it shares no elements with the target, or if its element (i,j) is the
    *       target's element (i,j) (e.g. A = A + B: each element is read before it is written).
    *       Otherwise (e.g. a view of other rows of the target) it is unsafe.
    *     - Other expressions check their operands (see unsafe_evaluation_target).
    ***********************************************************************************************/
    template<typename X>
    bool unsafe_operand(Footprint const& target, X const& operand)
    {
        if constexpr (Scalar<X> || is_product_v<X>) return false;
        else if constexpr (MatrixType<X> || is_view_v<X>)
        {
            const Footprint operand_footprint = footprint(operand);
            return overlap(target, operand_footprint) && !same_elements(target, operand_footprint);
        }
        else return operand.unsafe_evaluation_target(target);
    }

   /***********************************************************************************************
    * bool unsafe_elementwise_operand(Footprint const& target, X const& operand)
    *
    * Returns:
    *     'true' if the operand makes an elementwise expression unsafe to evaluate directly into
    *     the target: either it is unsafe on its own, or it reads the target through a transpose
    *     (element (i,j) of the result would read element (j,i) of the target, which may already
    *     have been overwritten).
    ***********************************************************************************************/
    template<typename X>
    bool unsafe_elementwise_operand(Footprint const& target, X const& operand)
    {
        return (is_transpose_v<X> && references(target, operand)) || unsafe_operand(target, operand);
    }


   /***********************************************************************************************
    * void check_dimensions(X const& operand)
    *
    * Checks that the dimensions of every operation in an expression fit (the operands of an
    * elementwise operation have the same dimensions; the inner dimensions of a product agree).
    * Expressions check their operands recursively. Matrices, scalars and expression types without
    * a check_dimensions() method have nothing to check. Errors are reported as matl::dimension_error
    * (see _matl_errors_.hpp).
    *
    * Dimensions known at compile time are already checked by the operators (see SizesMatch). This
    * run-time check covers dynamic dimensions, including matrices that were reallocated after a
    * stored expression was built. It runs when an expression is assigned, if MATL_DIMENSION_CHECK
    * is defined.
    ***********************************************************************************************/
    template<typename X>
    void check_dimensions(X const& operand)
    {
        if constexpr (requires { operand.check_dimensions(); })
            operand.check_dimensions();
    }


   /***********************************************************************************************
    * class materialized_operand<X>
    *
    * An expression evaluated into a temporary matrix, for operands whose elements are read many
    * times (operands of a matrix product).
    ***********************************************************************************************/
    template<typename X>
    class materialized_operand
    {
      private:
        EvaluationTemporary< result_matrix_t<X> > temp_;

      public:
        typedef dtype_t<X> dtype;

        // Evaluate the expression into a temporary (a fresh matrix, so it is a safe target)
        explicit materialized_operand(X const& expr) : temp_(expr.rows(), expr.cols())
        {expr.eval(this->temp_.matrix());}

        dtype operator() (index_t i, index_t j) const {return this->temp_.matrix()(i,j);}
        dtype operator() (index_t i)            const {return this->temp_.matrix()(i);}
    };

   /***********************************************************************************************
    * product_operand_t<X>
    *
    * How a matrix product reads an operand of type X. Each operand element is read many times
    * during a product, so operands that are matrices (or owned matrices) are read directly, while
    * any other expression is first evaluated into a temporary matrix.
    ***********************************************************************************************/
    template<typename X>
    using product_operand_t = std::conditional_t< MatrixType<X> || is_view_v<X> || is_owned_v<X>,
                                                  evaluator<X>,
                                                  materialized_operand<X> >;


   /***********************************************************************************************
    * TILE
    *
    * When elements are written in one storage order and read in the other (e.g. a row-major
    * target and a transposed operand), they are visited in TILE x TILE tiles: the rows and columns
    * of a tile both stay in cache, whereas visiting whole rows would read the operand one cache
    * line per element. (16 measured best: larger tiles of power-of-two sized matrices map their
    * columns to too few cache sets.)
    ***********************************************************************************************/
    inline constexpr index_t TILE = 16;

   /***********************************************************************************************
    * void assign_elementwise(MAT& target, EXPR const& expr)
    *
    * Evaluates an elementwise expression into the target (a matrix or a view). Elements are
    * visited in storage order when the target and the expression share a storage order (with a
    * single index when both allow it, see linear_index_v), and in tiles otherwise (see TILE).
    * Large evaluations run on several threads with OpenMP (see _matl_parallel_.hpp).
    *
    * The caller is responsible for ensuring the target is a safe evaluation target (see
    * unsafe_evaluation_target) and has the dimensions of the expression.
    ***********************************************************************************************/
    template<WritableMatrix MAT, typename EXPR>
    void assign_elementwise(MAT& target, EXPR const& expr)
    {
        const evaluator<EXPR> source(expr);
        const index_t M = expr.rows();
        const index_t N = expr.cols();

        // Each thread writes its own range of rows, columns or tiles (see parallel_ranges); the
        // result does not depend on the number of threads
        if constexpr (same_storage_order_v<MAT, EXPR> && linear_index_v<MAT> && linear_index_v<EXPR>)
        {
            const index_t K = expr.size();
            parallel_ranges(K, K, [&](index_t first, index_t last)
            {
                for (index_t k = first + 1; k <= last; ++k)
                    target(k) = source(k);
            });
        }
        else if constexpr (same_storage_order_v<MAT, EXPR>)
        {
            // Storage order by (row, column): row by row, or column by column
            if constexpr (traits<MAT>::storage_order == StorageOrder::RowMajor)
            {
                parallel_ranges(M, M * N, [&](index_t first, index_t last)
                {
                    for (index_t i = first + 1; i <= last; ++i)
                        for (index_t j = 1; j <= N; ++j)
                            target(i,j) = source(i,j);
                });
            }
            else
            {
                parallel_ranges(N, M * N, [&](index_t first, index_t last)
                {
                    for (index_t j = first + 1; j <= last; ++j)
                        for (index_t i = 1; i <= M; ++i)
                            target(i,j) = source(i,j);
                });
            }
        }
        else
        {
            // Rows of tiles
            parallel_ranges((M + TILE - 1) / TILE, M * N, [&](index_t first, index_t last)
            {
                for (index_t ii = first * TILE + 1; ii <= std::min(last * TILE, M); ii += TILE)
                {
                    const index_t imax = std::min(ii + TILE - 1, M);
                    for (index_t jj = 1; jj <= N; jj += TILE)
                    {
                        const index_t jmax = std::min(jj + TILE - 1, N);
                        for (index_t i = ii; i <= imax; ++i)
                            for (index_t j = jj; j <= jmax; ++j)
                                target(i,j) = source(i,j);
                    }
                }
            });
        }
    }
}
#endif
