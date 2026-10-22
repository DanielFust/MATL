/**************************************************************************************************
 * File: _matl_evaluation_.hpp
 *
 * Description: Evaluators and the shared routines used to evaluate expressions.
 *
 * Expression objects are pure formulas: they store their operands and compute nothing. Work
 * happens only when an expression is assigned to a matrix. The assignment builds a tree of
 * evaluator objects that mirrors the expression tree:
 *
 *   - evaluator<MATRIX> reads elements of a matrix directly
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
    * bool references(void const* target, X const& operand)
    *
    * Returns:
    *     'true' if the operand is, or reads, the matrix at address 'target'
    ***********************************************************************************************/
    template<typename X>
    bool references(void const* target, X const& operand)
    {
        if constexpr (Scalar<X>)          return false;
        else if constexpr (MatrixType<X>) return static_cast<void const*>(&operand) == target;
        else                              return operand.references_address(target);
    }

   /***********************************************************************************************
    * bool unsafe_operand(MAT const& target, X const& operand)
    *
    * Returns:
    *     'true' if evaluating the operand as part of an expression assigned directly to the target
    *     matrix would overwrite elements of the target that the operand still has to read.
    *     Scalars and matrices are always safe on their own, and so are products, which are
    *     computed into a temporary (see evaluator<expr::Product>) before any element of the
    *     target is written.
    ***********************************************************************************************/
    template<MatrixType MAT, typename X>
    bool unsafe_operand(MAT const& target, X const& operand)
    {
        if constexpr (Scalar<X> || MatrixType<X> || is_product_v<X>) return false;
        else                                                       return operand.unsafe_evaluation_target(target);
    }

   /***********************************************************************************************
    * bool unsafe_elementwise_operand(MAT const& target, X const& operand)
    *
    * Returns:
    *     'true' if the operand makes an elementwise expression unsafe to evaluate directly into
    *     the target matrix: either it is unsafe on its own, or it reads the target through a
    *     transpose (element (i,j) of the result would read element (j,i) of the target, which
    *     may already have been overwritten).
    ***********************************************************************************************/
    template<MatrixType MAT, typename X>
    bool unsafe_elementwise_operand(MAT const& target, X const& operand)
    {
        return (is_transpose_v<X> && references(&target, operand)) || unsafe_operand(target, operand);
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
    using product_operand_t = std::conditional_t< MatrixType<X> || is_owned_v<X>,
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
    * Evaluates an elementwise expression into the target matrix. Elements are visited in storage
    * order when the target and the expression share a storage order, and in tiles otherwise (see
    * TILE).
    *
    * The caller is responsible for ensuring the target is a safe evaluation target (see
    * unsafe_evaluation_target) and has the dimensions of the expression.
    ***********************************************************************************************/
    template<MatrixType MAT, typename EXPR>
    void assign_elementwise(MAT& target, EXPR const& expr)
    {
        const evaluator<EXPR> source(expr);

        if constexpr (same_storage_order_v<MAT, EXPR>)
        {
            const index_t N = expr.size();
            for (index_t i = 1; i <= N; ++i)
                target(i) = source(i);
        }
        else
        {
            const index_t M = expr.rows();
            const index_t N = expr.cols();
            for (index_t ii = 1; ii <= M; ii += TILE)
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
        }
    }
}
#endif
