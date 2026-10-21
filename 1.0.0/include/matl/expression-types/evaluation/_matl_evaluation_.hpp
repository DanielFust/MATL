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
 *   - evaluators of elementwise expressions combine the elements of their child evaluators
 *   - evaluator<_MATRIX_TIMES_MATRIX_> computes the product into a temporary matrix when it is
 *     constructed, and then reads elements of that temporary
 *
 * Evaluators live only for the duration of one assignment, so every temporary is released when
 * the assignment finishes, in the reverse order of allocation. Temporaries are held by
 * EvaluationTemporary objects, which take their memory from this thread's TemporaryArena (see
 * memory-management/temporary-arena) so that repeated evaluations recycle the same memory. Because all products are fully
 * computed before any element of the target is written, products nested inside an expression
 * never cause aliasing problems, even if they read the target matrix.
 *
 * The evaluator<> specialization for each expression type is defined alongside that expression
 * type's implementation.
 **************************************************************************************************/
#ifndef _MATL_EVALUATION__HPP
#define _MATL_EVALUATION__HPP
namespace matl
{
   /***********************************************************************************************
    * struct evaluator<MAT> (matrices)
    *
    * Reads elements of a matrix. Holds a reference to the matrix.
    ***********************************************************************************************/
    template<typename MAT> requires MatrixType<MAT>
    struct evaluator<MAT>
    {
        typedef typename type_traits<MAT>::dtype dtype;

        MAT const& mat_;

        explicit evaluator(MAT const& mat) : mat_(mat) {}

        dtype operator() (index_t i, index_t j) const {return mat_(i,j);}
        dtype operator() (index_t i)            const {return mat_(i);}
    };


   /***********************************************************************************************
    * is_owned_matrix<T>
    *
    * True for _OWNED_MATRIX_ expressions, which hold a complete matrix.
    ***********************************************************************************************/
    template<typename T>
    struct is_owned_matrix : std::false_type {};

    template<typename MAT>
    struct is_owned_matrix< _OWNED_MATRIX_<MAT> > : std::true_type {};


   /***********************************************************************************************
    * product_operand_t<X>
    *
    * How a matrix product reads an operand of type X. Each operand element is read many times
    * during a product, so operands that are matrices (or owned matrices) are read directly, while
    * any other expression is first evaluated into a temporary matrix.
    ***********************************************************************************************/
    template<typename X>
    class materialized_operand
    {
      private:
        EvaluationTemporary<typename is_fixed_size<X>::matrix_type> temp_;

      public:
        typedef typename type_traits<X>::dtype dtype;

        // Evaluate the expression into a temporary (a fresh matrix, so it is a safe target)
        explicit materialized_operand(X const& expr) : temp_(expr.rows(), expr.cols())
        {expr.eval(this->temp_.matrix());}

        dtype operator() (index_t i, index_t j) const {return this->temp_.matrix()(i,j);}
        dtype operator() (index_t i)            const {return this->temp_.matrix()(i);}
    };

    template<typename X>
    using product_operand_t = std::conditional_t< MatrixType<X> || is_owned_matrix<X>::value,
                                                  evaluator<X>,
                                                  materialized_operand<X> >;


   /***********************************************************************************************
    * void assign_elementwise(_MATRIX_<MAT>& target, EXPR const& expr)
    *
    * Evaluates an elementwise expression into the target matrix. Elements are visited in storage
    * order when the target and the expression share a storage order, and by (row, column)
    * otherwise.
    *
    * The caller is responsible for ensuring the target is a safe evaluation target (see
    * unsafe_evaluation_target) and has the dimensions of the expression.
    ***********************************************************************************************/
    template<typename MAT, typename EXPR>
    void assign_elementwise(_MATRIX_<MAT>& target, EXPR const& expr)
    {
        const evaluator<EXPR> source(expr);

        if constexpr (is_same_storage_order<MAT, EXPR>::value)
        {
            const index_t N = expr.size();
            for (index_t i = 1; i <= N; ++i)
                target(i) = source(i);
        }
        else
        {
            const index_t M = expr.rows();
            const index_t N = expr.cols();
            for (index_t i = 1; i <= M; ++i)
                for (index_t j = 1; j <= N; ++j)
                    target(i,j) = source(i,j);
        }
    }
}
#endif
