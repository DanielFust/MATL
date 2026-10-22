/**************************************************************************************************
 * File: _matl_product_implementation_.hpp
 *
 * Description: Implementation of the matrix product expression template.
 **************************************************************************************************/
#ifndef MATL_PRODUCT_IMPLEMENTATION_HPP
#define MATL_PRODUCT_IMPLEMENTATION_HPP
namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * void compute_product(MAT& target, expr::Product const& product)
        *
        * Computes a matrix product into the target matrix. Matrix operands are read directly;
        * expression operands are first evaluated into temporary matrices, which are released as
        * soon as the product has been computed.
        *
        * Products of arithmetic element types use the cache-blocked kernel (see
        * _matl_product_kernel_.hpp), except for small products, where the simple loop is faster.
        * Both add the products of each element in the same order, so they give the same result.
        *
        * Parameters:
        *     target       Target matrix (must not be read by the product, see
        *                  unsafe_evaluation_target)
        *     product      Product expression
        *******************************************************************************************/
        template<WritableMatrix MAT, typename LHS, typename RHS>
        void compute_product(MAT& target, expr::Product<LHS,RHS> const& product)
        {
            typedef typename expr::Product<LHS,RHS>::dtype dtype;

            const product_operand_t<LHS> lhs(product.lhs());
            const product_operand_t<RHS> rhs(product.rhs());

            const index_t M = product.rows();
            const index_t N = product.cols();
            const index_t K = product.lhs().cols();

            // Cache-blocked kernel (the target is written in the product's element type)
            if constexpr (std::is_arithmetic_v<dtype> && std::is_same_v<dtype_t<MAT>, dtype>)
            {
                if (M * N * K >= product_blocking<dtype>::MIN_BLOCKED_WORK)
                {
                    blocked_product<dtype>(target, lhs, rhs, M, N, K);
                    return;
                }
            }

            // Simple loop
            for (index_t i = 1; i <= M; i++)
            {
                for (index_t j = 1; j <= N; j++)
                {
                    dtype sum = 0;
                    for (index_t k = 1; k <= K; k++)
                        sum += lhs(i,k) * rhs(k,j);
                    target(i,j) = sum;
                }
            }
        }
    }

    namespace expr
    {
       /*******************************************************************************************
        * Constructor
        *
        * Stores the operands: expressions are copied, matrices are referenced. Nothing is
        * computed.
        *
        * Parameters:
        *     lhs          Left operand (matrix or expression)
        *     rhs          Right operand (matrix or expression)
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        Product<LHS,RHS>::Product(LHS const& lhs, RHS const& rhs) : lhs_(lhs), rhs_(rhs)
        {}

       /*******************************************************************************************
        * LHS const& lhs() const
        *
        * Returns:
        *     Left operand
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        LHS const& Product<LHS,RHS>::lhs() const
        {return this->lhs_;}

       /*******************************************************************************************
        * RHS const& rhs() const
        *
        * Returns:
        *     Right operand
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        RHS const& Product<LHS,RHS>::rhs() const
        {return this->rhs_;}

       /*******************************************************************************************
        * index_t rows() const
        *
        * Returns:
        *     index_t      Number of rows in resultant matrix
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        index_t Product<LHS,RHS>::rows() const
        {return this->lhs_.rows();}

       /*******************************************************************************************
        * index_t cols() const
        *
        * Returns:
        *     index_t      Number of columns in resultant matrix
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        index_t Product<LHS,RHS>::cols() const
        {return this->rhs_.cols();}

       /*******************************************************************************************
        * index_t size() const
        *
        * Returns:
        *     index_t      Number of elements in resultant matrix
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        index_t Product<LHS,RHS>::size() const
        {return this->lhs_.rows() * this->rhs_.cols();}

       /*******************************************************************************************
        * void check_dimensions() const
        *
        * Checks the operands (recursively), then that the number of columns of the left operand
        * equals the number of rows of the right operand. Reports a matl::dimension_error otherwise.
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        void Product<LHS,RHS>::check_dimensions() const
        {
            detail::check_dimensions(this->lhs_);
            detail::check_dimensions(this->rhs_);

            if (this->lhs_.cols() != this->rhs_.rows())
                detail::incompatible_dimensions("matrix product (columns of the left operand must equal rows of the right operand)",
                                        this->lhs_.rows(), this->lhs_.cols(),
                                        this->rhs_.rows(), this->rhs_.cols());
        }

       /*******************************************************************************************
        * bool unsafe_evaluation_target(detail::Footprint const&)
        *
        * Check if computing the product directly into the target matrix would corrupt memory.
        *
        * The product reads matrix and view operands directly while it writes the target, so it is
        * unsafe if either operand shares elements with the target (e.g. C = A*C, or
        * C(all, range(1, 2)) = C * B). Expression operands are evaluated into temporaries before
        * any element of the target is written, so they are always safe, even if they read the
        * target (e.g. C = (A*C)*B or C = transpose(C)*B).
        *
        * Parameters:
        *     target       Footprint of the target matrix or view
        * Returns:
        *     'true' if computing the product directly into the target would corrupt memory
        *     'false' otherwise
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        bool Product<LHS,RHS>::unsafe_evaluation_target(detail::Footprint const& target) const
        {
            bool unsafe = false;
            if constexpr (MatrixType<LHS> || detail::is_view_v<LHS>) unsafe = unsafe || detail::references(target, this->lhs_);
            if constexpr (MatrixType<RHS> || detail::is_view_v<RHS>) unsafe = unsafe || detail::references(target, this->rhs_);
            return unsafe;
        }

       /*******************************************************************************************
        * bool reads(detail::Footprint const&)
        *
        * Check if the expression may read elements of the target (see detail::Footprint).
        *
        * Parameters:
        *     target       Footprint of the target (where its elements are in memory)
        * Returns:
        *     'true' if the expression references the target matrix, 'false' otherwise
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        bool Product<LHS,RHS>::reads(detail::Footprint const& target) const
        {return detail::references(target, this->lhs_) || detail::references(target, this->rhs_);}

       /*******************************************************************************************
        * void eval(MAT& target) const
        *
        * Computes the product directly into the target matrix.
        *
        * Parameters:
        *     target       Target matrix (must be a safe evaluation target with matching dimensions)
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        template<WritableMatrix MAT>
        void Product<LHS,RHS>::eval(MAT& target) const
        {detail::compute_product(target, *this);}
    }

    namespace detail
    {
       /*******************************************************************************************
        * struct evaluator<expr::Product>
        *
        * Evaluator for a product nested inside a larger expression. The product is computed into
        * a temporary matrix when the evaluator is constructed (before any element of the
        * assignment target is written), and the temporary is released when the evaluator is
        * destroyed at the end of the assignment.
        *******************************************************************************************/
        template<typename LHS, typename RHS>
        struct evaluator< expr::Product<LHS,RHS> >
        {
            typedef expr::Product<LHS,RHS>             expression_type;
            typedef typename expression_type::dtype    dtype;
            typedef result_matrix_t<expression_type>   result_type;

            EvaluationTemporary<result_type> result_;

            explicit evaluator(expression_type const& expr) : result_(expr.rows(), expr.cols())
            {compute_product(result_.matrix(), expr);}

            dtype operator() (index_t i, index_t j) const {return result_.matrix()(i,j);}
            dtype operator() (index_t i)            const {return result_.matrix()(i);}
        };
    }
}
#endif
