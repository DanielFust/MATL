/**************************************************************************************************
 * File: _matl_sparse_product_.hpp
 *
 * Description: Products of a sparse matrix with a dense matrix or vector:
 *
 *     y = S * x;                  sparse matrix times dense matrix (or vector, an n x 1 matrix)
 *     y = transpose(S) * x;       the transpose of S times x, without forming the transpose
 *     r = b - S * x;              products are expressions, and nest like dense products
 *
 * The result is dense. Like a dense product, S * x is an expression: it is computed when it is
 * assigned (or, nested in a larger expression, into a temporary before the rest is evaluated), and
 * an assignment whose target is read by the product (x = S * x) goes through a temporary. The
 * sparse matrix is referenced by the expression, so a stored expression reads its current values.
 *
 * Cost: one multiply-add per stored entry and column of x. S * x visits the rows of S
 * independently, and runs on several threads with OpenMP (see _matl_parallel_.hpp);
 * transpose(S) * x scatters into the rows of the result, and runs serially (an explicit transpose,
 * CSRMatrix St = transpose(S), makes St * x parallel).
 *
 * The dimensions are always checked (matl::dimension_error).
 **************************************************************************************************/
#ifndef MATL_SPARSE_PRODUCT_HPP
#define MATL_SPARSE_PRODUCT_HPP
namespace matl
{
    namespace expr
    {
        // transpose(S) of a sparse matrix: a reference to S, used as the left operand of a product
        template<class T, class I>
        class SparseTransposed
        {
          private:
            CSRMatrix<T, I> const& matrix_;

          public:
            explicit SparseTransposed(CSRMatrix<T, I> const& matrix) : matrix_(matrix) {}

            // The sparse matrix that is transposed
            CSRMatrix<T, I> const& matrix() const {return this->matrix_;}

            index_t rows() const {return this->matrix_.cols();}
            index_t cols() const {return this->matrix_.rows();}
        };

        // Sparse matrix (transposed if TRANSPOSED) times a dense matrix, view or expression X
        template<class T, class I, bool TRANSPOSED, typename X> class SparseProduct;
    }

    namespace detail
    {
        // Traits of a sparse product (see _matl_type_traits_.hpp): dense, dynamic dimensions,
        // computed row by row
        template<class T, class I, bool TRANSPOSED, typename X>
        struct traits< expr::SparseProduct<T, I, TRANSPOSED, X> >
        {
            typedef promote_t<T, dtype_t<X>> dtype;
            static constexpr bool         fixed_size    = false;
            static constexpr index_t      rows          = 0;
            static constexpr index_t      cols          = 0;
            static constexpr StorageOrder storage_order = StorageOrder::RowMajor;
        };

        // Like a dense product, a nested sparse product is computed into a temporary before the
        // rest of the expression is evaluated (see unsafe_operand)
        template<class T, class I, bool TRANSPOSED, typename X>
        inline constexpr bool is_product_v< expr::SparseProduct<T, I, TRANSPOSED, X> > = true;
    }

    namespace expr
    {
       /*******************************************************************************************
        * class SparseProduct<T, I, TRANSPOSED, X>
        *
        * Expression for S * X (or transpose(S) * X), with S a CSRMatrix<T, I> and X a dense
        * matrix, view or expression (see the description above).
        *******************************************************************************************/
        template<class T, class I, bool TRANSPOSED, typename X>
        class SparseProduct : public ExpressionBase< SparseProduct<T, I, TRANSPOSED, X> >
        {
          public:
            typedef typename detail::traits<SparseProduct>::dtype dtype;

          private:
            CSRMatrix<T, I> const& matrix_;
            detail::stored_t<X>    x_;            // matrices by reference, expressions by value

          public:
            SparseProduct(CSRMatrix<T, I> const& matrix, X const& x) : matrix_(matrix), x_(x) {}

            CSRMatrix<T, I> const& matrix() const {return this->matrix_;}
            X const& operand() const {return this->x_;}

            // Dimensions of the result
            index_t rows() const {return TRANSPOSED ? this->matrix_.cols() : this->matrix_.rows();}
            index_t cols() const {return this->x_.cols();}
            index_t size() const {return this->rows() * this->cols();}

            // Rows of x the product needs (the inner dimension)
            index_t inner() const {return TRANSPOSED ? this->matrix_.rows() : this->matrix_.cols();}

            // Checks the operand (recursively) and the inner dimension
            void check_dimensions() const
            {
                detail::check_dimensions(this->x_);
                if (this->inner() != this->x_.rows())
                    detail::incompatible_dimensions(TRANSPOSED ? "sparse matrix product transpose(S) * X"
                                                               : "sparse matrix product S * X",
                                                    TRANSPOSED ? this->matrix_.cols() : this->matrix_.rows(),
                                                    this->inner(), this->x_.rows(), this->x_.cols());
            }

            // Unsafe if X is a matrix or view that shares elements with the target (X is read while
            // the target is written); expression operands are evaluated into temporaries first
            bool unsafe_evaluation_target(detail::Footprint const& target) const
            {
                if constexpr (MatrixType<X> || detail::is_view_v<X>) return detail::references(target, this->x_);
                else                                                 return false;
            }

            // Check if the expression may read elements of the target
            bool reads(detail::Footprint const& target) const
            {return detail::references(target, this->x_);}

            // Computes the product into the target (a safe target with the product's dimensions)
            template<WritableMatrix MAT>
            void eval(MAT& target) const;
        };
    }

    namespace detail
    {
       /*******************************************************************************************
        * void compute_sparse_product(MAT& target, expr::SparseProduct const& product)
        *
        * target = S * X (row by row: target(i,:) = sum over the stored entries (i,j) of
        * S(i,j) X(j,:)), or target = S^T X (scattering: target(j,:) += S(i,j) X(i,:)). X is read
        * directly if it is a matrix or view, and evaluated into a temporary otherwise.
        *******************************************************************************************/
        template<WritableMatrix MAT, class T, class I, bool TRANSPOSED, typename X>
        void compute_sparse_product(MAT& target, expr::SparseProduct<T, I, TRANSPOSED, X> const& product)
        {
            typedef typename expr::SparseProduct<T, I, TRANSPOSED, X>::dtype dtype;
            product.check_dimensions();                        // always: cheap next to the product

            CSRMatrix<T, I> const& s = product.matrix();
            const product_operand_t<X> x(product.operand());
            const std::span<const I> offsets = s.row_offsets();
            const std::span<const I> columns = s.column_indices();
            const std::span<const T> values  = s.values();
            const index_t M = s.rows();
            const index_t K = product.cols();

            if constexpr (!TRANSPOSED)
            {
                // Rows are independent: with OpenMP, each thread computes a range of rows holding
                // about the same number of stored entries (rows can differ in length)
                const index_t nnz = s.nnz();
                const auto by_entries = [offsets, nnz](index_t t, index_t threads)
                {
                    const I target = static_cast<I>(nnz / threads * t + std::min(t, nnz % threads));
                    return static_cast<index_t>(std::lower_bound(offsets.begin(), offsets.end(), target) - offsets.begin());
                };
                parallel_ranges(M, nnz * K, [&](index_t first, index_t last)
                {
                    if (K == 1)
                    {
                        for (index_t i = first + 1; i <= last; ++i)
                        {
                            dtype sum = 0;
                            for (I p = offsets[i - 1]; p < offsets[i]; ++p)
                                sum += values[p] * x(static_cast<index_t>(columns[p]) + 1, 1);
                            target(i, 1) = sum;
                        }
                    }
                    else
                    {
                        std::vector<dtype> row(K);
                        for (index_t i = first + 1; i <= last; ++i)
                        {
                            std::fill(row.begin(), row.end(), dtype(0));
                            for (I p = offsets[i - 1]; p < offsets[i]; ++p)
                            {
                                const T v = values[p];
                                const index_t j = static_cast<index_t>(columns[p]) + 1;
                                for (index_t c = 1; c <= K; ++c) row[c - 1] += v * x(j, c);
                            }
                            for (index_t c = 1; c <= K; ++c) target(i, c) = row[c - 1];
                        }
                    }
                }, by_entries);
            }
            else
            {
                // Serial: the scattered updates of different rows of S would write the same rows
                // of the result (form the transpose, CSRMatrix<T, I> St = transpose(S), for a
                // product that runs in parallel)
                for (index_t j = 1; j <= product.rows(); ++j)
                    for (index_t c = 1; c <= K; ++c) target(j, c) = dtype(0);
                for (index_t i = 1; i <= M; ++i)
                    for (I p = offsets[i - 1]; p < offsets[i]; ++p)
                    {
                        const T v = values[p];
                        const index_t j = static_cast<index_t>(columns[p]) + 1;
                        for (index_t c = 1; c <= K; ++c) target(j, c) += v * x(i, c);
                    }
            }
        }

       /*******************************************************************************************
        * struct evaluator<expr::SparseProduct>
        *
        * A sparse product nested in a larger expression is computed into a temporary when the
        * evaluator is constructed (before any element of the target is written).
        *******************************************************************************************/
        template<class T, class I, bool TRANSPOSED, typename X>
        struct evaluator< expr::SparseProduct<T, I, TRANSPOSED, X> >
        {
            typedef expr::SparseProduct<T, I, TRANSPOSED, X> expression_type;
            typedef typename expression_type::dtype          dtype;
            typedef result_matrix_t<expression_type>         result_type;

            EvaluationTemporary<result_type> result_;

            explicit evaluator(expression_type const& expr) : result_(expr.rows(), expr.cols())
            {compute_sparse_product(result_.matrix(), expr);}

            dtype operator() (index_t i, index_t j) const {return result_.matrix()(i,j);}
            dtype operator() (index_t i)            const {return result_.matrix()(i);}
        };
    }

    template<class T, class I, bool TRANSPOSED, typename X>
    template<WritableMatrix MAT>
    void expr::SparseProduct<T, I, TRANSPOSED, X>::eval(MAT& target) const
    {detail::compute_sparse_product(target, *this);}


    /******************************************* Operators ***************************************/

    // Sparse matrix times a dense matrix, view or expression
    template<class T, class I, MatrixAlgebraObject X>
    expr::SparseProduct<T, I, false, X> operator* (CSRMatrix<T, I> const& s, X const& x)
    {return expr::SparseProduct<T, I, false, X>(s, x);}

    // Sparse matrix times a temporary dense matrix (moved into the expression)
    template<class T, class I, TemporaryMatrix M>
    auto operator* (CSRMatrix<T, I> const& s, M&& m)
    {return s * own_temporary(std::forward<M>(m));}

    // transpose(S) of a sparse matrix, for products: transpose(S) * x
    template<class T, class I>
    expr::SparseTransposed<T, I> transpose(CSRMatrix<T, I> const& s)
    {return expr::SparseTransposed<T, I>(s);}

    // Transpose of a sparse matrix times a dense matrix, view or expression
    template<class T, class I, MatrixAlgebraObject X>
    expr::SparseProduct<T, I, true, X> operator* (expr::SparseTransposed<T, I> const& st, X const& x)
    {return expr::SparseProduct<T, I, true, X>(st.matrix(), x);}

    template<class T, class I, TemporaryMatrix M>
    auto operator* (expr::SparseTransposed<T, I> const& st, M&& m)
    {return st * own_temporary(std::forward<M>(m));}
}
#endif
