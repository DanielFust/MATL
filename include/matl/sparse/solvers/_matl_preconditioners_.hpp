/**************************************************************************************************
 * File: _matl_preconditioners_.hpp
 *
 * Description: Preconditioners for the iterative solvers (see _matl_iterative_solvers_.hpp).
 *
 * A preconditioner M approximates A, and is cheap to apply: z = M^-1 r. The iterative solvers then
 * work on a better-conditioned system and need fewer iterations.
 *
 *     IdentityPreconditioner         M = I (no preconditioning; the default)
 *     JacobiPreconditioner(A)        M = diag(A): cheap, and symmetric (for CG)
 *     ILU0Preconditioner(S)          M = L U, the incomplete LU factorization of a sparse matrix
 *                                    with the pattern of S (no fill-in), as in most PDE codes; for
 *                                    a symmetric S it is symmetric too (L D L^T), so CG can use it
 *
 * Any type with a member apply(r, z) that sets z = M^-1 r (both n x 1 matrices) can be passed as
 * a preconditioner.
 **************************************************************************************************/
#ifndef MATL_PRECONDITIONERS_HPP
#define MATL_PRECONDITIONERS_HPP
namespace matl
{
   /***********************************************************************************************
    * struct IdentityPreconditioner
    *
    * No preconditioning: z = r.
    ***********************************************************************************************/
    struct IdentityPreconditioner
    {
        template<typename R, typename Z>
        void apply(R const& r, Z& z) const {z = r;}
    };


   /***********************************************************************************************
    * class JacobiPreconditioner<T>
    *
    * M = diag(A): z(i) = r(i) / A(i,i). Made from a sparse or dense square matrix; a zero on the
    * diagonal is reported as matl::singular_matrix_error.
    ***********************************************************************************************/
    template<class T>
    class JacobiPreconditioner
    {
        static_assert(std::is_floating_point_v<T>, "JacobiPreconditioner: the element type must be a floating-point type");

      private:
        std::vector<T> inverse_diagonal_;

        void set(index_t n, auto const& diagonal)
        {
            this->inverse_diagonal_.resize(n);
            for (index_t i = 1; i <= n; ++i)
            {
                const T d = static_cast<T>(diagonal(i));
                if (d == T(0))
                    detail::raise<singular_matrix_error>("JacobiPreconditioner: diagonal element " + std::to_string(i) + " is 0");
                this->inverse_diagonal_[i - 1] = T(1) / d;
            }
        }

      public:
        // From a sparse matrix
        template<class U, class I>
        explicit JacobiPreconditioner(CSRMatrix<U, I> const& a)
        {
            detail::check_square("JacobiPreconditioner", a.rows(), a.cols());
            this->set(a.rows(), [&a](index_t i) {return a(i, i);});
        }

        // From a dense matrix, view or expression
        template<MatrixAlgebraObject A>
        explicit JacobiPreconditioner(A const& a)
        {
            detail::check_square("JacobiPreconditioner", a.rows(), a.cols());
            const detail::product_operand_t<A> source(a);
            this->set(a.rows(), [&source](index_t i) {return source(i, i);});
        }

        index_t size() const {return this->inverse_diagonal_.size();}

        // z = M^-1 r
        template<typename R, typename Z>
        void apply(R const& r, Z& z) const
        {
            for (index_t i = 1; i <= this->size(); ++i) z(i, 1) = r(i, 1) * this->inverse_diagonal_[i - 1];
        }
    };

    template<class U, class I> JacobiPreconditioner(CSRMatrix<U, I> const&) -> JacobiPreconditioner<detail::solver_t<U>>;
    template<MatrixAlgebraObject A> JacobiPreconditioner(A const&) -> JacobiPreconditioner<detail::solver_t<detail::dtype_t<A>>>;


   /***********************************************************************************************
    * class ILU0Preconditioner<T, I>
    *
    * Incomplete LU factorization with zero fill-in, ILU(0): L U ~ S, where L (unit lower
    * triangular) and U (upper triangular) have nonzeros only where S does. Computed row by row
    * (the IKJ variant, as in Saad, "Iterative Methods for Sparse Linear Systems", Algorithm 10.4),
    * in one CSR matrix holding L below the diagonal and U on and above it. Every diagonal entry
    * of S must be stored; a zero pivot is reported as matl::singular_matrix_error.
    *
    * apply(r, z) solves L U z = r by forward and back substitution on the pattern.
    ***********************************************************************************************/
    template<class T, class I = index_t>
    class ILU0Preconditioner
    {
        static_assert(std::is_floating_point_v<T>, "ILU0Preconditioner: the element type must be a floating-point type");

      private:
        CSRMatrix<T, I>      lu_;           // L below the diagonal (unit diagonal not stored), U on and above
        std::vector<index_t> diagonal_;     // position of each row's diagonal entry in lu_

      public:
        template<class U>
        explicit ILU0Preconditioner(CSRMatrix<U, I> const& s)
        {
            detail::check_square("ILU0Preconditioner", s.rows(), s.cols());
            const index_t n = s.rows();
            std::vector<T> values(s.values().begin(), s.values().end());
            const auto offsets = s.row_offsets();
            const auto columns = s.column_indices();

            // The diagonal entries, and for each row, where its stored columns are (for the updates)
            this->diagonal_.assign(n, 0);
            for (index_t i = 0; i < n; ++i)
            {
                bool found = false;
                for (I p = offsets[i]; p < offsets[i + 1]; ++p)
                    if (static_cast<index_t>(columns[p]) == i) {this->diagonal_[i] = static_cast<index_t>(p); found = true;}
                if (!found)
                    detail::raise<singular_matrix_error>("ILU0Preconditioner: diagonal entry " + std::to_string(i + 1) +
                                                         " is not stored");
            }

            // IKJ: for each row i, eliminate with the rows k < i that row i has entries in,
            // updating only entries in row i's pattern
            std::vector<index_t> where(n, static_cast<index_t>(-1));     // position of column j in row i
            for (index_t i = 0; i < n; ++i)
            {
                for (I p = offsets[i]; p < offsets[i + 1]; ++p) where[static_cast<index_t>(columns[p])] = static_cast<index_t>(p);
                for (I p = offsets[i]; p < offsets[i + 1]; ++p)
                {
                    const index_t k = static_cast<index_t>(columns[p]);
                    if (k >= i) break;                                     // columns are sorted
                    const T pivot = values[this->diagonal_[k]];
                    if (pivot == T(0))
                        detail::raise<singular_matrix_error>("ILU0Preconditioner: pivot " + std::to_string(k + 1) + " is 0");
                    const T l = values[p] / pivot;
                    values[p] = l;
                    for (I q = static_cast<I>(this->diagonal_[k] + 1); q < offsets[k + 1]; ++q)
                    {
                        const index_t w = where[static_cast<index_t>(columns[q])];
                        if (w != static_cast<index_t>(-1)) values[w] -= l * values[q];
                    }
                }
                if (values[this->diagonal_[i]] == T(0))
                    detail::raise<singular_matrix_error>("ILU0Preconditioner: pivot " + std::to_string(i + 1) + " is 0");
                for (I p = offsets[i]; p < offsets[i + 1]; ++p) where[static_cast<index_t>(columns[p])] = static_cast<index_t>(-1);
            }

            this->lu_ = CSRMatrix<T, I>(n, n, std::vector<I>(offsets.begin(), offsets.end()),
                                        std::vector<I>(columns.begin(), columns.end()), std::move(values));
        }

        index_t size() const {return this->lu_.rows();}

        // The factors: L below the diagonal (unit diagonal not stored), U on and above it
        CSRMatrix<T, I> const& factors() const {return this->lu_;}

        // z = (L U)^-1 r
        template<typename R, typename Z>
        void apply(R const& r, Z& z) const
        {
            const index_t n = this->size();
            const auto offsets = this->lu_.row_offsets();
            const auto columns = this->lu_.column_indices();
            const auto values  = this->lu_.values();
            // L y = r (unit lower)
            for (index_t i = 0; i < n; ++i)
            {
                T s = r(i + 1, 1);
                for (I p = offsets[i]; static_cast<index_t>(p) < this->diagonal_[i]; ++p)
                    s -= values[p] * z(static_cast<index_t>(columns[p]) + 1, 1);
                z(i + 1, 1) = s;
            }
            // U z = y
            for (index_t i = n; i-- > 0;)
            {
                T s = z(i + 1, 1);
                for (I p = static_cast<I>(this->diagonal_[i] + 1); p < offsets[i + 1]; ++p)
                    s -= values[p] * z(static_cast<index_t>(columns[p]) + 1, 1);
                z(i + 1, 1) = s / values[this->diagonal_[i]];
            }
        }
    };

    template<class U, class I> ILU0Preconditioner(CSRMatrix<U, I> const&) -> ILU0Preconditioner<detail::solver_t<U>, I>;
}
#endif
