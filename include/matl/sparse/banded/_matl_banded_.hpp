/**************************************************************************************************
 * File: _matl_banded_.hpp
 *
 * Description: Direct solvers for structured sparse matrices: tridiagonal and banded systems.
 *
 *     x = solve_tridiagonal(lower, diagonal, upper, b)
 *                                 the Thomas algorithm: O(n), no pivoting (for diagonally
 *                                 dominant or symmetric positive definite matrices, e.g. 1D
 *                                 implicit time steps; otherwise use the banded LU)
 *
 *     BandedMatrix<T> B(n, kl, ku)    an n x n matrix with kl subdiagonals and ku superdiagonals;
 *                                     B(i,j) reads, B.entry(i,j) = v writes (inside the band)
 *     auto F = lu(B);                 LU with partial pivoting, O(n kl (kl + ku))
 *     auto C = cholesky(B);           Cholesky (symmetric positive definite; reads the lower
 *                                     band), O(n kl^2)
 *     x = F.solve(b);  F.rcond();  F.determinant();     as for the dense factorizations
 *     banded(S), banded(D)            the banded matrix of a sparse or dense matrix (bandwidths
 *                                     found from its nonzeros)
 *
 * A banded solver is a direct solver for sparse matrices whose nonzeros are near the diagonal:
 * the 5-point Laplacian on a k x k grid (numbered row by row) has bandwidth k, so its banded
 * Cholesky factorization costs O(n k^2) = O(n^2) for n = k^2 unknowns, instead of O(n^3) dense.
 * The fill-in stays inside the band.
 *
 * Storage: the band is kept column by column as in LAPACK (element (i,j) of the band of column j
 * at row ku + i - j of a (kl + ku + 1) x n array), so the arrays can be passed to LAPACK's banded
 * routines. Errors: a zero pivot is matl::singular_matrix_error (a pivot that is not positive in
 * Cholesky matl::not_positive_definite_error); dimensions and writes outside the band are always
 * checked (matl::dimension_error, matl::index_error).
 **************************************************************************************************/
#ifndef MATL_BANDED_HPP
#define MATL_BANDED_HPP
namespace matl
{
    /**************************************** Tridiagonal ****************************************/

   /***********************************************************************************************
    * solve_tridiagonal(lower, diagonal, upper, b)
    *
    * Solves the tridiagonal system with subdiagonal lower (n-1 x 1), diagonal (n x 1) and
    * superdiagonal upper (n-1 x 1), for the right-hand sides b (n x k), by the Thomas algorithm
    * (Gaussian elimination without pivoting). Stable for diagonally dominant and for symmetric
    * positive definite matrices; a zero pivot is reported as matl::singular_matrix_error.
    ***********************************************************************************************/
    template<MatrixAlgebraObject L, MatrixAlgebraObject D, MatrixAlgebraObject U, MatrixAlgebraObject B>
    auto solve_tridiagonal(L const& lower, D const& diagonal, U const& upper, B const& b)
    {
        typedef detail::solver_t<detail::promote_t<detail::dtype_t<L>, detail::dtype_t<D>, detail::dtype_t<U>, detail::dtype_t<B>>> T;
        const index_t n = diagonal.rows();
        if (diagonal.cols() != 1 || ((n > 0) && (lower.rows() != n - 1 || upper.rows() != n - 1 || lower.cols() != 1 || upper.cols() != 1)))
            detail::raise<dimension_error>("solve_tridiagonal: the diagonal must be n x 1 and the sub- and superdiagonals "
                                           "(n-1) x 1 (they are " + detail::dimensions(lower.rows(), lower.cols()) + ", " +
                                           detail::dimensions(diagonal.rows(), diagonal.cols()) + " and " +
                                           detail::dimensions(upper.rows(), upper.cols()) + ")");
        detail::check_right_hand_side("solve_tridiagonal", n, b.rows(), b.cols());

        const DynamicMatrix<T, StorageOrder::ColMajor> a(lower), d(diagonal), c(upper);
        DynamicMatrix<T, detail::solution_order_v<B>> x(b);
        if (n == 0) return x;
        std::vector<T> c_prime(n);                                  // modified superdiagonal

        // Forward elimination: row i becomes (0, 1, c'(i)) with right-hand side x(i,:)
        const auto pivot_check = [](T p, index_t i)
        {
            if (p == T(0))
                detail::raise<singular_matrix_error>("solve_tridiagonal: pivot " + std::to_string(i) +
                                                     " is 0 (the matrix needs pivoting; use lu of a BandedMatrix)");
        };
        T pivot = d(1);
        pivot_check(pivot, 1);
        for (index_t i = 1; i <= n; ++i)
        {
            if (i > 1)
            {
                pivot = d(i) - a(i - 1) * c_prime[i - 2];
                pivot_check(pivot, i);
                for (index_t k = 1; k <= x.cols(); ++k) x(i, k) = x(i, k) - a(i - 1) * x(i - 1, k);
            }
            if (i < n) c_prime[i - 1] = c(i) / pivot;
            for (index_t k = 1; k <= x.cols(); ++k) x(i, k) /= pivot;
        }
        // Back substitution
        for (index_t i = n - 1; i >= 1; --i)
            for (index_t k = 1; k <= x.cols(); ++k) x(i, k) -= c_prime[i - 1] * x(i + 1, k);
        return x;
    }


    /****************************************** Banded *******************************************/

   /***********************************************************************************************
    * class BandedMatrix<T>
    *
    * n x n matrix with kl subdiagonals and ku superdiagonals (see the description above).
    * B(i,j) reads an element (0 outside the band); B.entry(i,j) = v writes one inside the band.
    ***********************************************************************************************/
    template<class T>
    class BandedMatrix
    {
        static_assert(std::is_arithmetic_v<T>, "BandedMatrix: the element type must be an arithmetic type");

      public:
        typedef T dtype;

      private:
        index_t        n_;
        index_t        kl_;
        index_t        ku_;
        std::vector<T> band_;      // (kl + ku + 1) x n, column by column

        bool in_band(index_t i, index_t j) const {return (j <= i + this->ku_) && (i <= j + this->kl_);}

        index_t position(index_t i, index_t j) const {return (j - 1) * (this->kl_ + this->ku_ + 1) + (this->ku_ + i - j);}

      public:
        // n x n zero matrix with kl subdiagonals and ku superdiagonals
        BandedMatrix(index_t n, index_t kl, index_t ku) :
            n_(n), kl_(std::min(kl, n ? n - 1 : 0)), ku_(std::min(ku, n ? n - 1 : 0)),
            band_(n * (std::min(kl, n ? n - 1 : 0) + std::min(ku, n ? n - 1 : 0) + 1), T(0))
        {}

        index_t rows() const {return this->n_;}
        index_t cols() const {return this->n_;}
        index_t lower_bandwidth() const {return this->kl_;}
        index_t upper_bandwidth() const {return this->ku_;}

        // Element (i,j) (1-based): 0 outside the band
        T operator() (index_t i, index_t j) const
        {
            #ifdef MATL_BOUNDS_CHECK
            if ((i < 1) || (i > this->n_) || (j < 1) || (j > this->n_)) detail::index_out_of_range(i, j, this->n_, this->n_);
            #endif
            return this->in_band(i, j) ? this->band_[this->position(i, j)] : T(0);
        }

        // Element (i,j) to write (inside the band; matl::index_error outside it, always checked),
        // like CSRMatrix::entry
        T& entry(index_t i, index_t j)
        {
            if ((i < 1) || (i > this->n_) || (j < 1) || (j > this->n_)) detail::index_out_of_range(i, j, this->n_, this->n_);
            if (!this->in_band(i, j))
                detail::raise<index_error>("BandedMatrix: element (" + std::to_string(i) + "," + std::to_string(j) +
                                           ") is outside the band (" + std::to_string(this->kl_) + " subdiagonals, " +
                                           std::to_string(this->ku_) + " superdiagonals)");
            return this->band_[this->position(i, j)];
        }

        // The band in LAPACK's layout: (kl + ku + 1) x n, column by column
        std::span<const T> band() const {return this->band_;}
        std::span<T>       band()       {return this->band_;}
    };

    // The dense matrix of a banded matrix
    template<class T>
    DynamicMatrix<T, StorageOrder::RowMajor> dense(BandedMatrix<T> const& b)
    {
        DynamicMatrix<T, StorageOrder::RowMajor> d(b.rows(), b.cols());
        for (index_t i = 1; i <= b.rows(); ++i)
            for (index_t j = 1; j <= b.cols(); ++j) d(i, j) = b(i, j);
        return d;
    }

    // The banded matrix of a square sparse matrix (the bandwidths of its stored entries)
    template<class T, class I>
    BandedMatrix<T> banded(CSRMatrix<T, I> const& s)
    {
        detail::check_square("banded", s.rows(), s.cols());
        const auto offsets = s.row_offsets(), columns = s.column_indices();
        const auto values = s.values();
        index_t kl = 0, ku = 0;
        for (index_t i = 1; i <= s.rows(); ++i)
            for (I p = offsets[i - 1]; p < offsets[i]; ++p)
            {
                const index_t j = static_cast<index_t>(columns[p]) + 1;
                if (i > j) kl = std::max(kl, i - j); else ku = std::max(ku, j - i);
            }
        BandedMatrix<T> b(s.rows(), kl, ku);
        for (index_t i = 1; i <= s.rows(); ++i)
            for (I p = offsets[i - 1]; p < offsets[i]; ++p) b.entry(i, static_cast<index_t>(columns[p]) + 1) = values[p];
        return b;
    }

    // The banded matrix of a square dense matrix, view or expression (the bandwidths of its nonzeros)
    template<MatrixAlgebraObject X>
    BandedMatrix<detail::dtype_t<X>> banded(X const& x)
    {return banded(sparse(x));}

    // Banded matrix times a dense matrix, view or expression (a dense result, computed now)
    template<class T, MatrixAlgebraObject X>
    DynamicMatrix<detail::promote_t<T, detail::dtype_t<X>>, StorageOrder::RowMajor> operator* (BandedMatrix<T> const& b, X const& x_in)
    {
        typedef detail::promote_t<T, detail::dtype_t<X>> R;
        if (b.cols() != x_in.rows())
            detail::incompatible_dimensions("banded matrix product", b.rows(), b.cols(), x_in.rows(), x_in.cols());
        const DynamicMatrix<R, StorageOrder::RowMajor> x(x_in);
        DynamicMatrix<R, StorageOrder::RowMajor> y(b.rows(), x.cols());
        y = R(0);
        const index_t kl = b.lower_bandwidth(), ku = b.upper_bandwidth();
        for (index_t i = 1; i <= b.rows(); ++i)
            for (index_t j = (i > kl ? i - kl : 1); j <= std::min(b.cols(), i + ku); ++j)
            {
                const R v = R(b(i, j));
                for (index_t c = 1; c <= x.cols(); ++c) y(i, c) += v * x(j, c);
            }
        return y;
    }


   /***********************************************************************************************
    * class BandedLUFactorization<T>
    *
    * LU factorization with partial pivoting of a banded matrix (LAPACK's xGBTRF, unblocked): the
    * row interchanges widen U to kl + ku superdiagonals, so the factors are kept in a band with kl
    * more rows. Made by lu(B).
    ***********************************************************************************************/
    template<class T>
    class BandedLUFactorization
    {
        static_assert(std::is_floating_point_v<T>, "BandedLUFactorization: the element type must be a floating-point type");

      private:
        index_t              n_, kl_, ku_, ld_;    // ld_ = 2 kl + ku + 1 rows of the band
        std::vector<T>       lu_;                  // U (kl + ku superdiagonals) and the multipliers of L
        std::vector<index_t> pivots_;
        T                    norm_1_;

        T&       at(index_t i, index_t j)       {return this->lu_[(j - 1) * this->ld_ + (this->kl_ + this->ku_ + i - j)];}
        T const& at(index_t i, index_t j) const {return this->lu_[(j - 1) * this->ld_ + (this->kl_ + this->ku_ + i - j)];}

        // x = A^-1 x (transposed: A^-T x), in place, for an n x k matrix x
        template<typename X>
        void solve_columns(X& x, bool transposed) const
        {
            const index_t n = this->n_, kl = this->kl_, kv = this->kl_ + this->ku_;
            if (n == 0) return;                                    // (the loops below count down from n - 1)
            for (index_t c = 1; c <= x.cols(); ++c)
            {
                if (!transposed)
                {
                    // L: interchanges and multipliers, step by step
                    for (index_t j = 1; j < n; ++j)
                    {
                        const index_t p = this->pivots_[j - 1];
                        if (p != j) std::swap(x(j, c), x(p, c));
                        for (index_t i = j + 1; i <= std::min(n, j + kl); ++i) x(i, c) -= this->at(i, j) * x(j, c);
                    }
                    // U: back substitution within kl + ku superdiagonals
                    for (index_t i = n; i >= 1; --i)
                    {
                        T s = x(i, c);
                        for (index_t j = i + 1; j <= std::min(n, i + kv); ++j) s -= this->at(i, j) * x(j, c);
                        x(i, c) = s / this->at(i, i);
                    }
                }
                else
                {
                    // U^T: forward substitution
                    for (index_t i = 1; i <= n; ++i)
                    {
                        T s = x(i, c);
                        for (index_t j = (i > kv ? i - kv : 1); j < i; ++j) s -= this->at(j, i) * x(j, c);
                        x(i, c) = s / this->at(i, i);
                    }
                    // L^T: multipliers and interchanges, in reverse
                    for (index_t j = n - 1; j >= 1; --j)
                    {
                        T s = x(j, c);
                        for (index_t i = j + 1; i <= std::min(n, j + kl); ++i) s -= this->at(i, j) * x(i, c);
                        x(j, c) = s;
                        const index_t p = this->pivots_[j - 1];
                        if (p != j) std::swap(x(j, c), x(p, c));
                    }
                }
            }
        }

      public:
        explicit BandedLUFactorization(BandedMatrix<T> const& b) :
            n_(b.rows()), kl_(b.lower_bandwidth()), ku_(b.upper_bandwidth()), ld_(2 * b.lower_bandwidth() + b.upper_bandwidth() + 1),
            lu_(b.rows() * (2 * b.lower_bandwidth() + b.upper_bandwidth() + 1), T(0)), pivots_(b.rows()), norm_1_(0)
        {
            const index_t n = this->n_, kl = this->kl_, ku = this->ku_;
            for (index_t j = 1; j <= n; ++j)
            {
                T column = 0;
                for (index_t i = (j > ku ? j - ku : 1); i <= std::min(n, j + kl); ++i)
                {
                    this->at(i, j) = b(i, j);
                    column += std::abs(b(i, j));
                }
                this->norm_1_ = std::max(this->norm_1_, column);
            }

            index_t ju = 1;                                        // last column U reaches so far
            for (index_t j = 1; j <= n; ++j)
            {
                const index_t last = std::min(n, j + kl);          // last row with a nonzero in column j
                index_t p = j;
                for (index_t i = j + 1; i <= last; ++i)
                    if (std::abs(this->at(i, j)) > std::abs(this->at(p, j))) p = i;
                this->pivots_[j - 1] = p;
                if (this->at(p, j) == T(0))
                    detail::raise<singular_matrix_error>("lu (banded): the matrix is singular (column " + std::to_string(j) +
                                                         " has no nonzero pivot)");
                ju = std::max(ju, std::min(n, p + ku));
                if (p != j)
                    for (index_t c = j; c <= ju; ++c) std::swap(this->at(j, c), this->at(p, c));
                const T reciprocal = T(1) / this->at(j, j);
                for (index_t i = j + 1; i <= last; ++i) this->at(i, j) *= reciprocal;
                for (index_t c = j + 1; c <= ju; ++c)
                {
                    const T u = this->at(j, c);
                    if (u != T(0))
                        for (index_t i = j + 1; i <= last; ++i) this->at(i, c) -= this->at(i, j) * u;
                }
            }
        }

        index_t size() const {return this->n_;}

        // Row interchanges: row k was interchanged with row pivots()[k-1] at step k (1-based)
        std::vector<index_t> const& pivots() const {return this->pivots_;}

        // Overwrites x (a matrix or view holding B, one system per column) with the solution of A X = B
        template<typename X> requires WritableMatrix<X>
        void solve_inplace(X&& x) const
        {
            static_assert(std::is_floating_point_v<detail::dtype_t<X>>, "solve_inplace: the solution must have a floating-point type");
            detail::check_right_hand_side("lu solve (banded)", this->n_, x.rows(), x.cols());
            this->solve_columns(x, false);
        }

        // Solution of A X = B
        template<MatrixAlgebraObject B>
        auto solve(B const& b) const
        {
            DynamicMatrix<detail::solver_t<detail::promote_t<T, detail::dtype_t<B>>>, detail::solution_order_v<B>> x(b);
            this->solve_inplace(x);
            return x;
        }

        T determinant() const
        {
            T d = 1;
            for (index_t k = 1; k <= this->n_; ++k)
            {
                d *= this->at(k, k);
                if (this->pivots_[k - 1] != k) d = -d;
            }
            return d;
        }

        // Estimated reciprocal condition number in the 1-norm (see LUFactorization::rcond)
        T rcond() const
        {
            const auto solve = [this](auto& x) {this->solve_columns(x, false);};
            const auto solve_transposed = [this](auto& x) {this->solve_columns(x, true);};
            return detail::reciprocal_condition(this->norm_1_, detail::estimate_inverse_norm_1<T>(this->n_, solve, solve_transposed), this->n_);
        }
    };


   /***********************************************************************************************
    * class BandedCholeskyFactorization<T>
    *
    * Cholesky factorization A = L L^T of a symmetric positive definite banded matrix (LAPACK's
    * xPBTRF, unblocked), reading the lower band of A; L has A's kl subdiagonals. Made by
    * cholesky(B).
    ***********************************************************************************************/
    template<class T>
    class BandedCholeskyFactorization
    {
        static_assert(std::is_floating_point_v<T>, "BandedCholeskyFactorization: the element type must be a floating-point type");

      private:
        index_t        n_, kd_;
        std::vector<T> l_;          // (kd + 1) x n: column j holds L(j, j) .. L(j + kd, j)
        T              norm_1_;

        T&       at(index_t i, index_t j)       {return this->l_[(j - 1) * (this->kd_ + 1) + (i - j)];}
        T const& at(index_t i, index_t j) const {return this->l_[(j - 1) * (this->kd_ + 1) + (i - j)];}

      public:
        explicit BandedCholeskyFactorization(BandedMatrix<T> const& b) :
            n_(b.rows()), kd_(b.lower_bandwidth()), l_(b.rows() * (b.lower_bandwidth() + 1), T(0)), norm_1_(0)
        {
            const index_t n = this->n_, kd = this->kd_;
            for (index_t j = 1; j <= n; ++j)
                for (index_t i = j; i <= std::min(n, j + kd); ++i) this->at(i, j) = b(i, j);

            // 1-norm of the symmetric matrix from its lower band
            for (index_t j = 1; j <= n; ++j)
            {
                T column = 0;
                for (index_t i = (j > kd ? j - kd : 1); i < j; ++i) column += std::abs(this->at(j, i));
                for (index_t i = j; i <= std::min(n, j + kd); ++i) column += std::abs(this->at(i, j));
                this->norm_1_ = std::max(this->norm_1_, column);
            }

            for (index_t j = 1; j <= n; ++j)
            {
                const index_t first = (j > kd ? j - kd : 1);
                T d = this->at(j, j);
                for (index_t p = first; p < j; ++p) d -= this->at(j, p) * this->at(j, p);
                if (!(d > T(0)))
                    detail::raise<not_positive_definite_error>("cholesky (banded): the matrix is not positive definite (pivot " +
                                                               std::to_string(j) + " is not positive)");
                const T l_jj = std::sqrt(d);
                this->at(j, j) = l_jj;
                for (index_t i = j + 1; i <= std::min(n, j + kd); ++i)
                {
                    T s = this->at(i, j);
                    for (index_t p = (i > kd ? i - kd : 1); p < j; ++p) s -= this->at(i, p) * this->at(j, p);
                    this->at(i, j) = s / l_jj;
                }
            }
        }

        index_t size() const {return this->n_;}

        // Overwrites x (holding B) with the solution of A X = B: L Y = B, then L^T X = Y
        template<typename X> requires WritableMatrix<X>
        void solve_inplace(X&& x) const
        {
            static_assert(std::is_floating_point_v<detail::dtype_t<X>>, "solve_inplace: the solution must have a floating-point type");
            detail::check_right_hand_side("cholesky solve (banded)", this->n_, x.rows(), x.cols());
            const index_t n = this->n_, kd = this->kd_;
            for (index_t c = 1; c <= x.cols(); ++c)
            {
                for (index_t i = 1; i <= n; ++i)
                {
                    T s = x(i, c);
                    for (index_t p = (i > kd ? i - kd : 1); p < i; ++p) s -= this->at(i, p) * x(p, c);
                    x(i, c) = s / this->at(i, i);
                }
                for (index_t i = n; i >= 1; --i)
                {
                    T s = x(i, c);
                    for (index_t p = i + 1; p <= std::min(n, i + kd); ++p) s -= this->at(p, i) * x(p, c);
                    x(i, c) = s / this->at(i, i);
                }
            }
        }

        template<MatrixAlgebraObject B>
        auto solve(B const& b) const
        {
            DynamicMatrix<detail::solver_t<detail::promote_t<T, detail::dtype_t<B>>>, detail::solution_order_v<B>> x(b);
            this->solve_inplace(x);
            return x;
        }

        T determinant() const
        {
            T d = 1;
            for (index_t k = 1; k <= this->n_; ++k) d *= this->at(k, k);
            return d * d;
        }

        T rcond() const
        {
            const auto solve = [this](auto& x) {this->solve_inplace(x);};
            return detail::reciprocal_condition(this->norm_1_, detail::estimate_inverse_norm_1<T>(this->n_, solve, solve), this->n_);
        }
    };


    // LU factorization (partial pivoting) of a banded matrix, in a floating-point type
    template<class T>
    auto lu(BandedMatrix<T> const& b)
    {
        typedef detail::solver_t<T> F;
        if constexpr (std::is_same_v<F, T>) return BandedLUFactorization<F>(b);
        else
        {
            BandedMatrix<F> f(b.rows(), b.lower_bandwidth(), b.upper_bandwidth());
            for (index_t k = 0; k < b.band().size(); ++k) f.band()[k] = static_cast<F>(b.band()[k]);
            return BandedLUFactorization<F>(f);
        }
    }

    // Cholesky factorization of a symmetric positive definite banded matrix (its lower band is read)
    template<class T>
    auto cholesky(BandedMatrix<T> const& b)
    {
        typedef detail::solver_t<T> F;
        if constexpr (std::is_same_v<F, T>) return BandedCholeskyFactorization<F>(b);
        else
        {
            BandedMatrix<F> f(b.rows(), b.lower_bandwidth(), b.upper_bandwidth());
            for (index_t k = 0; k < b.band().size(); ++k) f.band()[k] = static_cast<F>(b.band()[k]);
            return BandedCholeskyFactorization<F>(f);
        }
    }

    template<class T, MatrixAlgebraObject B>
    auto lu_solve(BandedMatrix<T> const& a, B const& b) {return lu(a).solve(b);}

    template<class T, MatrixAlgebraObject B>
    auto cholesky_solve(BandedMatrix<T> const& a, B const& b) {return cholesky(a).solve(b);}
}
#endif
