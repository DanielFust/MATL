/**************************************************************************************************
 * File: _matl_lu_.hpp
 *
 * Description: LU factorization with partial pivoting, P A = L U, for square matrices.
 *
 *     auto F = lu(A);                 factorizes a copy of A (A: matrix, view or expression)
 *     auto F = lu_inplace(A);         factorizes A in its own memory (A: a matrix of a
 *                                     floating-point type), which then holds the factors
 *     X = F.solve(B);                 solution of A X = B (B: one system per column)
 *     F.solve_inplace(X);             overwrites X (holding B) with the solution
 *     X = lu_solve(A, B);             factorizes and solves once
 *
 *     F.rcond()                       estimated reciprocal condition number (1-norm)
 *     F.determinant()                 determinant of A
 *     F.L(), F.U(), F.P()             the factors as separate matrices (P A = L U)
 *     F.factors(), F.pivots()         the packed factors and the row interchanges
 *
 * The factors are packed in one matrix, as in LAPACK: U on and above the diagonal, L below it
 * (L's diagonal is 1 and is not stored). Row k was interchanged with row pivots()[k-1] at step k
 * (both 1-based).
 *
 * An exactly zero pivot (a column with no nonzero element left to pivot on) means the matrix is
 * singular, and is reported as matl::singular_matrix_error. A matrix that is only close to
 * singular is factorized; rcond() then tells how close (see CONDITION below).
 *
 * Algorithm: right-looking, blocked. Each panel of BLOCK columns is factorized column by column
 * (pivot search, row interchange within the panel, elimination within the panel), the panel's
 * row interchanges are applied to the other columns, the rows of U to its right are solved for,
 * and the rest of the matrix is updated with one matrix product,
 * A22 -= A21 * U12, which runs on the cache-blocked product kernel through views. For large
 * matrices this product is almost all the work, so the factorization runs close to the speed of
 * the product.
 *
 * CONDITION: rcond() estimates 1 / (||A||_1 ||A^-1||_1) with a few triangular solves (see
 * detail::estimate_inverse_norm_1), not with the inverse. It is close to 1 for a well-conditioned
 * matrix and close to 0 for a nearly singular one: a solution then loses about
 * log10(1 / rcond) of the type's significant digits (about 16 for double).
 **************************************************************************************************/
#ifndef MATL_LU_HPP
#define MATL_LU_HPP
namespace matl
{
    namespace detail
    {
        // Applies the row interchanges of steps k0..k1 to columns c0..c1 of a. Row-major matrices
        // swap row segments; column-major ones apply all the interchanges to one column at a
        // time, so that each column is read once from memory (as LAPACK's xLASWP).
        template<class T, StorageOrder ORDER>
        void interchange_rows(DynamicMatrix<T, ORDER>& a, std::vector<index_t> const& pivots,
                              index_t k0, index_t k1, index_t c0, index_t c1)
        {
            if constexpr (ORDER == StorageOrder::RowMajor)
            {
                for (index_t k = k0; k <= k1; ++k)
                    if (pivots[k - 1] != k)
                        for (index_t c = c0; c <= c1; ++c) std::swap(a(k,c), a(pivots[k - 1],c));
            }
            else
            {
                for (index_t c = c0; c <= c1; ++c)
                    for (index_t k = k0; k <= k1; ++k)
                        if (pivots[k - 1] != k) std::swap(a(k,c), a(pivots[k - 1],c));
            }
        }

       /*******************************************************************************************
        * void lu_factorize(DynamicMatrix<T, ORDER>& a, std::vector<index_t>& pivots)
        *
        * Overwrites the square matrix a with its LU factors (see the description above) and
        * fills in the row interchanges. Reports matl::singular_matrix_error for a zero pivot.
        *******************************************************************************************/
        template<class T, StorageOrder ORDER>
        void lu_factorize(DynamicMatrix<T, ORDER>& a, std::vector<index_t>& pivots)
        {
            const index_t n = a.rows();
            pivots.resize(n);

            for (index_t k0 = 1; k0 <= n; k0 += FACTORIZATION_BLOCK)
            {
                const index_t k1 = std::min(k0 + FACTORIZATION_BLOCK - 1, n);    // last column of the panel

                // Panel: columns k0..k1, rows k0..n
                for (index_t j = k0; j <= k1; ++j)
                {
                    // Pivot: the largest element (in magnitude) of column j on or below the diagonal
                    index_t p = j;
                    T largest = std::abs(a(j,j));
                    for (index_t i = j + 1; i <= n; ++i)
                        if (std::abs(a(i,j)) > largest) {largest = std::abs(a(i,j)); p = i;}
                    pivots[j - 1] = p;
                    if (largest == T(0))
                        raise<singular_matrix_error>("lu: the matrix is singular (column " + std::to_string(j) +
                                                     " has no nonzero pivot)");

                    // Interchange rows j and p within the panel (the other columns follow below)
                    if (p != j)
                        for (index_t c = k0; c <= k1; ++c) std::swap(a(j,c), a(p,c));

                    // Multipliers: column j of L
                    const T reciprocal = T(1) / a(j,j);
                    for (index_t i = j + 1; i <= n; ++i) a(i,j) *= reciprocal;

                    // Eliminate within the panel: a(i,c) -= l(i,j) * u(j,c) for the panel's columns
                    // right of j (in the storage order of a)
                    if constexpr (ORDER == StorageOrder::RowMajor)
                    {
                        for (index_t i = j + 1; i <= n; ++i)
                        {
                            const T l = a(i,j);
                            for (index_t c = j + 1; c <= k1; ++c) a(i,c) -= l * a(j,c);
                        }
                    }
                    else
                    {
                        for (index_t c = j + 1; c <= k1; ++c)
                        {
                            const T u = a(j,c);
                            for (index_t i = j + 1; i <= n; ++i) a(i,c) -= a(i,j) * u;
                        }
                    }
                }

                // The panel's row interchanges, in the columns left and right of it
                interchange_rows(a, pivots, k0, k1, 1, k0 - 1);
                interchange_rows(a, pivots, k0, k1, k1 + 1, n);

                if (k1 == n) break;

                // U12 = L11^-1 A12: rows k0..k1 of U right of the panel (L11 is unit lower
                // triangular), along the rows or down the columns of A12
                if constexpr (ORDER == StorageOrder::RowMajor)
                {
                    for (index_t i = k0 + 1; i <= k1; ++i)
                        for (index_t p = k0; p < i; ++p)
                        {
                            const T l = a(i,p);
                            for (index_t c = k1 + 1; c <= n; ++c) a(i,c) -= l * a(p,c);
                        }
                }
                else
                {
                    for (index_t c = k1 + 1; c <= n; ++c)
                        for (index_t i = k0 + 1; i <= k1; ++i)
                        {
                            T s = a(i,c);
                            for (index_t p = k0; p < i; ++p) s -= a(i,p) * a(p,c);
                            a(i,c) = s;
                        }
                }

                // A22 -= L21 * U12 (the product kernel, through views)
                a(range(k1 + 1, n), range(k1 + 1, n)) -= a(range(k1 + 1, n), range(k0, k1)) *
                                                          a(range(k0, k1), range(k1 + 1, n));
            }
        }
    }


   /***********************************************************************************************
    * class LUFactorization<T, ORDER>
    *
    * LU factorization with partial pivoting, P A = L U (see the description above). Made by lu
    * and lu_inplace.
    *
    * Template Parameters:
    *     T         Floating-point element type
    *     ORDER     Storage order of the factors
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    class LUFactorization
    {
        static_assert(std::is_floating_point_v<T>, "LUFactorization: the element type must be a floating-point type");

      private:
        DynamicMatrix<T, ORDER> lu_;        // U on and above the diagonal, L below it
        std::vector<index_t>    pivots_;    // row k was interchanged with row pivots_[k-1]
        T                       norm_1_;    // 1-norm of A (for rcond)

        // Applies the row interchanges to x (P x), or undoes them (P^T x)
        template<typename X>
        void permute(X& x, bool inverse) const
        {
            const index_t n = this->size();
            const auto interchange = [&x](index_t k, index_t p)
            {
                if (p != k)
                    for (index_t j = 1; j <= x.cols(); ++j) std::swap(x(k,j), x(p,j));
            };
            if (!inverse) for (index_t k = 1; k <= n; ++k) interchange(k, this->pivots_[k - 1]);
            else          for (index_t k = n; k >= 1; --k) interchange(k, this->pivots_[k - 1]);
        }

      public:
        // Factorizes the square matrix a (its memory may be a matrix's, see lu_inplace)
        explicit LUFactorization(DynamicMatrix<T, ORDER> a) : lu_(std::move(a)), pivots_(), norm_1_(0)
        {
            detail::check_square("lu", this->lu_.rows(), this->lu_.cols());
            this->norm_1_ = norm_1(this->lu_);
            detail::lu_factorize(this->lu_, this->pivots_);
        }

        // Order of the system (A is size() x size())
        index_t size() const {return this->lu_.rows();}

        // Packed factors: U on and above the diagonal, L (unit diagonal not stored) below it
        DynamicMatrix<T, ORDER> const& factors() const {return this->lu_;}

        // Row interchanges: row k was interchanged with row pivots()[k-1] at step k (1-based)
        std::vector<index_t> const& pivots() const {return this->pivots_;}

        // Unit lower triangular factor L
        DynamicMatrix<T, ORDER> L() const
        {
            const index_t n = this->size();
            DynamicMatrix<T, ORDER> l(n, n);
            for (index_t i = 1; i <= n; ++i)
                for (index_t j = 1; j <= n; ++j)
                    l(i,j) = (i > j) ? this->lu_(i,j) : (i == j) ? T(1) : T(0);
            return l;
        }

        // Upper triangular factor U
        DynamicMatrix<T, ORDER> U() const
        {
            const index_t n = this->size();
            DynamicMatrix<T, ORDER> u(n, n);
            for (index_t i = 1; i <= n; ++i)
                for (index_t j = 1; j <= n; ++j)
                    u(i,j) = (i <= j) ? this->lu_(i,j) : T(0);
            return u;
        }

        // Permutation matrix P, with P A = L U
        DynamicMatrix<T, ORDER> P() const
        {
            DynamicMatrix<T, ORDER> p = eye<T>(this->size());
            this->permute(p, false);
            return p;
        }

        // Overwrites x (a matrix or view holding B, one system per column) with the solution of
        // A X = B
        template<typename X> requires WritableMatrix<X>
        void solve_inplace(X&& x) const
        {
            static_assert(std::is_floating_point_v<detail::dtype_t<X>>,
                          "LUFactorization::solve_inplace: the solution must have a floating-point type");
            detail::check_right_hand_side("lu solve", this->size(), x.rows(), x.cols());
            this->permute(x, false);
            detail::triangular_solve<true, true>(this->lu_, this->size(), x, "lu solve");
            detail::triangular_solve<false, false>(this->lu_, this->size(), x, "lu solve");
        }

        // Solution of A X = B (B: matrix, view or expression; one system per column)
        template<MatrixAlgebraObject B>
        auto solve(B const& b) const
        {
            DynamicMatrix<detail::solver_t<detail::promote_t<T, detail::dtype_t<B>>>, detail::solution_order_v<B>> x(b);
            this->solve_inplace(x);
            return x;
        }

        // Determinant of A: the product of U's diagonal, with the sign of the permutation
        T determinant() const
        {
            T d = 1;
            for (index_t k = 1; k <= this->size(); ++k)
            {
                d *= this->lu_(k,k);
                if (this->pivots_[k - 1] != k) d = -d;
            }
            return d;
        }

        // Estimated reciprocal condition number in the 1-norm, 1 / (||A||_1 ||A^-1||_1): near 1
        // for a well-conditioned A, near 0 (or 0) for a nearly singular one. Costs a few
        // triangular solves (O(n^2)).
        T rcond() const
        {
            const index_t n = this->size();
            // A^-1 x = U^-1 L^-1 P x, and A^-T x = P^T L^-T U^-T x
            const auto solve = [this](auto& x) {this->solve_inplace(x);};
            const auto solve_transposed = [this, n](auto& x)
            {
                const auto ut = [this](index_t i, index_t p) {return this->lu_(p,i);};
                detail::triangular_solve<true, false>(ut, n, x, "lu rcond");
                detail::triangular_solve<false, true>(ut, n, x, "lu rcond");
                this->permute(x, true);
            };
            return detail::reciprocal_condition(this->norm_1_,
                                                detail::estimate_inverse_norm_1<T>(n, solve, solve_transposed), n);
        }
    };


    // LU factorization of a copy of a (a square matrix, view or expression), in a floating-point
    // type (double for integer matrices) and the storage order of a
    template<MatrixAlgebraObject A>
    auto lu(A const& a)
    {
        typedef detail::solver_t<detail::dtype_t<A>> T;
        constexpr StorageOrder order = detail::solution_order_v<A>;
        return LUFactorization<T, order>(DynamicMatrix<T, order>(a));
    }

    // LU factorization in the memory of a (a square matrix of a floating-point type): a then holds
    // the packed factors, and must outlive the factorization
    template<class T, StorageOrder ORDER>
    LUFactorization<T, ORDER> lu_inplace(DynamicMatrix<T, ORDER>& a)
    {return LUFactorization<T, ORDER>(DynamicMatrix<T, ORDER>(a.rows(), a.cols(), a.data()));}

    template<class T, index_t M, index_t N, StorageOrder ORDER>
    LUFactorization<T, ORDER> lu_inplace(FixedSizeMatrix<T, M, N, ORDER>& a)
    {return LUFactorization<T, ORDER>(DynamicMatrix<T, ORDER>(M, N, a.data()));}

    // Solution of A X = B by LU factorization (A square; B: one system per column)
    template<MatrixAlgebraObject A, MatrixAlgebraObject B>
    auto lu_solve(A const& a, B const& b)
    {return lu(a).solve(b);}
}
#endif
