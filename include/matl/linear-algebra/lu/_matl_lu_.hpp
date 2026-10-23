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
 * Algorithm: recursive (Toledo). The columns are split in halves: the left half is factorized
 * (recursively), its row interchanges are applied to the right half, the right half's rows of U
 * are solved for (a recursive triangular solve), the rest of the right half is updated with one
 * matrix product, A22 -= L21 * U12, and then factorized (recursively). Splitting all the way down
 * to RECURSION_BASE columns turns almost all the arithmetic into matrix products of every size,
 * which run on the cache-blocked product kernel through views.
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
        template<class MAT, class PIVOTS>
        void interchange_rows(MAT& a, PIVOTS const& pivots, index_t k0, index_t k1, index_t c0, index_t c1)
        {
            if constexpr (traits<MAT>::storage_order == StorageOrder::RowMajor)
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
        * void lu_columns(a, pivots, j0, j1)
        *
        * Factorizes columns j0..j1 of a (rows j0..n) column by column: pivot search, row
        * interchange within these columns, multipliers, and elimination within these columns.
        * The base case of lu_recursive.
        *******************************************************************************************/
        template<class MAT, class PIVOTS>
        void lu_columns(MAT& a, PIVOTS& pivots, index_t j0, index_t j1)
        {
            typedef dtype_t<MAT> T;
            constexpr StorageOrder ORDER = traits<MAT>::storage_order;
            const index_t n = a.rows();
            for (index_t j = j0; j <= j1; ++j)
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

                // Interchange rows j and p within these columns (lu_recursive does the others)
                if (p != j)
                    for (index_t c = j0; c <= j1; ++c) std::swap(a(j,c), a(p,c));

                // Multipliers: column j of L
                const T reciprocal = T(1) / a(j,j);
                for (index_t i = j + 1; i <= n; ++i) a(i,j) *= reciprocal;

                // Eliminate: a(i,c) -= l(i,j) * u(j,c) for the columns right of j (in the storage
                // order of a)
                if constexpr (ORDER == StorageOrder::RowMajor)
                {
                    for (index_t i = j + 1; i <= n; ++i)
                    {
                        const T l = a(i,j);
                        for (index_t c = j + 1; c <= j1; ++c) a(i,c) -= l * a(j,c);
                    }
                }
                else
                {
                    for (index_t c = j + 1; c <= j1; ++c)
                    {
                        const T u = a(j,c);
                        for (index_t i = j + 1; i <= n; ++i) a(i,c) -= a(i,j) * u;
                    }
                }
            }
        }

       /*******************************************************************************************
        * void solve_unit_lower(a, l0, l1, c0, c1)
        *
        * Solves L X = B in place, for the unit lower triangular L = a(l0:l1, l0:l1) and
        * B = a(l0:l1, c0:c1) (rows of U right of a factorized block). L is split in halves:
        * X1 = L11^-1 B1, B2 -= L21 X1 (a matrix product), X2 = L22^-1 B2, so that most of the work
        * runs on the product kernel.
        *******************************************************************************************/
        template<class MAT>
        void solve_unit_lower(MAT& a, index_t l0, index_t l1, index_t c0, index_t c1)
        {
            typedef dtype_t<MAT> T;
            constexpr StorageOrder ORDER = traits<MAT>::storage_order;
            if (c0 > c1 || l0 >= l1) return;                  // nothing to solve, or L = 1
            const index_t k = l1 - l0 + 1;
            if (k <= RECURSION_BASE)
            {
                if constexpr (ORDER == StorageOrder::RowMajor)
                {
                    for (index_t i = l0 + 1; i <= l1; ++i)
                        for (index_t p = l0; p < i; ++p)
                        {
                            const T l = a(i,p);
                            for (index_t c = c0; c <= c1; ++c) a(i,c) -= l * a(p,c);
                        }
                }
                else
                {
                    for (index_t c = c0; c <= c1; ++c)
                        for (index_t i = l0 + 1; i <= l1; ++i)
                        {
                            T s = a(i,c);
                            for (index_t p = l0; p < i; ++p) s -= a(i,p) * a(p,c);
                            a(i,c) = s;
                        }
                }
                return;
            }
            const index_t lm = l0 + k / 2;                    // first row of the second half
            solve_unit_lower(a, l0, lm - 1, c0, c1);
            auto B2 = a(range(lm, l1), range(c0, c1));
            subtract_product(B2, a(range(lm, l1), range(l0, lm - 1)), a(range(l0, lm - 1), range(c0, c1)));
            solve_unit_lower(a, lm, l1, c0, c1);
        }

       /*******************************************************************************************
        * void lu_recursive(a, pivots, j0, j1)
        *
        * Factorizes columns j0..j1 of a (rows j0..n), applying their row interchanges within those
        * columns (Toledo's recursive LU): factorize the left half of the columns, apply its
        * interchanges to the right half, solve for the right half's rows of U, update the rest of
        * the right half with one matrix product, factorize it, and apply its interchanges to the
        * left half. Narrow column ranges are factorized column by column (lu_columns).
        *******************************************************************************************/
        template<class MAT, class PIVOTS>
        void lu_recursive(MAT& a, PIVOTS& pivots, index_t j0, index_t j1)
        {
            const index_t n = a.rows();
            if (j1 - j0 + 1 <= RECURSION_BASE)
            {
                lu_columns(a, pivots, j0, j1);
                return;
            }
            const index_t jm = j0 + (j1 - j0 + 1) / 2;       // first column of the right half

            lu_recursive(a, pivots, j0, jm - 1);
            interchange_rows(a, pivots, j0, jm - 1, jm, j1);
            solve_unit_lower(a, j0, jm - 1, jm, j1);         // U12 = L11^-1 A12
            auto A22 = a(range(jm, n), range(jm, j1));       // A22 -= L21 * U12
            subtract_product(A22, a(range(jm, n), range(j0, jm - 1)), a(range(j0, jm - 1), range(jm, j1)));
            lu_recursive(a, pivots, jm, j1);
            interchange_rows(a, pivots, jm, j1, j0, jm - 1);
        }

       /*******************************************************************************************
        * void lu_factorize(MAT& a, PIVOTS& pivots)
        *
        * Overwrites the square matrix a (dynamic or fixed-size) with its LU factors (see the
        * description above) and fills in the row interchanges (a std::vector, or a std::array of
        * a fixed size). Reports matl::singular_matrix_error for a zero pivot.
        *******************************************************************************************/
        template<class MAT, class PIVOTS>
        void lu_factorize(MAT& a, PIVOTS& pivots)
        {
            if constexpr (requires {pivots.resize(a.rows());}) pivots.resize(a.rows());
            if (a.rows() > 0) lu_recursive(a, pivots, 1, a.rows());
        }
    }


   /***********************************************************************************************
    * class BasicLUFactorization<MAT>
    *
    * LU factorization with partial pivoting, P A = L U (see the description above), with the
    * factors in a matrix of type MAT. Made by lu and lu_inplace, and usually named through:
    *
    *     LUFactorization<T, ORDER>          dynamic factors (any n)
    *     FixedLUFactorization<T, N, ORDER>  N x N fixed-size factors and a std::array of pivots:
    *                                        no heap memory (lu of a fixed-size matrix, e.g. an
    *                                        element matrix); solving a fixed-size right-hand side
    *                                        gives a fixed-size solution
    *
    * Template Parameters:
    *     MAT       Matrix of the factors: DynamicMatrix<T, ORDER> or FixedSizeMatrix<T, N, N, ORDER>
    *               (T a floating-point type)
    ***********************************************************************************************/
    template<class MAT>
    class BasicLUFactorization
    {
        typedef detail::factor_storage<MAT> storage;
        typedef typename storage::dtype     T;
        static_assert(std::is_floating_point_v<T>, "LUFactorization: the element type must be a floating-point type");

      public:
        typedef T dtype;

      private:
        MAT                         lu_;        // U on and above the diagonal, L below it
        typename storage::pivots    pivots_;    // row k was interchanged with row pivots_[k-1]
        T                           norm_1_;    // 1-norm of A (for rcond)

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
        explicit BasicLUFactorization(MAT a) : lu_(std::move(a)), pivots_(), norm_1_(0)
        {
            detail::check_square("lu", this->lu_.rows(), this->lu_.cols());
            this->norm_1_ = norm_1(this->lu_);
            detail::lu_factorize(this->lu_, this->pivots_);
        }

        // Order of the system (A is size() x size())
        index_t size() const {return this->lu_.rows();}

        // Packed factors: U on and above the diagonal, L (unit diagonal not stored) below it
        MAT const& factors() const {return this->lu_;}

        // Row interchanges: row k was interchanged with row pivots()[k-1] at step k (1-based)
        typename storage::pivots const& pivots() const {return this->pivots_;}

        // Unit lower triangular factor L
        MAT L() const
        {
            const index_t n = this->size();
            MAT l = storage::square(n);
            for (index_t i = 1; i <= n; ++i)
                for (index_t j = 1; j <= n; ++j)
                    l(i,j) = (i > j) ? this->lu_(i,j) : (i == j) ? T(1) : T(0);
            return l;
        }

        // Upper triangular factor U
        MAT U() const
        {
            const index_t n = this->size();
            MAT u = storage::square(n);
            for (index_t i = 1; i <= n; ++i)
                for (index_t j = 1; j <= n; ++j)
                    u(i,j) = (i <= j) ? this->lu_(i,j) : T(0);
            return u;
        }

        // Permutation matrix P, with P A = L U
        MAT P() const
        {
            MAT p = storage::square(this->size());
            p = T(0);
            for (index_t k = 1; k <= this->size(); ++k) p(k,k) = T(1);
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

        // Solution of A X = B (B: matrix, view or expression; one system per column). The
        // solution is fixed-size if B is.
        template<MatrixAlgebraObject B>
        auto solve(B const& b) const
        {
            detail::solution_t<detail::solver_t<detail::promote_t<T, detail::dtype_t<B>>>, B> x(b);
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


    // LU factorization with dynamic factors
    template<class T, StorageOrder ORDER>
    using LUFactorization = BasicLUFactorization< DynamicMatrix<T, ORDER> >;

    // LU factorization of an N x N fixed-size matrix: fixed-size factors, no heap memory
    template<class T, index_t N, StorageOrder ORDER>
    using FixedLUFactorization = BasicLUFactorization< FixedSizeMatrix<T, N, N, ORDER> >;


    // LU factorization of a copy of a (a square matrix, view or expression), in a floating-point
    // type (double for integer matrices) and the storage order of a. The factors are fixed-size
    // (and the factorization uses no heap memory) if a's dimensions are known at compile time.
    template<MatrixAlgebraObject A>
    auto lu(A const& a)
    {
        typedef detail::solver_t<detail::dtype_t<A>> T;
        constexpr StorageOrder order = detail::solution_order_v<A>;
        if constexpr (detail::traits<A>::fixed_size)
        {
            constexpr index_t N = detail::traits<A>::rows;
            static_assert(detail::traits<A>::cols == N, "lu: the matrix must be square");
            return FixedLUFactorization<T, N, order>(FixedSizeMatrix<T, N, N, order>(a));
        }
        else
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
