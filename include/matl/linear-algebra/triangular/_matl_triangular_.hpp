/**************************************************************************************************
 * File: _matl_triangular_.hpp
 *
 * Description: Triangular solves (forward and back substitution), and the routines the
 *              factorizations share.
 *
 *     X = solve_lower(L, B)              solves L X = B for a lower triangular L
 *     X = solve_upper(U, B)              solves U X = B for an upper triangular U
 *     solve_lower_inplace(L, X)          overwrites X (a matrix or view holding B) with the solution
 *     solve_upper_inplace(U, X)
 *
 * Only the triangle that is named is read: the other one may hold anything (e.g. the other
 * factor of an LU factorization). With Diagonal::Unit the diagonal is taken to be 1 and is not
 * read either. B may have any number of columns (one system per column).
 *
 * The dimensions are always checked (matl::dimension_error), and a zero on the diagonal is
 * reported as matl::singular_matrix_error.
 **************************************************************************************************/
#ifndef MATL_TRIANGULAR_HPP
#define MATL_TRIANGULAR_HPP
namespace matl
{
    // Diagonal of a triangular matrix: as stored (NonUnit), or taken to be all ones (Unit, e.g. the
    // L factor of an LU factorization, whose diagonal is not stored)
    enum class Diagonal
    {
        NonUnit,
        Unit
    };

    namespace detail
    {
        // Columns per panel of the blocked QR factorization, and the size below which Cholesky's
        // symmetric update stops splitting
        inline constexpr index_t FACTORIZATION_BLOCK = 64;

        // Columns (or rows) at which the recursive factorizations and solves stop splitting and
        // work column by column
        inline constexpr index_t RECURSION_BASE = 16;

        // Floating-point element type in which a system with elements of type T is solved (double
        // for integers)
        template<typename T>
        using solver_t = real_t<T>;

        // Storage order of the result of a solve with right-hand side B (that of B, row-major
        // when B has no single storage order)
        template<typename B>
        inline constexpr StorageOrder solution_order_v =
            (traits<B>::storage_order == StorageOrder::ColMajor) ? StorageOrder::ColMajor : StorageOrder::RowMajor;

        // Storage of a factorization's factors: an n x n dynamic matrix, or a fixed-size N x N one
        // (whose factorization then needs no heap memory: its pivots are a std::array)
        template<class MAT> struct factor_storage;

        template<class T, StorageOrder ORDER>
        struct factor_storage< DynamicMatrix<T, ORDER> >
        {
            typedef T                    dtype;
            typedef std::vector<index_t> pivots;
            static DynamicMatrix<T, ORDER> square(index_t n) {return DynamicMatrix<T, ORDER>(n, n);}
        };

        template<class T, index_t N, StorageOrder ORDER>
        struct factor_storage< FixedSizeMatrix<T, N, N, ORDER> >
        {
            typedef T                       dtype;
            typedef std::array<index_t, N>  pivots;
            static FixedSizeMatrix<T, N, N, ORDER> square(index_t) {return FixedSizeMatrix<T, N, N, ORDER>{};}
        };

        // Type of the solution of a system with right-hand side B, in element type R: fixed-size if
        // B's dimensions are known at compile time, dynamic otherwise; the storage order of B
        template<typename R, typename B>
        using solution_t = std::conditional_t< traits<B>::fixed_size,
                                               FixedSizeMatrix<R, traits<B>::rows, traits<B>::cols, solution_order_v<B>>,
                                               DynamicMatrix<R, solution_order_v<B>> >;

        // A matrix of ROWS x COLS elements: fixed-size if both are known at compile time (> 0),
        // dynamic otherwise (0 stands for a dimension known only at run time)
        template<class T, index_t ROWS, index_t COLS, StorageOrder ORDER>
        using sized_matrix_t = std::conditional_t< (ROWS > 0 && COLS > 0),
                                                   FixedSizeMatrix<T, (ROWS > 0 ? ROWS : 1), (COLS > 0 ? COLS : 1), ORDER>,
                                                   DynamicMatrix<T, ORDER> >;

        // A matrix of type MAT (see sized_matrix_t) with the given dimensions, elements uninitialized
        template<class MAT>
        MAT make_sized(index_t rows, index_t cols)
        {
            if constexpr (traits<MAT>::fixed_size) return MAT{};
            else                                   return MAT(rows, cols);
        }

        // Reports a matrix that is not square
        inline void check_square(char const* operation, index_t rows, index_t cols)
        {
            if (rows != cols)
                raise<dimension_error>(std::string(operation) + " of a non-square " + dimensions(rows, cols) +
                                       " matrix");
        }

        // Reports a right-hand side whose number of rows is not the order of the system
        inline void check_right_hand_side(char const* operation, index_t n, index_t rows, index_t cols)
        {
            if (rows != n)
                raise<dimension_error>(std::string(operation) + ": the right-hand side is " + dimensions(rows, cols) +
                                       ", but the system has " + std::to_string(n) + " equations");
        }

       /*******************************************************************************************
        * void triangular_solve<LOWER, UNIT>(coef, n, x, operation)
        *
        * Solves T X = X in place for the n x n triangular matrix T whose element (i,p) is
        * coef(i,p) (1-based; only the triangle that is solved with is read). X is a matrix or view
        * with n rows; each column is one system.
        *
        * Rows of X are combined as whole rows, X(i,:) -= T(i,p) * X(p,:), so a right-hand side
        * with several columns is processed in one pass.
        *
        * Template Parameters:
        *     LOWER        'true' for forward substitution (lower triangular), 'false' for back
        *                  substitution (upper triangular)
        *     UNIT         'true' if the diagonal is 1 (not read)
        * Parameters:
        *     operation    Name of the operation, for the error message of a zero diagonal
        *******************************************************************************************/
        template<bool LOWER, bool UNIT, typename COEF, typename X>
        void triangular_solve(COEF const& coef, index_t n, X& x, char const* operation)
        {
            const index_t K = x.cols();

            // X(i,:) -= c * X(p,:)
            const auto eliminate = [&x, K](index_t i, index_t p, auto c)
            {
                if (c != 0)
                    for (index_t j = 1; j <= K; ++j)
                        x(i,j) -= c * x(p,j);
            };

            // X(i,:) /= T(i,i)
            const auto divide = [&](index_t i)
            {
                if constexpr (!UNIT)
                {
                    const auto d = coef(i,i);
                    if (d == 0)
                        raise<singular_matrix_error>(std::string(operation) + ": the triangular matrix is singular "
                                                     "(diagonal element " + std::to_string(i) + " is 0)");
                    for (index_t j = 1; j <= K; ++j)
                        x(i,j) /= d;
                }
            };

            if constexpr (LOWER)
            {
                for (index_t i = 1; i <= n; ++i)
                {
                    for (index_t p = 1; p < i; ++p) eliminate(i, p, coef(i,p));
                    divide(i);
                }
            }
            else
            {
                for (index_t i = n; i >= 1; --i)
                {
                    for (index_t p = i + 1; p <= n; ++p) eliminate(i, p, coef(i,p));
                    divide(i);
                }
            }
        }

        // Solves with the lower (LOWER) or upper triangle of the matrix tri, into x
        template<bool LOWER, MatrixAlgebraObject TRI, typename X>
        void solve_triangle(TRI const& tri, X& x, Diagonal diagonal, char const* operation)
        {
            check_square(operation, tri.rows(), tri.cols());
            check_right_hand_side(operation, tri.rows(), x.rows(), x.cols());

            // Matrices and views are read directly; other expressions are evaluated first
            const product_operand_t<TRI> coef(tri);
            if (diagonal == Diagonal::Unit) triangular_solve<LOWER, true>(coef, tri.rows(), x, operation);
            else                            triangular_solve<LOWER, false>(coef, tri.rows(), x, operation);
        }

        // Solution of a triangular system with right-hand side b (a new matrix, fixed-size if b is)
        template<bool LOWER, MatrixAlgebraObject TRI, MatrixAlgebraObject B>
        auto solve_triangle_copy(TRI const& tri, B const& b, Diagonal diagonal, char const* operation)
        {
            typedef solver_t<promote_t<dtype_t<TRI>, dtype_t<B>>> T;
            solution_t<T, B> x(b);
            solve_triangle<LOWER>(tri, x, diagonal, operation);
            return x;
        }

       /*******************************************************************************************
        * T estimate_inverse_norm_1<T>(n, solve, solve_transposed)
        *
        * Estimates the 1-norm of the inverse of an n x n matrix A, given functions that overwrite
        * an n x 1 column x with A^-1 x and with A^-T x. This is Hager's method with Higham's
        * refinements (as in LAPACK's xLACN2): it usually finds the norm exactly, and never
        * overestimates it. It needs a few solves (O(n^2) each), not the inverse (O(n^3)).
        *******************************************************************************************/
        template<typename T, typename SOLVE, typename SOLVE_T>
        T estimate_inverse_norm_1(index_t n, SOLVE const& solve, SOLVE_T const& solve_transposed)
        {
            if (n == 0) return T(0);
            typedef DynamicMatrix<T, StorageOrder::ColMajor> Vector;

            const auto norm_1 = [n](Vector const& v)
            {T s = 0; for (index_t i = 1; i <= n; ++i) s += std::abs(v(i)); return s;};

            // Start from x = (1/n, ..., 1/n)
            Vector x(n, 1);
            x = T(1) / static_cast<T>(n);
            T estimate = 0;
            for (int iteration = 1; iteration <= 5; ++iteration)
            {
                Vector y = x;
                solve(y);
                const T y_norm = norm_1(y);
                if (iteration > 1 && y_norm <= estimate) break;     // no improvement
                estimate = y_norm;

                // z = A^-T sign(y): the largest |z_j| points to the column of A^-1 to try next
                Vector z(n, 1);
                for (index_t i = 1; i <= n; ++i) z(i) = (y(i) >= 0) ? T(1) : T(-1);
                solve_transposed(z);
                index_t j_max = 1;
                T z_x = 0;
                for (index_t i = 1; i <= n; ++i)
                {
                    if (std::abs(z(i)) > std::abs(z(j_max))) j_max = i;
                    z_x += z(i) * x(i);
                }
                if (std::abs(z(j_max)) <= z_x) break;               // converged
                x = T(0);
                x(j_max) = T(1);
            }

            // Higham's alternative vector, which catches matrices the iteration underestimates
            Vector v(n, 1);
            for (index_t i = 1; i <= n; ++i)
            {
                const T magnitude = (n > 1) ? T(1) + static_cast<T>(i - 1) / static_cast<T>(n - 1) : T(1);
                v(i) = (i % 2 == 1) ? magnitude : -magnitude;
            }
            solve(v);
            return std::max(estimate, T(2) * norm_1(v) / (T(3) * static_cast<T>(n)));
        }

        // Reciprocal condition number in the 1-norm from the norm of A and the estimated norm of
        // A^-1 (0 for a matrix that is singular to working precision)
        template<typename T>
        T reciprocal_condition(T norm, T inverse_norm, index_t n)
        {
            if (n == 0) return T(1);
            if (norm == 0 || !(inverse_norm > 0) || !std::isfinite(inverse_norm)) return T(0);
            return (T(1) / inverse_norm) / norm;
        }
    }

    // Solves L X = B for a lower triangular L (only the lower triangle of L is read)
    template<MatrixAlgebraObject TRI, MatrixAlgebraObject B>
    auto solve_lower(TRI const& L, B const& b, Diagonal diagonal = Diagonal::NonUnit)
    {return detail::solve_triangle_copy<true>(L, b, diagonal, "solve_lower");}

    // Solves U X = B for an upper triangular U (only the upper triangle of U is read)
    template<MatrixAlgebraObject TRI, MatrixAlgebraObject B>
    auto solve_upper(TRI const& U, B const& b, Diagonal diagonal = Diagonal::NonUnit)
    {return detail::solve_triangle_copy<false>(U, b, diagonal, "solve_upper");}

    // Overwrites x (a matrix or view holding the right-hand side) with the solution of L X = B
    template<MatrixAlgebraObject TRI, typename X> requires WritableMatrix<X>
    void solve_lower_inplace(TRI const& L, X&& x, Diagonal diagonal = Diagonal::NonUnit)
    {detail::solve_triangle<true>(L, x, diagonal, "solve_lower_inplace");}

    // Overwrites x (a matrix or view holding the right-hand side) with the solution of U X = B
    template<MatrixAlgebraObject TRI, typename X> requires WritableMatrix<X>
    void solve_upper_inplace(TRI const& U, X&& x, Diagonal diagonal = Diagonal::NonUnit)
    {detail::solve_triangle<false>(U, x, diagonal, "solve_upper_inplace");}
}
#endif
