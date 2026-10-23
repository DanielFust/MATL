/**************************************************************************************************
 * File: _matl_iterative_solvers_.hpp
 *
 * Description: Iterative solvers for A x = b, mainly for large sparse systems:
 *
 *     cg(A, b, x)          conjugate gradients, for symmetric positive definite A
 *     bicgstab(A, b, x)    BiCGSTAB, for general (nonsymmetric) A
 *     gmres(A, b, x)       restarted GMRES(m), for general A (the most robust; more memory)
 *
 *     auto result = cg(A, b, x, {.tolerance = 1e-8}, ILU0Preconditioner(A));
 *     if (!result.converged) ...
 *
 * x holds the initial guess (e.g. zeros, or the previous time step's solution) and receives the
 * approximate solution. A is a sparse matrix, a dense matrix, or a linear operator given as a
 * function op(in, out) that sets out = A in (both n x 1 matrices), so that a matrix need not be
 * stored. An optional preconditioner (see _matl_preconditioners_.hpp) speeds up convergence.
 *
 * Each solver stops when the relative residual ||b - A x|| / ||b|| (2-norm) is at most
 * options.tolerance, or after options.max_iterations iterations (0: the order of the system), and
 * returns an IterativeResult: whether it converged, the number of iterations, and the final
 * relative residual (recomputed from x, not the solver's running estimate). Not converging is not
 * an error: check result.converged. Wrong dimensions are reported as matl::dimension_error, and a
 * matrix that CG finds is not positive definite (p^T A p <= 0) as
 * matl::not_positive_definite_error.
 *
 * The vectors are n x 1 column-major matrices of the floating-point type of b (double for
 * integers); x may be any n x 1 matrix or view of a floating-point type.
 **************************************************************************************************/
#ifndef MATL_ITERATIVE_SOLVERS_HPP
#define MATL_ITERATIVE_SOLVERS_HPP
namespace matl
{
    // Stopping criteria of the iterative solvers
    struct IterativeOptions
    {
        double  tolerance      = 1e-10;   // stop at ||b - A x|| / ||b|| <= tolerance
        index_t max_iterations = 0;       // at most this many iterations (0: the order of the system)
        index_t restart        = 30;      // GMRES: Krylov vectors per cycle (memory: restart + 1 vectors)
    };

    // Outcome of an iterative solve
    template<typename T>
    struct IterativeResult
    {
        bool    converged  = false;       // the tolerance was reached
        index_t iterations = 0;           // iterations taken (for GMRES, over all cycles)
        T       residual   = 0;           // final ||b - A x|| / ||b||
    };

    namespace detail
    {
        template<typename T>
        using krylov_vector = DynamicMatrix<T, StorageOrder::ColMajor>;

        // out = A in, for a sparse or dense matrix, or a linear operator function op(in, out)
        template<typename A, typename V, typename W>
        void apply_operator(A const& a, V const& in, W& out)
        {
            if constexpr (SparseMatrixType<A> || MatrixAlgebraObject<A>) out = a * in;
            else                                                           a(in, out);
        }

        // Order of the system, checked against b and x
        template<typename A, typename B, typename X>
        index_t check_system(char const* solver, A const& a, B const& b, X const& x)
        {
            const index_t n = b.rows();
            if constexpr (SparseMatrixType<A> || MatrixAlgebraObject<A>)
            {
                check_square(solver, a.rows(), a.cols());
                if (a.rows() != n)
                    raise<dimension_error>(std::string(solver) + ": the matrix is " + dimensions(a.rows(), a.cols()) +
                                           ", but b is " + dimensions(b.rows(), b.cols()));
            }
            if (b.cols() != 1 || x.rows() != n || x.cols() != 1)
                raise<dimension_error>(std::string(solver) + ": b and x must be n x 1 vectors of the system's order (b is " +
                                       dimensions(b.rows(), b.cols()) + ", x " + dimensions(x.rows(), x.cols()) + ")");
            return n;
        }

        // 2-norm of a vector
        template<typename V>
        auto norm_2(V const& v) {return norm_fro(v);}

        // The result with the residual recomputed from x
        template<typename T, typename A, typename V, typename X>
        IterativeResult<T> finish(A const& a, V const& b, X const& x, T b_norm, index_t iterations, T tolerance)
        {
            krylov_vector<T> ax(b.rows(), 1);
            krylov_vector<T> xv(x);
            apply_operator(a, xv, ax);
            const T residual = norm_2(b - ax) / b_norm;
            return IterativeResult<T>{residual <= tolerance, iterations, residual};
        }
    }

   /***********************************************************************************************
    * IterativeResult<T> cg(A, b, x, options, preconditioner)
    *
    * Preconditioned conjugate gradients for a symmetric positive definite A (and a symmetric
    * positive definite preconditioner, e.g. Jacobi, or ILU(0) of a symmetric matrix). Each
    * iteration costs one product with A, one preconditioner application, and a few vector
    * operations; in exact arithmetic CG converges in at most n iterations.
    ***********************************************************************************************/
    template<typename A, MatrixAlgebraObject B, typename X, typename P = IdentityPreconditioner>
        requires WritableMatrix<X>
    auto cg(A const& a, B const& b_in, X&& x, IterativeOptions const& options = {}, P const& preconditioner = {})
    {
        typedef detail::solver_t<detail::dtype_t<B>> T;
        const index_t n = detail::check_system("cg", a, b_in, x);
        const index_t max_iterations = options.max_iterations ? options.max_iterations : n;
        const T tolerance = static_cast<T>(options.tolerance);

        const detail::krylov_vector<T> b(b_in);
        const T b_norm = detail::norm_2(b);
        if (b_norm == T(0)) {x = T(0); return IterativeResult<T>{true, 0, T(0)};}

        detail::krylov_vector<T> r(n, 1), z(n, 1), p(n, 1), q(n, 1), xv(x);
        detail::apply_operator(a, xv, q);
        r = b - q;
        preconditioner.apply(r, z);
        p = z;
        T rz = dot(r, z);

        index_t k = 0;
        while (k < max_iterations && detail::norm_2(r) / b_norm > tolerance)
        {
            ++k;
            detail::apply_operator(a, p, q);
            const T pq = dot(p, q);
            if (!(pq > T(0)))
                detail::raise<not_positive_definite_error>("cg: the matrix (or the preconditioner) is not positive definite "
                                                           "(p^T A p <= 0 at iteration " + std::to_string(k) + ")");
            const T alpha = rz / pq;
            xv += alpha * p;
            r -= alpha * q;
            preconditioner.apply(r, z);
            const T rz_next = dot(r, z);
            p = z + (rz_next / rz) * p;
            rz = rz_next;
        }
        x = xv;
        return detail::finish(a, b, x, b_norm, k, tolerance);
    }

   /***********************************************************************************************
    * IterativeResult<T> bicgstab(A, b, x, options, preconditioner)
    *
    * Preconditioned BiCGSTAB (van der Vorst) for a general A, with right preconditioning. Each
    * iteration costs two products with A and two preconditioner applications. BiCGSTAB can break
    * down (a zero denominator); it then stops and reports that it did not converge.
    ***********************************************************************************************/
    template<typename A, MatrixAlgebraObject B, typename X, typename P = IdentityPreconditioner>
        requires WritableMatrix<X>
    auto bicgstab(A const& a, B const& b_in, X&& x, IterativeOptions const& options = {}, P const& preconditioner = {})
    {
        typedef detail::solver_t<detail::dtype_t<B>> T;
        const index_t n = detail::check_system("bicgstab", a, b_in, x);
        const index_t max_iterations = options.max_iterations ? options.max_iterations : n;
        const T tolerance = static_cast<T>(options.tolerance);

        const detail::krylov_vector<T> b(b_in);
        const T b_norm = detail::norm_2(b);
        if (b_norm == T(0)) {x = T(0); return IterativeResult<T>{true, 0, T(0)};}

        detail::krylov_vector<T> r(n, 1), r0(n, 1), p(n, 1), v(n, 1), s(n, 1), t(n, 1), ph(n, 1), sh(n, 1), xv(x);
        detail::apply_operator(a, xv, v);
        r = b - v;
        r0 = r;                                                    // the shadow residual
        p = T(0);
        v = T(0);
        T rho = 1, alpha = 1, omega = 1;

        index_t k = 0;
        while (k < max_iterations && detail::norm_2(r) / b_norm > tolerance)
        {
            ++k;
            const T rho_next = dot(r0, r);
            if (rho_next == T(0)) break;                           // breakdown
            p = r + ((rho_next / rho) * (alpha / omega)) * (p - omega * v);
            preconditioner.apply(p, ph);
            detail::apply_operator(a, ph, v);
            const T r0v = dot(r0, v);
            if (r0v == T(0)) break;                                // breakdown
            alpha = rho_next / r0v;
            s = r - alpha * v;
            if (detail::norm_2(s) / b_norm <= tolerance) {xv += alpha * ph; r = s; break;}
            preconditioner.apply(s, sh);
            detail::apply_operator(a, sh, t);
            const T tt = dot(t, t);
            if (tt == T(0)) {xv += alpha * ph; r = s; break;}
            omega = dot(t, s) / tt;
            xv += alpha * ph + omega * sh;
            r = s - omega * t;
            if (omega == T(0)) break;                              // breakdown
            rho = rho_next;
        }
        x = xv;
        return detail::finish(a, b, x, b_norm, k, tolerance);
    }

   /***********************************************************************************************
    * IterativeResult<T> gmres(A, b, x, options, preconditioner)
    *
    * Restarted GMRES(m) (Saad and Schultz) for a general A, with right preconditioning, modified
    * Gram-Schmidt and Givens rotations: each cycle builds an orthonormal basis of up to
    * m = options.restart Krylov vectors and takes the x that minimizes the residual over it. Each
    * iteration costs one product with A, one preconditioner application, and orthogonalization
    * against the cycle's vectors (O(n m)).
    ***********************************************************************************************/
    template<typename A, MatrixAlgebraObject B, typename X, typename P = IdentityPreconditioner>
        requires WritableMatrix<X>
    auto gmres(A const& a, B const& b_in, X&& x, IterativeOptions const& options = {}, P const& preconditioner = {})
    {
        typedef detail::solver_t<detail::dtype_t<B>> T;
        const index_t n = detail::check_system("gmres", a, b_in, x);
        const index_t max_iterations = options.max_iterations ? options.max_iterations : n;
        const index_t m = std::max<index_t>(1, std::min(options.restart, n));
        const T tolerance = static_cast<T>(options.tolerance);

        const detail::krylov_vector<T> b(b_in);
        const T b_norm = detail::norm_2(b);
        if (b_norm == T(0)) {x = T(0); return IterativeResult<T>{true, 0, T(0)};}

        detail::krylov_vector<T> xv(x), r(n, 1), w(n, 1), z(n, 1), y(m, 1);
        DynamicMatrix<T, StorageOrder::ColMajor> V(n, m + 1), H(m + 1, m);
        std::vector<T> cs(m), sn(m), g(m + 1);

        index_t k = 0;
        bool converged = false;
        while (k < max_iterations && !converged)
        {
            // r = b - A x, and the first basis vector
            detail::apply_operator(a, xv, w);
            r = b - w;
            const T beta = detail::norm_2(r);
            if (beta / b_norm <= tolerance) {converged = true; break;}
            V(all, 1) = r / beta;
            std::fill(g.begin(), g.end(), T(0));
            g[0] = beta;
            H = T(0);

            index_t j = 0;
            while (j < m && k < max_iterations)
            {
                ++j;
                ++k;
                // w = A M^-1 v_j, orthogonalized against v_1..v_j (modified Gram-Schmidt)
                preconditioner.apply(V(all, j), z);
                detail::apply_operator(a, z, w);
                for (index_t i = 1; i <= j; ++i)
                {
                    H(i, j) = dot(w, V(all, i));
                    w -= H(i, j) * V(all, i);
                }
                H(j + 1, j) = detail::norm_2(w);
                if (H(j + 1, j) != T(0)) V(all, j + 1) = w / H(j + 1, j);

                // Rotate column j of H with the previous rotations, then eliminate H(j+1,j)
                for (index_t i = 1; i < j; ++i)
                {
                    const T h1 = H(i, j), h2 = H(i + 1, j);
                    H(i, j)     =  cs[i - 1] * h1 + sn[i - 1] * h2;
                    H(i + 1, j) = -sn[i - 1] * h1 + cs[i - 1] * h2;
                }
                const T d = std::hypot(H(j, j), H(j + 1, j));
                cs[j - 1] = (d == T(0)) ? T(1) : H(j, j) / d;
                sn[j - 1] = (d == T(0)) ? T(0) : H(j + 1, j) / d;
                H(j, j) = d;
                H(j + 1, j) = T(0);
                g[j]     = -sn[j - 1] * g[j - 1];
                g[j - 1] =  cs[j - 1] * g[j - 1];

                if (std::abs(g[j]) / b_norm <= tolerance) {converged = true; break;}
                if (d == T(0)) break;                               // the Krylov space is exhausted
            }

            // x += M^-1 (V(:, 1:j) y), with R y = g(1:j) (R: the rotated, upper triangular H)
            if (j > 0)
            {
                auto yj = y(range(1, j), all);
                for (index_t i = 1; i <= j; ++i) yj(i, 1) = g[i - 1];
                detail::triangular_solve<false, false>(H, j, yj, "gmres");
                w = V(all, range(1, j)) * yj;
                preconditioner.apply(w, z);
                xv += z;
            }
        }
        x = xv;
        return detail::finish(a, b, x, b_norm, k, tolerance);
    }
}
#endif
