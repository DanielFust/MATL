/**************************************************************************************************
 * File: _matl_svd_.hpp
 *
 * Description: Singular value decomposition, A = U diag(sigma) V^T, and the pseudo-inverse.
 *
 *     auto F = svd(A);                decomposes a copy of A (m x n, any shape)
 *     F.U(), F.V()                    m x k and n x k with orthonormal columns (k = min(m, n))
 *     F.singular_values()             sigma, k x 1, in descending order (all >= 0)
 *     F.S()                           diag(sigma), k x k, so that F.U() * F.S() * transpose(F.V()) is A
 *     F.rank(), F.cond()              numerical rank, and the 2-norm condition number sigma_1 / sigma_k
 *     F.pinv()                        the pseudo-inverse A+ (n x m)
 *     X = F.solve(B);                 the minimum-norm least-squares solution, X = A+ B
 *
 *     auto s = singular_values(A);    the singular values only (less work)
 *     auto X = pinv(A);               the pseudo-inverse (also pinv(A, tol))
 *     index_t r = rank(A);            the numerical rank (also rank(A, tol))
 *
 * The singular values of A are the square roots of the eigenvalues of A^T A; the columns of U and
 * V are orthonormal bases of the column space and the row space of A (for the k largest singular
 * values). A singular vector's sign is arbitrary (the signs of U's and V's columns go together).
 *
 * TOLERANCE: singular values at or below tol are treated as zero by rank, pinv and solve. The
 * default, as in MATLAB, is max(m, n) * eps * sigma_1 (eps: the machine epsilon of the type, about
 * 2.2e-16 for double), so that singular values that are only rounding errors do not count. A larger
 * tol truncates small singular values deliberately (regularization of an ill-conditioned problem).
 *
 * Pseudo-inverse: A+ = V diag(1/sigma) U^T over the singular values above tol. For a square
 * nonsingular A it is the inverse; for a full-rank tall A, (A^T A)^-1 A^T; in general,
 * X = A+ B minimizes ||A X - B||, and is the solution of smallest norm among the minimizers (for
 * an underdetermined or rank-deficient system).
 *
 * A fixed-size matrix gives fixed-size results without heap memory. Integer matrices are
 * decomposed in double; float stays float.
 *
 * Algorithm: one-sided Jacobi (Hestenes), applied to A if m >= n and to A^T otherwise. Plane
 * rotations of pairs of columns, accumulated in V, make the columns mutually orthogonal (a pair is
 * rotated while the cosine of its angle exceeds eps); the column norms are then the singular
 * values, and the normalized columns are U. Like the symmetric Jacobi method (_matl_eigen_.hpp), it
 * is accurate (also for small singular values) and simple, and suits small and moderate matrices:
 * a sweep costs about 6 m n^2 operations, and it usually needs 5 to 10 sweeps.
 **************************************************************************************************/
#ifndef MATL_SVD_HPP
#define MATL_SVD_HPP
namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * void one_sided_jacobi(w, v)
        *
        * Rotates pairs of columns of w (rows >= cols) until they are orthogonal to working
        * precision, and applies the same rotations to the columns of v (cols x cols), if not null.
        *******************************************************************************************/
        template<class W, class V>
        void one_sided_jacobi(W& w, V* v)
        {
            typedef dtype_t<W> T;
            constexpr T eps = std::numeric_limits<T>::epsilon();
            const index_t m = w.rows(), n = w.cols();
            for (int sweep = 0; sweep < JACOBI_MAX_SWEEPS; ++sweep)
            {
                bool rotated = false;
                for (index_t p = 1; p < n; ++p)
                    for (index_t q = p + 1; q <= n; ++q)
                    {
                        // The 2 x 2 Gram matrix of columns p and q: [alpha gamma; gamma beta]
                        T alpha = 0, beta = 0, gamma = 0;
                        for (index_t i = 1; i <= m; ++i)
                        {
                            alpha += w(i,p) * w(i,p);
                            beta  += w(i,q) * w(i,q);
                            gamma += w(i,p) * w(i,q);
                        }
                        if (!(std::abs(gamma) > eps * std::sqrt(alpha) * std::sqrt(beta))) continue;
                        const jacobi_rotation<T> r(alpha, beta, gamma);
                        rotate_columns(w, p, q, r);
                        if (v) rotate_columns(*v, p, q, r);
                        rotated = true;
                    }
                if (!rotated) break;
            }
        }

        // Column norms of w (the singular values once one_sided_jacobi has run)
        template<class W, class COL>
        void column_norms(W const& w, COL& sigma)
        {
            typedef dtype_t<W> T;
            for (index_t j = 1; j <= w.cols(); ++j)
            {
                T s = 0;
                for (index_t i = 1; i <= w.rows(); ++i) s += w(i,j) * w(i,j);
                sigma(j,1) = std::sqrt(s);
            }
        }

        // Sorts sigma in descending order, and the columns of w and v (if not null) with it
        template<class COL, class W, class V>
        void sort_descending(COL& sigma, W& w, V* v)
        {
            const index_t k = sigma.rows();
            for (index_t j = 1; j < k; ++j)
            {
                index_t largest = j;
                for (index_t i = j + 1; i <= k; ++i)
                    if (sigma(i,1) > sigma(largest,1)) largest = i;
                if (largest != j)
                {
                    std::swap(sigma(j,1), sigma(largest,1));
                    swap_columns(w, j, largest);
                    if (v) swap_columns(*v, j, largest);
                }
            }
        }

       /*******************************************************************************************
        * void orthonormalize_columns(w, sigma)
        *
        * Divides the columns of w by their norms sigma (descending). Columns with a zero norm
        * (exactly zero singular values, which come last) are replaced by unit vectors orthogonal
        * to all the columns before them (Gram-Schmidt, twice, on e_1, e_2, ...), so that the
        * columns of w are orthonormal.
        *******************************************************************************************/
        template<class W, class COL>
        void orthonormalize_columns(W& w, COL const& sigma)
        {
            typedef dtype_t<W> T;
            const index_t m = w.rows(), k = w.cols();
            for (index_t j = 1; j <= k; ++j)
            {
                if (sigma(j,1) > T(0))
                {
                    const T reciprocal = T(1) / sigma(j,1);
                    for (index_t i = 1; i <= m; ++i) w(i,j) *= reciprocal;
                    continue;
                }
                for (index_t e = 1; e <= m; ++e)
                {
                    for (index_t i = 1; i <= m; ++i) w(i,j) = (i == e) ? T(1) : T(0);
                    for (int pass = 0; pass < 2; ++pass)
                        for (index_t l = 1; l < j; ++l)
                        {
                            T d = 0;
                            for (index_t i = 1; i <= m; ++i) d += w(i,l) * w(i,j);
                            for (index_t i = 1; i <= m; ++i) w(i,j) -= d * w(i,l);
                        }
                    T s = 0;
                    for (index_t i = 1; i <= m; ++i) s += w(i,j) * w(i,j);
                    if (s > T(0.25))                              // e_e was far from the others' span
                    {
                        const T reciprocal = T(1) / std::sqrt(s);
                        for (index_t i = 1; i <= m; ++i) w(i,j) *= reciprocal;
                        break;
                    }
                }
            }
        }
    }


   /***********************************************************************************************
    * class BasicSVD<T, M, N, ORDER>
    *
    * Singular value decomposition A = U diag(sigma) V^T of an m x n matrix (see the description
    * above). Made by svd, and usually named through:
    *
    *     SVD<T, ORDER>              dynamic results (any m x n)
    *     FixedSVD<T, M, N, ORDER>   fixed-size results: no heap memory
    *
    * Template Parameters:
    *     T         Floating-point element type
    *     M, N      Dimensions of A, or 0 and 0 for a dynamic one
    *     ORDER     Storage order of the results
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    class BasicSVD
    {
        static_assert(std::is_floating_point_v<T>, "SVD: the element type must be a floating-point type");
        static_assert((M == 0) == (N == 0), "SVD: both dimensions are fixed, or both are dynamic");

        static constexpr bool    FIXED = M > 0;
        static constexpr index_t K     = std::min(M, N);           // min(m, n), or 0 if dynamic

      public:
        typedef T                                          dtype;
        typedef detail::sized_matrix_t<T, M, N, ORDER>     input_type;      // m x n
        typedef detail::sized_matrix_t<T, M, K, ORDER>     u_type;          // m x k
        typedef detail::sized_matrix_t<T, N, K, ORDER>     v_type;          // n x k
        typedef detail::sized_matrix_t<T, K, 1, ORDER>     vector_type;     // k x 1
        typedef detail::sized_matrix_t<T, K, K, ORDER>     diagonal_type;   // k x k
        typedef detail::sized_matrix_t<T, N, M, ORDER>     pinv_type;       // n x m

      private:
        // The rotations work on columns, so they run on column-major copies
        typedef detail::sized_matrix_t<T, M, K, StorageOrder::ColMajor> u_work;
        typedef detail::sized_matrix_t<T, N, K, StorageOrder::ColMajor> v_work;

        u_type      u_;
        v_type      v_;
        vector_type sigma_;

        // Overwrites x (identity of the given order)
        template<class X>
        static void set_identity(X& x)
        {
            x = T(0);
            for (index_t k = 1; k <= x.rows(); ++k) x(k,k) = T(1);
        }

        // 1 / sigma_j for the singular values above tol, 0 for the others
        T inverse_sigma(index_t j, T tol) const
        {return (this->sigma_(j,1) > tol) ? T(1) / this->sigma_(j,1) : T(0);}

        T resolve(T tol) const {return (tol < T(0)) ? this->tolerance() : tol;}

      public:
        // Decomposes the matrix a
        explicit BasicSVD(input_type const& a)
            : u_(detail::make_sized<u_type>(a.rows(), std::min(a.rows(), a.cols()))),
              v_(detail::make_sized<v_type>(a.cols(), std::min(a.rows(), a.cols()))),
              sigma_(detail::make_sized<vector_type>(std::min(a.rows(), a.cols()), 1))
        {
            const index_t k = std::min(a.rows(), a.cols());
            u_work u = detail::make_sized<u_work>(a.rows(), k);
            v_work v = detail::make_sized<v_work>(a.cols(), k);
            const bool tall = a.rows() >= a.cols();
            if constexpr (!FIXED || M >= N)
                if (tall)
                {
                    // A = W V^T after the rotations; then W = U diag(sigma)
                    u = a;
                    set_identity(v);
                    detail::one_sided_jacobi(u, &v);
                    detail::column_norms(u, this->sigma_);
                    detail::sort_descending(this->sigma_, u, &v);
                    detail::orthonormalize_columns(u, this->sigma_);
                }
            if constexpr (!FIXED || M < N)
                if (!tall)
                {
                    // The same for A^T = V diag(sigma) U^T
                    v = transpose(a);
                    set_identity(u);
                    detail::one_sided_jacobi(v, &u);
                    detail::column_norms(v, this->sigma_);
                    detail::sort_descending(this->sigma_, v, &u);
                    detail::orthonormalize_columns(v, this->sigma_);
                }
            this->u_ = u;
            this->v_ = v;
        }

        // Dimensions of A
        index_t rows() const {return this->u_.rows();}
        index_t cols() const {return this->v_.rows();}

        // Left singular vectors (m x k, orthonormal columns)
        u_type const& U() const {return this->u_;}

        // Right singular vectors (n x k, orthonormal columns)
        v_type const& V() const {return this->v_;}

        // Singular values in descending order (k x 1)
        vector_type const& singular_values() const {return this->sigma_;}

        // diag(sigma), k x k
        diagonal_type S() const
        {
            const index_t k = this->sigma_.rows();
            diagonal_type s = detail::make_sized<diagonal_type>(k, k);
            s = T(0);
            for (index_t j = 1; j <= k; ++j) s(j,j) = this->sigma_(j,1);
            return s;
        }

        // Default tolerance (see TOLERANCE above): max(m, n) * eps * sigma_1
        T tolerance() const
        {
            if (this->sigma_.rows() == 0) return T(0);
            return static_cast<T>(std::max(this->rows(), this->cols())) * std::numeric_limits<T>::epsilon() * this->sigma_(1,1);
        }

        // Number of singular values above tol (default: tolerance())
        index_t rank(T tol = T(-1)) const
        {
            tol = this->resolve(tol);
            index_t r = 0;
            while (r < this->sigma_.rows() && this->sigma_(r + 1, 1) > tol) ++r;
            return r;
        }

        // 2-norm condition number sigma_1 / sigma_k (infinite if sigma_k is 0; 1 for an empty matrix)
        T cond() const
        {
            const index_t k = this->sigma_.rows();
            if (k == 0) return T(1);
            if (this->sigma_(k,1) == T(0)) return std::numeric_limits<T>::infinity();
            return this->sigma_(1,1) / this->sigma_(k,1);
        }

        // Pseudo-inverse A+ = V diag(1/sigma) U^T over the singular values above tol (n x m)
        pinv_type pinv(T tol = T(-1)) const
        {
            tol = this->resolve(tol);
            v_type scaled = this->v_;
            for (index_t j = 1; j <= this->sigma_.rows(); ++j)
            {
                const T f = this->inverse_sigma(j, tol);
                for (index_t i = 1; i <= scaled.rows(); ++i) scaled(i,j) *= f;
            }
            pinv_type x = scaled * transpose(this->u_);
            return x;
        }

        // Minimum-norm least-squares solution X = A+ B (B: m rows, one system per column; tol as
        // for pinv). Fixed-size if A and B are.
        template<MatrixAlgebraObject B>
        auto solve(B const& b, T tol = T(-1)) const
        {
            typedef detail::solver_t<detail::promote_t<T, detail::dtype_t<B>>> R;
            constexpr StorageOrder order = detail::solution_order_v<B>;
            constexpr index_t C = detail::traits<B>::fixed_size ? detail::traits<B>::cols : 0;
            detail::check_right_hand_side("svd solve", this->rows(), b.rows(), b.cols());
            tol = this->resolve(tol);

            detail::sized_matrix_t<R, (C > 0 ? K : 0), C, order> y = transpose(this->u_) * b;     // U^T B
            for (index_t j = 1; j <= this->sigma_.rows(); ++j)
            {
                const R f = this->inverse_sigma(j, tol);
                for (index_t c = 1; c <= y.cols(); ++c) y(j,c) *= f;
            }
            detail::sized_matrix_t<R, (C > 0 ? N : 0), C, order> x = this->v_ * y;               // V diag(1/sigma) U^T B
            return x;
        }
    };

    // Singular value decomposition with dynamic results
    template<class T, StorageOrder ORDER>
    using SVD = BasicSVD<T, 0, 0, ORDER>;

    // Singular value decomposition of an M x N fixed-size matrix: fixed-size results, no heap memory
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    using FixedSVD = BasicSVD<T, M, N, ORDER>;


    // Singular value decomposition of a copy of a (a matrix, view or expression of any shape), in a
    // floating-point type (double for integer matrices) and the storage order of a; fixed-size if
    // a is
    template<MatrixAlgebraObject A>
    auto svd(A const& a)
    {
        typedef detail::solver_t<detail::dtype_t<A>> T;
        constexpr StorageOrder order = detail::solution_order_v<A>;
        if constexpr (detail::traits<A>::fixed_size)
            return FixedSVD<T, detail::traits<A>::rows, detail::traits<A>::cols, order>(
                       FixedSizeMatrix<T, detail::traits<A>::rows, detail::traits<A>::cols, order>(a));
        else
            return SVD<T, order>(DynamicMatrix<T, order>(a));
    }

    // Singular values of a matrix, view or expression, in descending order (min(m, n) x 1;
    // fixed-size if a is). The same values as svd(a).singular_values().
    template<MatrixAlgebraObject A>
    auto singular_values(A const& a)
    {
        typedef detail::solver_t<detail::dtype_t<A>> T;
        constexpr StorageOrder order = detail::solution_order_v<A>;
        constexpr bool    fixed = detail::traits<A>::fixed_size;
        constexpr index_t M = fixed ? detail::traits<A>::rows : 0;
        constexpr index_t N = fixed ? detail::traits<A>::cols : 0;
        typedef detail::sized_matrix_t<T, std::max(M, N), std::min(M, N), StorageOrder::ColMajor> W;   // as svd
        typedef detail::sized_matrix_t<T, std::min(M, N), 1, order> COL;

        const index_t m = a.rows(), n = a.cols();
        W w = detail::make_sized<W>(std::max(m, n), std::min(m, n));
        if constexpr (fixed)                                   // the same work as svd(a)
        {
            if constexpr (M >= N) w = a;
            else                  w = transpose(a);
        }
        else
        {
            if (m >= n) w = a;
            else        w = transpose(a);
        }
        COL sigma = detail::make_sized<COL>(std::min(m, n), 1);
        detail::one_sided_jacobi(w, static_cast<W*>(nullptr));
        detail::column_norms(w, sigma);
        detail::sort_descending(sigma, w, static_cast<W*>(nullptr));
        return sigma;
    }

    // Pseudo-inverse of a matrix, view or expression (see BasicSVD::pinv; tol < 0: the default)
    template<MatrixAlgebraObject A>
    auto pinv(A const& a, detail::solver_t<detail::dtype_t<A>> tol = -1)
    {return svd(a).pinv(tol);}

    // Numerical rank of a matrix, view or expression (see BasicSVD::rank; tol < 0: the default)
    template<MatrixAlgebraObject A>
    index_t rank(A const& a, detail::solver_t<detail::dtype_t<A>> tol = -1)
    {return svd(a).rank(tol);}
}
#endif
