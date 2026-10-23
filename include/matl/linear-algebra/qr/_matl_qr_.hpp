/**************************************************************************************************
 * File: _matl_qr_.hpp
 *
 * Description: QR factorization, A = Q R, by Householder reflections, and least squares.
 *
 *     auto F = qr(A);                 factorizes a copy of A (m x n, any shape)
 *     auto F = qr_inplace(A);         factorizes A in its own memory (A: a matrix of a
 *                                     floating-point type), which then holds the factors
 *     X = F.solve(B);                 least-squares solution: X minimizes ||A X - B|| (m >= n;
 *                                     for a square A, the solution of A X = B)
 *     X = qr_solve(A, B);             factorizes and solves once
 *
 *     F.Q(), F.R()                    the thin factors: Q is m x k with orthonormal columns, R is
 *                                     k x n upper triangular (k = min(m, n)), and Q R = A
 *     F.full_Q()                      the m x m orthogonal Q
 *     F.Qt_times(B), F.Q_times(B)     Q^T B and Q B, without forming Q
 *     F.rcond()                       estimated reciprocal condition number of R (1-norm)
 *     F.factors(), F.tau()            the packed factors and the reflector coefficients
 *
 * Q is the product of k Householder reflections, Q = H_1 H_2 ... H_k, with H_j = I - tau_j v_j v_j^T.
 * As in LAPACK, the factors are packed in one matrix: R on and above the diagonal, and each v_j
 * below the diagonal of column j (its first element, 1, is not stored).
 *
 * Least squares: with A = Q R and m >= n, the X that minimizes ||A X - B|| solves
 * R X = (Q^T B)(1:n, :). Unlike the normal equations (A^T A X = A^T B), this does not square the
 * condition number of A, so it is the accurate way to fit a model to data.
 *
 * Rank: QR exists for every matrix, so the factorization itself never fails. A solve needs R to be
 * nonsingular (A of full column rank): an exactly zero diagonal element of R is reported as
 * matl::singular_matrix_error, and rcond() shows how close R is to singular. (Column pivoting,
 * which reveals the rank, is not implemented.)
 *
 * Algorithm: blocked like LAPACK's xGEQRF. Each panel of BLOCK columns is factorized column by
 * column; its reflectors are then combined into one block reflector I - V T V^T (the "compact WY"
 * form), which is applied to the rest of the matrix with three matrix products on the product
 * kernel.
 **************************************************************************************************/
#ifndef MATL_QR_HPP
#define MATL_QR_HPP
namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * T householder(DynamicMatrix<T, ORDER>& a, index_t j)
        *
        * Computes the Householder reflection H = I - tau v v^T that maps a(j:m, j) to
        * (beta, 0, ..., 0), as LAPACK's xLARFG: a(j,j) becomes beta and a(j+1:m, j) the elements
        * of v after its first (which is 1). Returns tau (0 when a(j+1:m, j) is already zero, H = I).
        *******************************************************************************************/
        template<class T, StorageOrder ORDER>
        T householder(DynamicMatrix<T, ORDER>& a, index_t j)
        {
            const index_t m = a.rows();

            // Norm of a(j+1:m, j), scaled to avoid overflow and underflow
            T scale = 0;
            for (index_t i = j + 1; i <= m; ++i) scale = std::max(scale, std::abs(a(i,j)));
            if (scale == T(0)) return T(0);
            T sum = 0;
            for (index_t i = j + 1; i <= m; ++i) {const T t = a(i,j) / scale; sum += t * t;}
            const T norm = scale * std::sqrt(sum);

            // beta has the opposite sign of alpha, so that alpha - beta does not cancel
            const T alpha = a(j,j);
            const T beta  = -std::copysign(std::hypot(alpha, norm), alpha);
            const T tau   = (beta - alpha) / beta;
            const T reciprocal = T(1) / (alpha - beta);
            for (index_t i = j + 1; i <= m; ++i) a(i,j) *= reciprocal;
            a(j,j) = beta;
            return tau;
        }

       /*******************************************************************************************
        * void apply_householder(v, j, tau, x, c0, c1)
        *
        * Applies H = I - tau v v^T from the left to rows j..m of columns c0..c1 of x, where v is
        * stored below the diagonal of column j of 'v' (first element 1, not stored):
        * x(:,c) -= tau v (v^T x(:,c)). Visits x along its rows or down its columns, to match its
        * storage order.
        *******************************************************************************************/
        template<typename V, typename X>
        void apply_householder(V const& v, index_t j, dtype_t<V> tau, X& x, index_t c0, index_t c1)
        {
            typedef dtype_t<X> T;
            if (tau == 0 || c0 > c1) return;
            const index_t m = x.rows();

            if constexpr (traits<X>::storage_order == StorageOrder::RowMajor)
            {
                // w(c) = tau * v^T x(:,c), accumulated row by row
                std::vector<T> w(c1 - c0 + 1);
                for (index_t c = c0; c <= c1; ++c) w[c - c0] = x(j,c);
                for (index_t i = j + 1; i <= m; ++i)
                {
                    const T vi = v(i,j);
                    for (index_t c = c0; c <= c1; ++c) w[c - c0] += vi * x(i,c);
                }
                for (T& wc : w) wc *= tau;
                for (index_t c = c0; c <= c1; ++c) x(j,c) -= w[c - c0];
                for (index_t i = j + 1; i <= m; ++i)
                {
                    const T vi = v(i,j);
                    for (index_t c = c0; c <= c1; ++c) x(i,c) -= vi * w[c - c0];
                }
            }
            else
            {
                for (index_t c = c0; c <= c1; ++c)
                {
                    T w = x(j,c);
                    for (index_t i = j + 1; i <= m; ++i) w += v(i,j) * x(i,c);
                    w *= tau;
                    x(j,c) -= w;
                    for (index_t i = j + 1; i <= m; ++i) x(i,c) -= v(i,j) * w;
                }
            }
        }

       /*******************************************************************************************
        * void qr_factorize(DynamicMatrix<T, ORDER>& a, std::vector<T>& tau)
        *
        * Overwrites a with its packed QR factors (see the description above) and fills in the
        * reflector coefficients tau (min(m, n) of them).
        *******************************************************************************************/
        template<class T, StorageOrder ORDER>
        void qr_factorize(DynamicMatrix<T, ORDER>& a, std::vector<T>& tau)
        {
            const index_t m = a.rows();
            const index_t n = a.cols();
            const index_t k = std::min(m, n);
            tau.assign(k, T(0));

            for (index_t k0 = 1; k0 <= k; k0 += FACTORIZATION_BLOCK)
            {
                const index_t k1 = std::min(k0 + FACTORIZATION_BLOCK - 1, k);   // last column of the panel
                const index_t kb = k1 - k0 + 1;

                // Panel: one reflector per column, each applied to the panel's columns right of it
                for (index_t j = k0; j <= k1; ++j)
                {
                    tau[j - 1] = householder(a, j);
                    apply_householder(a, j, tau[j - 1], a, j + 1, k1);
                }

                if (k1 == n) break;

                // The panel's reflectors H_k0 ... H_k1 as one block reflector I - V T V^T:
                // V (rows k0..m) holds the vectors (unit diagonal, zeros above it), and T is
                // upper triangular (LAPACK's xLARFT)
                const index_t mv = m - k0 + 1;
                DynamicMatrix<T, ORDER> V(mv, kb), Tb(kb, kb);
                for (index_t i = 1; i <= mv; ++i)
                    for (index_t c = 1; c <= kb; ++c)
                        V(i,c) = (i < c) ? T(0) : (i == c) ? T(1) : a(k0 + i - 1, k0 + c - 1);
                Tb = T(0);
                for (index_t c = 1; c <= kb; ++c)
                {
                    const T t = tau[k0 + c - 2];
                    Tb(c,c) = t;
                    // T(1:c-1, c) = -t T(1:c-1, 1:c-1) V(:, 1:c-1)^T v_c
                    for (index_t p = 1; p < c; ++p)
                    {
                        T s = 0;
                        for (index_t i = c; i <= mv; ++i) s += V(i,p) * V(i,c);
                        Tb(p,c) = -t * s;
                    }
                    for (index_t p = 1; p < c; ++p)          // multiply by the upper triangular T(1:c-1, 1:c-1)
                    {
                        T s = 0;
                        for (index_t q = p; q < c; ++q) s += Tb(p,q) * Tb(q,c);
                        Tb(p,c) = s;
                    }
                }

                // Apply Q_panel^T = I - V T^T V^T to the trailing columns: A2 -= V (T^T (V^T A2))
                auto A2 = a(range(k0, m), range(k1 + 1, n));
                DynamicMatrix<T, ORDER> W(kb, n - k1);
                W = transpose(V) * A2;
                W = transpose(Tb) * W;
                A2 -= V * W;
            }
        }
    }


   /***********************************************************************************************
    * class QRFactorization<T, ORDER>
    *
    * QR factorization A = Q R of an m x n matrix by Householder reflections (see the description
    * above). Made by qr and qr_inplace.
    *
    * Template Parameters:
    *     T         Floating-point element type
    *     ORDER     Storage order of the factors
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    class QRFactorization
    {
        static_assert(std::is_floating_point_v<T>, "QRFactorization: the element type must be a floating-point type");

      private:
        DynamicMatrix<T, ORDER> qr_;     // R on and above the diagonal, the reflectors below it
        std::vector<T>          tau_;    // reflector coefficients

        // Applies the reflectors to x (m rows): Q^T x (H_1 first), or Q x (H_k first)
        template<typename X>
        void apply_q(X& x, bool transposed) const
        {
            const index_t k = this->tau_.size();
            if (transposed) for (index_t j = 1; j <= k; ++j) detail::apply_householder(this->qr_, j, this->tau_[j - 1], x, 1, x.cols());
            else            for (index_t j = k; j >= 1; --j) detail::apply_householder(this->qr_, j, this->tau_[j - 1], x, 1, x.cols());
        }

        template<typename B>
        using result_t = DynamicMatrix<detail::solver_t<detail::promote_t<T, detail::dtype_t<B>>>, detail::solution_order_v<B>>;

      public:
        // Factorizes the matrix a (its memory may be a matrix's, see qr_inplace)
        explicit QRFactorization(DynamicMatrix<T, ORDER> a) : qr_(std::move(a)), tau_()
        {detail::qr_factorize(this->qr_, this->tau_);}

        // Dimensions of A
        index_t rows() const {return this->qr_.rows();}
        index_t cols() const {return this->qr_.cols();}

        // Packed factors: R on and above the diagonal, the Householder vectors below it
        DynamicMatrix<T, ORDER> const& factors() const {return this->qr_;}

        // Reflector coefficients: H_j = I - tau()[j-1] v_j v_j^T
        std::vector<T> const& tau() const {return this->tau_;}

        // Upper triangular (trapezoidal) factor R, min(m, n) x n
        DynamicMatrix<T, ORDER> R() const
        {
            const index_t k = std::min(this->rows(), this->cols());
            DynamicMatrix<T, ORDER> r(k, this->cols());
            for (index_t i = 1; i <= k; ++i)
                for (index_t j = 1; j <= this->cols(); ++j)
                    r(i,j) = (i <= j) ? this->qr_(i,j) : T(0);
            return r;
        }

        // Thin Q: m x min(m, n) with orthonormal columns, Q() * R() = A
        DynamicMatrix<T, ORDER> Q() const
        {
            const index_t k = std::min(this->rows(), this->cols());
            DynamicMatrix<T, ORDER> q(this->rows(), k);
            q = T(0);
            for (index_t i = 1; i <= k; ++i) q(i,i) = T(1);
            this->apply_q(q, false);
            return q;
        }

        // Full Q: the m x m orthogonal matrix
        DynamicMatrix<T, ORDER> full_Q() const
        {
            DynamicMatrix<T, ORDER> q = eye<T>(this->rows());
            this->apply_q(q, false);
            return q;
        }

        // Q^T B (B with m rows), without forming Q
        template<MatrixAlgebraObject B>
        auto Qt_times(B const& b) const
        {
            detail::check_right_hand_side("qr Qt_times", this->rows(), b.rows(), b.cols());
            result_t<B> x(b);
            this->apply_q(x, true);
            return x;
        }

        // Q B (B with m rows), without forming Q
        template<MatrixAlgebraObject B>
        auto Q_times(B const& b) const
        {
            detail::check_right_hand_side("qr Q_times", this->rows(), b.rows(), b.cols());
            result_t<B> x(b);
            this->apply_q(x, false);
            return x;
        }

        // Least-squares solution of A X = B: the n x (columns of B) matrix X that minimizes
        // ||A X - B|| in each column (A must have at least as many rows as columns, and full
        // column rank)
        template<MatrixAlgebraObject B>
        auto solve(B const& b) const
        {
            const index_t m = this->rows();
            const index_t n = this->cols();
            if (m < n)
                detail::raise<dimension_error>("qr solve: the system is underdetermined (" + detail::dimensions(m, n) +
                                       "); least squares needs at least as many rows as columns");
            for (index_t j = 1; j <= n; ++j)
                if (this->qr_(j,j) == T(0))
                    detail::raise<singular_matrix_error>("qr solve: the matrix is rank deficient (R(" + std::to_string(j) + "," +
                                                 std::to_string(j) + ") is 0)");

            result_t<B> y = this->Qt_times(b);                          // Q^T B
            auto top = y(range(1, n), all);                            // R X = (Q^T B)(1:n, :)
            detail::triangular_solve<false, false>(this->qr_, n, top, "qr solve");
            return result_t<B>(top);
        }

        // Estimated reciprocal condition number (1-norm) of the leading min(m, n) x min(m, n)
        // block of R; for a matrix of full column rank, near 0 means nearly rank deficient
        T rcond() const
        {
            const index_t k = std::min(this->rows(), this->cols());
            T norm = 0;                                                   // 1-norm of the triangle
            for (index_t j = 1; j <= k; ++j)
            {
                T s = 0;
                for (index_t i = 1; i <= j; ++i) s += std::abs(this->qr_(i,j));
                norm = std::max(norm, s);
            }
            const auto solve = [this, k](auto& x) {detail::triangular_solve<false, false>(this->qr_, k, x, "qr rcond");};
            const auto solve_transposed = [this, k](auto& x)
            {
                const auto rt = [this](index_t i, index_t p) {return this->qr_(p,i);};
                detail::triangular_solve<true, false>(rt, k, x, "qr rcond");
            };
            return detail::reciprocal_condition(norm, detail::estimate_inverse_norm_1<T>(k, solve, solve_transposed), k);
        }
    };


    // QR factorization of a copy of a (a matrix, view or expression of any shape), in a
    // floating-point type (double for integer matrices) and the storage order of a
    template<MatrixAlgebraObject A>
    auto qr(A const& a)
    {
        typedef detail::solver_t<detail::dtype_t<A>> T;
        constexpr StorageOrder order = detail::solution_order_v<A>;
        return QRFactorization<T, order>(DynamicMatrix<T, order>(a));
    }

    // QR factorization in the memory of a (a matrix of a floating-point type): a then holds the
    // packed factors, and must outlive the factorization
    template<class T, StorageOrder ORDER>
    QRFactorization<T, ORDER> qr_inplace(DynamicMatrix<T, ORDER>& a)
    {return QRFactorization<T, ORDER>(DynamicMatrix<T, ORDER>(a.rows(), a.cols(), a.data()));}

    template<class T, index_t M, index_t N, StorageOrder ORDER>
    QRFactorization<T, ORDER> qr_inplace(FixedSizeMatrix<T, M, N, ORDER>& a)
    {return QRFactorization<T, ORDER>(DynamicMatrix<T, ORDER>(M, N, a.data()));}

    // Least-squares solution of A X = B by QR factorization (A m x n with m >= n; for a square A,
    // the solution of A X = B)
    template<MatrixAlgebraObject A, MatrixAlgebraObject B>
    auto qr_solve(A const& a, B const& b)
    {return qr(a).solve(b);}
}
#endif
