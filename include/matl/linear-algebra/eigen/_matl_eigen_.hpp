/**************************************************************************************************
 * File: _matl_eigen_.hpp
 *
 * Description: Eigenvalues and eigenvectors of symmetric matrices, A = V diag(lambda) V^T.
 *
 *     auto E = eig_sym(A);            eigenvalues and eigenvectors of a symmetric matrix (only its
 *                                     lower triangle is read)
 *     E.values()                      the eigenvalues, n x 1, in ascending order
 *     E.vectors()                     n x n with orthonormal columns: column k is an eigenvector
 *                                     for values()(k), and A V = V diag(values)
 *     auto lambda = eigenvalues_sym(A);   the eigenvalues only (about half the work)
 *
 * The eigenvalues of a symmetric matrix are real, and its eigenvectors can be chosen orthonormal.
 * An eigenvector's sign is arbitrary, and for a repeated eigenvalue any orthonormal basis of its
 * eigenspace may be returned.
 *
 * A fixed-size matrix gives fixed-size results without heap memory (e.g. the principal stresses
 * and directions of a 3 x 3 stress tensor). Integer matrices are decomposed in double; float stays
 * float.
 *
 * Algorithm: cyclic Jacobi. Each step is a plane rotation, V^T A V, that zeroes one off-diagonal
 * pair (p, q); a sweep visits every pair, and the off-diagonal part shrinks quadratically once it
 * is small (usually 5 to 10 sweeps). A pair is rotated while |a(p,q)| > eps sqrt(|a(p,p) a(q,q)|),
 * the criterion of Demmel and Veselic, which keeps small eigenvalues accurate relative to their own
 * size. Each sweep costs about 4 n^3 operations (2 n^3 without the eigenvectors), so the method
 * suits small and moderate matrices (up to a few hundred): it is simple, robust and accurate, not
 * as fast as LAPACK's tridiagonal methods for large ones.
 **************************************************************************************************/
#ifndef MATL_EIGEN_HPP
#define MATL_EIGEN_HPP
namespace matl
{
    namespace detail
    {
        // Upper limit on the sweeps of the Jacobi methods (they need about 5 to 10)
        inline constexpr int JACOBI_MAX_SWEEPS = 60;

        // Plane rotation [c s; -s c] that diagonalizes the symmetric 2 x 2 matrix
        // [alpha gamma; gamma beta] (gamma != 0), and t = s / c, the smaller of the two possible
        // angles' tangents (|t| <= 1)
        template<typename T>
        struct jacobi_rotation
        {
            T c, s, t;

            jacobi_rotation(T alpha, T beta, T gamma)
            {
                const T zeta = (beta - alpha) / (T(2) * gamma);
                t = std::copysign(T(1), zeta) / (std::abs(zeta) + std::hypot(T(1), zeta));
                c = T(1) / std::sqrt(T(1) + t * t);
                s = c * t;
            }
        };

        // Rotates columns p and q of x: x(:,p) = c x(:,p) - s x(:,q), x(:,q) = s x(:,p) + c x(:,q)
        template<class X, typename T>
        void rotate_columns(X& x, index_t p, index_t q, jacobi_rotation<T> const& r)
        {
            for (index_t i = 1; i <= x.rows(); ++i)
            {
                const T g = x(i,p), h = x(i,q);
                x(i,p) = r.c * g - r.s * h;
                x(i,q) = r.s * g + r.c * h;
            }
        }

        // Swaps columns p and q of x
        template<class X>
        void swap_columns(X& x, index_t p, index_t q)
        {
            for (index_t i = 1; i <= x.rows(); ++i) std::swap(x(i,p), x(i,q));
        }

       /*******************************************************************************************
        * void symmetric_eigen(a, values, vectors)
        *
        * Overwrites values (n x 1) with the eigenvalues of the symmetric matrix a (n x n, only its
        * lower triangle is read; a is overwritten), in ascending order, and, if vectors is not
        * null, vectors (n x n) with orthonormal eigenvectors in the same order (cyclic Jacobi, see
        * the description above). Fastest for column-major a and vectors, whose columns the
        * rotations update.
        *******************************************************************************************/
        template<class MAT, class COL>
        void symmetric_eigen(MAT& a, COL& values, MAT* vectors)
        {
            typedef dtype_t<MAT> T;
            constexpr T eps = std::numeric_limits<T>::epsilon();
            const index_t n = a.rows();

            for (index_t i = 1; i <= n; ++i)                  // the upper triangle from the lower
                for (index_t j = i + 1; j <= n; ++j) a(i,j) = a(j,i);
            if (vectors)
            {
                *vectors = T(0);
                for (index_t k = 1; k <= n; ++k) (*vectors)(k,k) = T(1);
            }

            for (int sweep = 0; sweep < JACOBI_MAX_SWEEPS; ++sweep)
            {
                bool rotated = false;
                for (index_t p = 1; p < n; ++p)
                    for (index_t q = p + 1; q <= n; ++q)
                    {
                        const T apq = a(p,q);
                        if (!(std::abs(apq) > eps * std::sqrt(std::abs(a(p,p))) * std::sqrt(std::abs(a(q,q))))) continue;
                        const jacobi_rotation<T> r(a(p,p), a(q,q), apq);

                        // A = J^T A J, for the rotation J in the (p, q) plane: the pair (p, q)
                        // becomes zero, and rows and columns p and q are rotated
                        a(p,p) -= r.t * apq;
                        a(q,q) += r.t * apq;
                        a(p,q) = a(q,p) = T(0);
                        for (index_t i = 1; i <= n; ++i)
                        {
                            if (i == p || i == q) continue;
                            const T g = a(i,p), h = a(i,q);
                            a(i,p) = a(p,i) = r.c * g - r.s * h;
                            a(i,q) = a(q,i) = r.s * g + r.c * h;
                        }
                        if (vectors) rotate_columns(*vectors, p, q, r);
                        rotated = true;
                    }
                if (!rotated) break;
            }

            // Ascending order (selection sort: at most n - 1 column swaps)
            for (index_t k = 1; k <= n; ++k) values(k,1) = a(k,k);
            for (index_t k = 1; k < n; ++k)
            {
                index_t smallest = k;
                for (index_t i = k + 1; i <= n; ++i)
                    if (values(i,1) < values(smallest,1)) smallest = i;
                if (smallest != k)
                {
                    std::swap(values(k,1), values(smallest,1));
                    if (vectors) swap_columns(*vectors, k, smallest);
                }
            }
        }
    }


   /***********************************************************************************************
    * class BasicSymmetricEigen<T, N, ORDER>
    *
    * Eigenvalues and orthonormal eigenvectors of a symmetric matrix (see the description above).
    * Made by eig_sym, and usually named through:
    *
    *     SymmetricEigen<T, ORDER>           dynamic results (any n)
    *     FixedSymmetricEigen<T, N, ORDER>   fixed-size results: no heap memory
    *
    * Template Parameters:
    *     T         Floating-point element type
    *     N         Order of the matrix, or 0 for a dynamic one
    *     ORDER     Storage order of the results
    ***********************************************************************************************/
    template<class T, index_t N, StorageOrder ORDER>
    class BasicSymmetricEigen
    {
        static_assert(std::is_floating_point_v<T>, "SymmetricEigen: the element type must be a floating-point type");

      public:
        typedef T                                         dtype;
        typedef detail::sized_matrix_t<T, N, N, ORDER>    matrix_type;     // n x n
        typedef detail::sized_matrix_t<T, N, 1, ORDER>    vector_type;     // n x 1

      private:
        // The rotations work on columns, so they run on column-major copies
        typedef detail::sized_matrix_t<T, N, N, StorageOrder::ColMajor> work_type;

        vector_type values_;
        matrix_type vectors_;

      public:
        // Decomposes the symmetric matrix a (only its lower triangle is read)
        explicit BasicSymmetricEigen(matrix_type const& a)
            : values_(detail::make_sized<vector_type>(a.rows(), 1)),
              vectors_(detail::make_sized<matrix_type>(a.rows(), a.rows()))
        {
            detail::check_square("eig_sym", a.rows(), a.cols());
            work_type work(a);
            work_type vectors = detail::make_sized<work_type>(a.rows(), a.rows());
            detail::symmetric_eigen(work, this->values_, &vectors);
            this->vectors_ = vectors;
        }

        // Order of the matrix
        index_t size() const {return this->values_.rows();}

        // Eigenvalues in ascending order (n x 1)
        vector_type const& values() const {return this->values_;}

        // Orthonormal eigenvectors (n x n): column k belongs to values()(k)
        matrix_type const& vectors() const {return this->vectors_;}
    };

    // Eigen-decomposition of a symmetric matrix with dynamic results
    template<class T, StorageOrder ORDER>
    using SymmetricEigen = BasicSymmetricEigen<T, 0, ORDER>;

    // Eigen-decomposition of an N x N fixed-size symmetric matrix: fixed-size results, no heap memory
    template<class T, index_t N, StorageOrder ORDER>
    using FixedSymmetricEigen = BasicSymmetricEigen<T, N, ORDER>;


    // Eigenvalues and eigenvectors of a symmetric matrix, view or expression (only its lower
    // triangle is read), in a floating-point type (double for integer matrices) and the storage
    // order of a; fixed-size if a is
    template<MatrixAlgebraObject A>
    auto eig_sym(A const& a)
    {
        typedef detail::solver_t<detail::dtype_t<A>> T;
        constexpr StorageOrder order = detail::solution_order_v<A>;
        if constexpr (detail::traits<A>::fixed_size)
        {
            constexpr index_t N = detail::traits<A>::rows;
            static_assert(detail::traits<A>::cols == N, "eig_sym: the matrix must be square");
            return FixedSymmetricEigen<T, N, order>(FixedSizeMatrix<T, N, N, order>(a));
        }
        else
            return SymmetricEigen<T, order>(DynamicMatrix<T, order>(a));
    }

    // Eigenvalues of a symmetric matrix, view or expression (only its lower triangle is read), in
    // ascending order (n x 1; fixed-size if a is). The same values as eig_sym(a).values().
    template<MatrixAlgebraObject A>
    auto eigenvalues_sym(A const& a)
    {
        typedef detail::solver_t<detail::dtype_t<A>> T;
        constexpr StorageOrder order = detail::solution_order_v<A>;
        constexpr index_t N = detail::traits<A>::fixed_size ? detail::traits<A>::rows : 0;
        static_assert(!detail::traits<A>::fixed_size || detail::traits<A>::cols == N, "eigenvalues_sym: the matrix must be square");
        typedef detail::sized_matrix_t<T, N, N, StorageOrder::ColMajor> MAT;      // as eig_sym
        MAT work(a);
        detail::check_square("eigenvalues_sym", work.rows(), work.cols());
        auto values = detail::make_sized<detail::sized_matrix_t<T, N, 1, order>>(work.rows(), 1);
        detail::symmetric_eigen(work, values, static_cast<MAT*>(nullptr));
        return values;
    }
}
#endif
