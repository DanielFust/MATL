/**************************************************************************************************
 * File: _matl_cholesky_.hpp
 *
 * Description: Cholesky factorization, A = L L^T, for symmetric positive definite matrices.
 *
 *     auto C = cholesky(A);           factorizes a copy of A (A: matrix, view or expression)
 *     auto C = cholesky_inplace(A);   factorizes A in its own memory (A: a matrix of a
 *                                     floating-point type), which then holds L
 *     X = C.solve(B);                 solution of A X = B (B: one system per column)
 *     C.solve_inplace(X);             overwrites X (holding B) with the solution
 *     X = cholesky_solve(A, B);       factorizes and solves once
 *
 *     C.rcond()                       estimated reciprocal condition number (1-norm)
 *     C.determinant()                 determinant of A
 *     C.L(), C.factor()               the factor L (a copy, and the stored factor)
 *
 * Only the lower triangle of A (with the diagonal) is read: A is taken to be symmetric, as in
 * LAPACK's xPOTRF with uplo = 'L'. The factor is stored with zeros above the diagonal, so after
 * cholesky_inplace(A), A is L.
 *
 * Cholesky needs no pivoting and half the arithmetic of LU, but only works for symmetric positive
 * definite matrices: a pivot that is not positive means A is not positive definite (or is too
 * close to singular to tell), and is reported as matl::not_positive_definite_error.
 *
 * Algorithm: right-looking, blocked like the LU factorization (see _matl_lu_.hpp): each panel's
 * columns of L are computed column by column, and the rest of the matrix is updated with
 * A22 -= L21 * L21^T on the product kernel, one product per block column of the lower triangle
 * (the upper triangle is not needed, so this does half the arithmetic of the full update).
 **************************************************************************************************/
#ifndef MATL_CHOLESKY_HPP
#define MATL_CHOLESKY_HPP
namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * void cholesky_factorize(DynamicMatrix<T, ORDER>& a)
        *
        * Overwrites the square matrix a with its Cholesky factor L (zeros above the diagonal),
        * reading only the lower triangle of a. Reports matl::not_positive_definite_error for a
        * pivot that is not positive.
        *******************************************************************************************/
        template<class T, StorageOrder ORDER>
        void cholesky_factorize(DynamicMatrix<T, ORDER>& a)
        {
            const index_t n = a.rows();

            for (index_t k0 = 1; k0 <= n; k0 += FACTORIZATION_BLOCK)
            {
                const index_t k1 = std::min(k0 + FACTORIZATION_BLOCK - 1, n);    // last column of the panel

                // Columns k0..k1 of L (the diagonal block and the rows below it):
                //     l(i,j) = (a(i,j) - sum_{p=k0}^{j-1} l(i,p) l(j,p)) / l(j,j)
                // (columns left of the panel were already subtracted by the earlier panels'
                // updates). Row-major matrices form each l(i,j) with a sum along rows i and j;
                // column-major ones subtract whole columns, down column j.
                for (index_t j = k0; j <= k1; ++j)
                {
                    if constexpr (ORDER == StorageOrder::ColMajor)
                        for (index_t p = k0; p < j; ++p)
                        {
                            const T l = a(j,p);
                            for (index_t i = j; i <= n; ++i) a(i,j) -= a(i,p) * l;
                        }

                    T d = a(j,j);
                    if constexpr (ORDER == StorageOrder::RowMajor)
                        for (index_t p = k0; p < j; ++p) d -= a(j,p) * a(j,p);
                    if (!(d > T(0)))
                        raise<not_positive_definite_error>("cholesky: the matrix is not positive definite (pivot " +
                                                           std::to_string(j) + " is not positive)");
                    const T l_jj = std::sqrt(d);
                    a(j,j) = l_jj;

                    for (index_t i = j + 1; i <= n; ++i)
                    {
                        T s = a(i,j);
                        if constexpr (ORDER == StorageOrder::RowMajor)
                            for (index_t p = k0; p < j; ++p) s -= a(i,p) * a(j,p);
                        a(i,j) = s / l_jj;
                    }
                }

                if (k1 == n) break;

                // A22 -= L21 * L21^T, lower triangle only: one product per block column of A22
                // (each block column from its diagonal block down), through views
                for (index_t c0 = k1 + 1; c0 <= n; c0 += FACTORIZATION_BLOCK)
                {
                    const index_t c1 = std::min(c0 + FACTORIZATION_BLOCK - 1, n);
                    a(range(c0, n), range(c0, c1)) -= a(range(c0, n), range(k0, k1)) *
                                                      transpose(a(range(c0, c1), range(k0, k1)));
                }
            }

            // Zeros above the diagonal (the upper triangle was not read, but the updates wrote it)
            for (index_t i = 1; i <= n; ++i)
                for (index_t j = i + 1; j <= n; ++j)
                    a(i,j) = T(0);
        }
    }


   /***********************************************************************************************
    * class CholeskyFactorization<T, ORDER>
    *
    * Cholesky factorization A = L L^T of a symmetric positive definite matrix (see the description
    * above). Made by cholesky and cholesky_inplace.
    *
    * Template Parameters:
    *     T         Floating-point element type
    *     ORDER     Storage order of the factor
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    class CholeskyFactorization
    {
        static_assert(std::is_floating_point_v<T>, "CholeskyFactorization: the element type must be a floating-point type");

      private:
        DynamicMatrix<T, ORDER> l_;         // L, with zeros above the diagonal
        T                       norm_1_;    // 1-norm of A (for rcond)

        // 1-norm of the symmetric matrix whose lower triangle a holds
        static T symmetric_norm_1(DynamicMatrix<T, ORDER> const& a)
        {
            const index_t n = a.rows();
            T largest = 0;
            for (index_t j = 1; j <= n; ++j)
            {
                T s = 0;
                for (index_t i = 1; i < j; ++i)  s += std::abs(a(j,i));    // a(i,j) = a(j,i) above the diagonal
                for (index_t i = j; i <= n; ++i) s += std::abs(a(i,j));
                largest = std::max(largest, s);
            }
            return largest;
        }

      public:
        // Factorizes the square matrix a (its memory may be a matrix's, see cholesky_inplace)
        explicit CholeskyFactorization(DynamicMatrix<T, ORDER> a) : l_(std::move(a)), norm_1_(0)
        {
            detail::check_square("cholesky", this->l_.rows(), this->l_.cols());
            this->norm_1_ = symmetric_norm_1(this->l_);
            detail::cholesky_factorize(this->l_);
        }

        // Order of the system (A is size() x size())
        index_t size() const {return this->l_.rows();}

        // The stored factor L (zeros above the diagonal)
        DynamicMatrix<T, ORDER> const& factor() const {return this->l_;}

        // A copy of the factor L
        DynamicMatrix<T, ORDER> L() const {return this->l_;}

        // Overwrites x (a matrix or view holding B, one system per column) with the solution of
        // A X = B: L Y = B, then L^T X = Y
        template<typename X> requires WritableMatrix<X>
        void solve_inplace(X&& x) const
        {
            static_assert(std::is_floating_point_v<detail::dtype_t<X>>,
                          "CholeskyFactorization::solve_inplace: the solution must have a floating-point type");
            detail::check_right_hand_side("cholesky solve", this->size(), x.rows(), x.cols());
            const auto lt = [this](index_t i, index_t p) {return this->l_(p,i);};
            detail::triangular_solve<true, false>(this->l_, this->size(), x, "cholesky solve");
            detail::triangular_solve<false, false>(lt, this->size(), x, "cholesky solve");
        }

        // Solution of A X = B (B: matrix, view or expression; one system per column)
        template<MatrixAlgebraObject B>
        auto solve(B const& b) const
        {
            DynamicMatrix<detail::solver_t<detail::promote_t<T, detail::dtype_t<B>>>, detail::solution_order_v<B>> x(b);
            this->solve_inplace(x);
            return x;
        }

        // Determinant of A: the square of the product of L's diagonal
        T determinant() const
        {
            T d = 1;
            for (index_t k = 1; k <= this->size(); ++k) d *= this->l_(k,k);
            return d * d;
        }

        // Estimated reciprocal condition number in the 1-norm (see LUFactorization::rcond). A is
        // symmetric, so a solve with A^T is a solve with A.
        T rcond() const
        {
            const auto solve = [this](auto& x) {this->solve_inplace(x);};
            return detail::reciprocal_condition(this->norm_1_,
                                                detail::estimate_inverse_norm_1<T>(this->size(), solve, solve),
                                                this->size());
        }
    };


    // Cholesky factorization of a copy of a (a symmetric positive definite matrix, view or
    // expression; only its lower triangle is read), in a floating-point type (double for integer
    // matrices) and the storage order of a
    template<MatrixAlgebraObject A>
    auto cholesky(A const& a)
    {
        typedef detail::solver_t<detail::dtype_t<A>> T;
        constexpr StorageOrder order = detail::solution_order_v<A>;
        return CholeskyFactorization<T, order>(DynamicMatrix<T, order>(a));
    }

    // Cholesky factorization in the memory of a (a square matrix of a floating-point type): a then
    // holds L (zeros above the diagonal), and must outlive the factorization
    template<class T, StorageOrder ORDER>
    CholeskyFactorization<T, ORDER> cholesky_inplace(DynamicMatrix<T, ORDER>& a)
    {return CholeskyFactorization<T, ORDER>(DynamicMatrix<T, ORDER>(a.rows(), a.cols(), a.data()));}

    template<class T, index_t M, index_t N, StorageOrder ORDER>
    CholeskyFactorization<T, ORDER> cholesky_inplace(FixedSizeMatrix<T, M, N, ORDER>& a)
    {return CholeskyFactorization<T, ORDER>(DynamicMatrix<T, ORDER>(M, N, a.data()));}

    // Solution of A X = B by Cholesky factorization (A symmetric positive definite; only its
    // lower triangle is read)
    template<MatrixAlgebraObject A, MatrixAlgebraObject B>
    auto cholesky_solve(A const& a, B const& b)
    {return cholesky(a).solve(b);}
}
#endif
