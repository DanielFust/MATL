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
 * Algorithm: recursive, like the LU factorization (see _matl_lu_.hpp): the matrix is split in
 * halves, A11 = L11 L11^T (recursively), L21 = A21 L11^-T (a recursive triangular solve),
 * A22 -= L21 L21^T on and below the diagonal (a recursive symmetric update, see
 * detail::subtract_lower_gram), and A22 = L22 L22^T (recursively). Almost all the arithmetic is
 * matrix products on the product kernel, which subtract from the matrix without temporaries.
 **************************************************************************************************/
#ifndef MATL_CHOLESKY_HPP
#define MATL_CHOLESKY_HPP
namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * void cholesky_columns(a, j0, j1)
        *
        * Factorizes the diagonal block a(j0:j1, j0:j1) column by column (its columns left of j0
        * have already been subtracted): l(j,j) = sqrt(a(j,j) - sum l(j,p)^2) and
        * l(i,j) = (a(i,j) - sum l(i,p) l(j,p)) / l(j,j) for the rows i of the block below j. The
        * base case of cholesky_recursive.
        *******************************************************************************************/
        template<class T, StorageOrder ORDER>
        void cholesky_columns(DynamicMatrix<T, ORDER>& a, index_t j0, index_t j1)
        {
            for (index_t j = j0; j <= j1; ++j)
            {
                T d = a(j,j);
                for (index_t p = j0; p < j; ++p) d -= a(j,p) * a(j,p);
                if (!(d > T(0)))
                    raise<not_positive_definite_error>("cholesky: the matrix is not positive definite (pivot " +
                                                       std::to_string(j) + " is not positive)");
                const T l_jj = std::sqrt(d);
                a(j,j) = l_jj;
                for (index_t i = j + 1; i <= j1; ++i)
                {
                    T s = a(i,j);
                    for (index_t p = j0; p < j; ++p) s -= a(i,p) * a(j,p);
                    a(i,j) = s / l_jj;
                }
            }
        }

       /*******************************************************************************************
        * void solve_lower_transposed_right(a, l0, l1, r0, r1)
        *
        * Solves X L^T = B in place for the lower triangular L = a(l0:l1, l0:l1) and
        * B = a(r0:r1, l0:l1) (the rows of L below a factorized diagonal block). L is split in
        * halves, X1 L11^T = B1, B2 -= X1 L21^T (a matrix product), X2 L22^T = B2, so that most of
        * the work runs on the product kernel.
        *******************************************************************************************/
        template<class T, StorageOrder ORDER>
        void solve_lower_transposed_right(DynamicMatrix<T, ORDER>& a, index_t l0, index_t l1, index_t r0, index_t r1)
        {
            if (r0 > r1) return;
            const index_t k = l1 - l0 + 1;
            if (k <= RECURSION_BASE)
            {
                // x(i,c) = (b(i,c) - sum_{p<c} x(i,p) l(c,p)) / l(c,c), row by row
                for (index_t i = r0; i <= r1; ++i)
                    for (index_t c = l0; c <= l1; ++c)
                    {
                        T s = a(i,c);
                        for (index_t p = l0; p < c; ++p) s -= a(i,p) * a(c,p);
                        a(i,c) = s / a(c,c);
                    }
                return;
            }
            const index_t lm = l0 + k / 2;                    // first column of the second half
            solve_lower_transposed_right(a, l0, lm - 1, r0, r1);
            auto B2 = a(range(r0, r1), range(lm, l1));
            subtract_product(B2, a(range(r0, r1), range(l0, lm - 1)), transpose(a(range(lm, l1), range(l0, lm - 1))));
            solve_lower_transposed_right(a, lm, l1, r0, r1);
        }

       /*******************************************************************************************
        * void subtract_lower_gram(a, c0, c1, k0, k1)
        *
        * C -= X X^T on and below the diagonal of C = a(c0:c1, c0:c1), with X = a(c0:c1, k0:k1)
        * (the symmetric update of the trailing matrix, as BLAS's xSYRK). C is split in halves:
        * the two diagonal blocks recursively, and the block below them with one matrix product, so
        * that the upper triangle costs no arithmetic except within small diagonal blocks.
        *******************************************************************************************/
        template<class T, StorageOrder ORDER>
        void subtract_lower_gram(DynamicMatrix<T, ORDER>& a, index_t c0, index_t c1, index_t k0, index_t k1)
        {
            if (c0 > c1) return;
            const index_t size = c1 - c0 + 1;
            if (size <= FACTORIZATION_BLOCK)
            {
                // A small diagonal block: the whole block (its upper triangle is zeroed at the end)
                auto C = a(range(c0, c1), range(c0, c1));
                subtract_product(C, a(range(c0, c1), range(k0, k1)), transpose(a(range(c0, c1), range(k0, k1))));
                return;
            }
            const index_t cm = c0 + size / 2;                 // first row and column of the second half
            subtract_lower_gram(a, c0, cm - 1, k0, k1);
            auto C21 = a(range(cm, c1), range(c0, cm - 1));
            subtract_product(C21, a(range(cm, c1), range(k0, k1)), transpose(a(range(c0, cm - 1), range(k0, k1))));
            subtract_lower_gram(a, cm, c1, k0, k1);
        }

       /*******************************************************************************************
        * void cholesky_recursive(a, j0, j1)
        *
        * Factorizes the trailing block a(j0:j1, j0:j1) in halves: A11 = L11 L11^T (recursively),
        * L21 = A21 L11^-T (a recursive triangular solve), A22 -= L21 L21^T (lower triangle only),
        * and A22 = L22 L22^T (recursively). Small blocks are factorized column by column.
        *******************************************************************************************/
        template<class T, StorageOrder ORDER>
        void cholesky_recursive(DynamicMatrix<T, ORDER>& a, index_t j0, index_t j1)
        {
            const index_t size = j1 - j0 + 1;
            if (size <= RECURSION_BASE)
            {
                cholesky_columns(a, j0, j1);
                return;
            }
            const index_t jm = j0 + size / 2;                 // first row and column of the second half
            cholesky_recursive(a, j0, jm - 1);
            solve_lower_transposed_right(a, j0, jm - 1, jm, j1);
            subtract_lower_gram(a, jm, j1, j0, jm - 1);
            cholesky_recursive(a, jm, j1);
        }

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
            if (n > 0) cholesky_recursive(a, 1, n);

            // Zeros above the diagonal (the upper triangle was not read, but the updates of small
            // diagonal blocks wrote it)
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
