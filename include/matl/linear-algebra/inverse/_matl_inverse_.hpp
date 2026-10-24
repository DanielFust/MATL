/**************************************************************************************************
 * File: _matl_inverse_.hpp
 *
 * Description: Determinants and inverses of square matrices.
 *
 *     double d = det(A);              determinant (0 for a singular matrix)
 *     auto X = inv(A);                inverse (matl::singular_matrix_error if A is singular)
 *
 * Fixed-size matrices of order 1, 2 and 3 use closed forms (the cofactor expansion, and the
 * adjugate divided by the determinant), e.g. for the Jacobian of a finite-element mapping at each
 * integration point. Other matrices use the LU factorization (partial pivoting). Neither needs
 * heap memory for a fixed-size matrix, whose inverse is fixed-size.
 *
 * To solve a system, prefer lu(A).solve(B) to inv(A) * B: it is about three times cheaper and
 * more accurate. A is singular for inv if its determinant (closed forms) or a pivot (LU) is
 * exactly zero; a nearly singular matrix has an inaccurate inverse (lu(A).rcond() tells how
 * close A is to singular).
 *
 * Integer matrices give double results; float stays float.
 **************************************************************************************************/
#ifndef MATL_INVERSE_HPP
#define MATL_INVERSE_HPP
namespace matl
{
    namespace detail
    {
        // Determinant of a fixed-size matrix of order 1, 2 or 3 by cofactor expansion
        template<class T, index_t N, StorageOrder ORDER>
        T closed_form_determinant(FixedSizeMatrix<T, N, N, ORDER> const& a)
        {
            const auto x = [&a](index_t i, index_t j) {return a(i,j);};
            if constexpr (N == 1) return x(1,1);
            else if constexpr (N == 2) return x(1,1) * x(2,2) - x(1,2) * x(2,1);
            else return x(1,1) * (x(2,2) * x(3,3) - x(2,3) * x(3,2))
                      - x(1,2) * (x(2,1) * x(3,3) - x(2,3) * x(3,1))
                      + x(1,3) * (x(2,1) * x(3,2) - x(2,2) * x(3,1));
        }

        // Determinant from packed LU factors (U's diagonal, with the sign of the permutation)
        template<class MAT, class PIVOTS>
        dtype_t<MAT> lu_determinant(MAT const& lu, PIVOTS const& pivots)
        {
            dtype_t<MAT> d = 1;
            for (index_t k = 1; k <= lu.rows(); ++k)
            {
                d *= lu(k,k);
                if (pivots[k - 1] != k) d = -d;
            }
            return d;
        }
    }

    // Determinant of a square matrix, view or expression (0 if it is singular)
    template<MatrixAlgebraObject A>
    auto det(A const& a)
    {
        typedef detail::solver_t<detail::dtype_t<A>> T;
        constexpr StorageOrder order = detail::solution_order_v<A>;
        if constexpr (detail::traits<A>::fixed_size)
        {
            constexpr index_t N = detail::traits<A>::rows;
            static_assert(detail::traits<A>::cols == N, "det: the matrix must be square");
            if constexpr (N <= 3) return detail::closed_form_determinant<T>(FixedSizeMatrix<T, N, N, order>(a));
            else
            {
                FixedSizeMatrix<T, N, N, order> lu(a);
                std::array<index_t, N> pivots;
                index_t zero_pivot = 0;
                detail::lu_factorize(lu, pivots, &zero_pivot);
                return zero_pivot ? T(0) : detail::lu_determinant(lu, pivots);
            }
        }
        else
        {
            detail::check_square("det", a.rows(), a.cols());
            DynamicMatrix<T, order> lu(a);
            std::vector<index_t> pivots;
            index_t zero_pivot = 0;
            detail::lu_factorize(lu, pivots, &zero_pivot);
            return zero_pivot ? T(0) : detail::lu_determinant(lu, pivots);
        }
    }

    // Inverse of a square matrix, view or expression: fixed-size if a is, in the storage order of
    // a. Reports matl::singular_matrix_error for a singular matrix.
    template<MatrixAlgebraObject A>
    auto inv(A const& a)
    {
        typedef detail::solver_t<detail::dtype_t<A>> T;
        constexpr StorageOrder order = detail::solution_order_v<A>;
        if constexpr (detail::traits<A>::fixed_size && detail::traits<A>::rows <= 3)
        {
            constexpr index_t N = detail::traits<A>::rows;
            static_assert(detail::traits<A>::cols == N, "inv: the matrix must be square");
            const FixedSizeMatrix<T, N, N, order> x(a);
            const T d = detail::closed_form_determinant<T>(x);
            if (d == T(0)) detail::raise<singular_matrix_error>("inv: the matrix is singular (its determinant is 0)");
            FixedSizeMatrix<T, N, N, order> y;
            if constexpr (N == 1) y(1,1) = T(1) / d;
            else if constexpr (N == 2)
            {
                y(1,1) =  x(2,2) / d;  y(1,2) = -x(1,2) / d;
                y(2,1) = -x(2,1) / d;  y(2,2) =  x(1,1) / d;
            }
            else
            {
                // The adjugate (transposed cofactors) divided by the determinant
                const auto cofactor = [&x](index_t i1, index_t i2, index_t j1, index_t j2)
                {return x(i1,j1) * x(i2,j2) - x(i1,j2) * x(i2,j1);};
                y(1,1) =  cofactor(2, 3, 2, 3) / d;  y(1,2) = -cofactor(1, 3, 2, 3) / d;  y(1,3) =  cofactor(1, 2, 2, 3) / d;
                y(2,1) = -cofactor(2, 3, 1, 3) / d;  y(2,2) =  cofactor(1, 3, 1, 3) / d;  y(2,3) = -cofactor(1, 2, 1, 3) / d;
                y(3,1) =  cofactor(2, 3, 1, 2) / d;  y(3,2) = -cofactor(1, 3, 1, 2) / d;  y(3,3) =  cofactor(1, 2, 1, 2) / d;
            }
            return y;
        }
        else
        {
            const auto F = lu(a);                                 // checks that a is square
            typedef std::remove_cvref_t<decltype(F.factors())> MAT;
            MAT y = detail::factor_storage<MAT>::square(F.size());
            y = T(0);
            for (index_t k = 1; k <= F.size(); ++k) y(k,k) = T(1);
            F.solve_inplace(y);
            return y;
        }
    }
}
#endif
