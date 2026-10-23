/**************************************************************************************************
 * File: _matl_reductions_.hpp
 *
 * Description: Reductions: scalars computed from all the elements of a matrix or expression.
 *
 *     sum(A), prod(A), mean(A)    sum, product and mean of all the elements
 *     min(A), max(A)              smallest and largest element
 *     trace(A)                    sum of the diagonal (A must be square)
 *     dot(A, B)                   sum of the elementwise products (the Frobenius inner product;
 *                                 for vectors, the dot product)
 *     norm_fro(A)                 Frobenius norm, sqrt(dot(A, A)) (for vectors, the 2-norm)
 *     norm_1(A)                   largest column sum of absolute values (for a column vector,
 *                                 the 1-norm)
 *     norm_inf(A)                 largest row sum of absolute values (for a column vector, the
 *                                 largest absolute value)
 *
 * Note: unlike MATLAB, where sum(A) of a matrix gives the column sums, every reduction here gives
 * one scalar for the whole matrix.
 *
 * Expressions are reduced element by element without evaluating them into a temporary matrix
 * (products inside them are computed into temporaries, as in any evaluation).
 **************************************************************************************************/
#ifndef MATL_REDUCTIONS_HPP
#define MATL_REDUCTIONS_HPP
namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * void for_each_element(mao, f)
        *
        * Calls f(value) for every element of a matrix or expression. When the expression has one
        * storage order, elements are visited in storage order: with a single index if it allows
        * one (see linear_index_v), otherwise (views) by (row, column). When its operands have
        * different storage orders (e.g. A + transpose(B)), elements are visited row by row,
        * because a single index then refers to different elements of the operands.
        *******************************************************************************************/
        template<MatrixAlgebraObject MAO, typename F>
        void for_each_element(MAO const& mao, F&& f)
        {
            #ifdef MATL_DIMENSION_CHECK
            check_dimensions(mao);
            #endif
            const evaluator<MAO> source(mao);
            const index_t M = mao.rows();
            const index_t N = mao.cols();
            constexpr StorageOrder order = traits<MAO>::storage_order;
            if constexpr ((order == StorageOrder::RowMajor || order == StorageOrder::ColMajor) && linear_index_v<MAO>)
            {
                const index_t K = mao.size();
                for (index_t k = 1; k <= K; ++k) f(source(k));
            }
            else if constexpr (order == StorageOrder::ColMajor)
            {
                for (index_t j = 1; j <= N; ++j)
                    for (index_t i = 1; i <= M; ++i)
                        f(source(i, j));
            }
            else
            {
                for (index_t i = 1; i <= M; ++i)
                    for (index_t j = 1; j <= N; ++j)
                        f(source(i, j));
            }
        }

        // Floating-point type for results that are not integers (mean, Frobenius norm)
        template<typename T>
        using real_t = std::conditional_t<std::is_integral_v<T>, double, T>;

        // Absolute value (also for unsigned types, which std::abs does not accept)
        template<typename T>
        T magnitude(T x)
        {
            if constexpr (std::is_unsigned_v<T>) return x;
            else                                 return (x < T(0)) ? T(-x) : x;
        }

        // Reports a reduction of an empty matrix that has no value
        [[noreturn]] inline void empty_reduction(char const* name)
        {raise<dimension_error>(std::string(name) + " of an empty matrix is not defined");}
    }

    // Sum of all the elements (0 for an empty matrix)
    template<MatrixAlgebraObject MAO>
    detail::dtype_t<MAO> sum(MAO const& mao)
    {
        detail::dtype_t<MAO> s = 0;
        detail::for_each_element(mao, [&](auto x) {s += x;});
        return s;
    }

    // Product of all the elements (1 for an empty matrix)
    template<MatrixAlgebraObject MAO>
    detail::dtype_t<MAO> prod(MAO const& mao)
    {
        detail::dtype_t<MAO> p = 1;
        detail::for_each_element(mao, [&](auto x) {p *= x;});
        return p;
    }

    // Mean of all the elements (double for integer matrices)
    template<MatrixAlgebraObject MAO>
    detail::real_t<detail::dtype_t<MAO>> mean(MAO const& mao)
    {
        typedef detail::real_t<detail::dtype_t<MAO>> real;
        if (mao.size() == 0) detail::empty_reduction("mean");
        real s = 0;
        detail::for_each_element(mao, [&](auto x) {s += static_cast<real>(x);});
        return s / static_cast<real>(mao.size());
    }

    // Smallest element
    template<MatrixAlgebraObject MAO>
    detail::dtype_t<MAO> min(MAO const& mao)
    {
        if (mao.size() == 0) detail::empty_reduction("min");
        detail::dtype_t<MAO> m{};
        bool first = true;
        detail::for_each_element(mao, [&](auto x) {if (first || x < m) m = x; first = false;});
        return m;
    }

    // Largest element
    template<MatrixAlgebraObject MAO>
    detail::dtype_t<MAO> max(MAO const& mao)
    {
        if (mao.size() == 0) detail::empty_reduction("max");
        detail::dtype_t<MAO> m{};
        bool first = true;
        detail::for_each_element(mao, [&](auto x) {if (first || m < x) m = x; first = false;});
        return m;
    }

    // Sum of the diagonal elements (the matrix must be square)
    template<MatrixAlgebraObject MAO>
    detail::dtype_t<MAO> trace(MAO const& mao)
    {
        #ifdef MATL_DIMENSION_CHECK
        if (mao.rows() != mao.cols())
            detail::raise<dimension_error>("trace of a non-square " + detail::dimensions(mao.rows(), mao.cols()) + " matrix");
        detail::check_dimensions(mao);
        #endif
        const detail::evaluator<MAO> source(mao);
        detail::dtype_t<MAO> s = 0;
        for (index_t i = 1; i <= mao.rows(); ++i) s += source(i, i);
        return s;
    }

    // Sum of the elementwise products of A and B (A and B must have the same dimensions)
    template<MatrixAlgebraObject A, MatrixAlgebraObject B> requires SizesMatch<A, B>
    auto dot(A const& a, B const& b)
    {
        return sum(expr::Elementwise<expr::Multiply, A, B>(a, b));
    }

    // Frobenius norm: square root of the sum of the squared elements
    template<MatrixAlgebraObject MAO>
    detail::real_t<detail::dtype_t<MAO>> norm_fro(MAO const& mao)
    {
        typedef detail::real_t<detail::dtype_t<MAO>> real;
        real s = 0;
        detail::for_each_element(mao, [&](auto x) {s += static_cast<real>(x) * static_cast<real>(x);});
        return std::sqrt(s);
    }

    // 1-norm: largest sum of absolute values of a column (0 for an empty matrix). The elements are
    // visited in storage order: a row-major matrix keeps one running sum per column.
    template<MatrixAlgebraObject MAO>
    detail::dtype_t<MAO> norm_1(MAO const& mao)
    {
        #ifdef MATL_DIMENSION_CHECK
        detail::check_dimensions(mao);
        #endif
        typedef detail::dtype_t<MAO> dtype;
        const detail::evaluator<MAO> source(mao);
        const index_t M = mao.rows();
        const index_t N = mao.cols();
        dtype largest = 0;
        if constexpr (detail::traits<MAO>::storage_order == StorageOrder::RowMajor)
        {
            std::vector<dtype> sums(N, dtype(0));
            for (index_t i = 1; i <= M; ++i)
                for (index_t j = 1; j <= N; ++j) sums[j - 1] += detail::magnitude(source(i, j));
            for (dtype s : sums) largest = std::max(largest, s);
        }
        else
        {
            for (index_t j = 1; j <= N; ++j)
            {
                dtype s = 0;
                for (index_t i = 1; i <= M; ++i) s += detail::magnitude(source(i, j));
                largest = std::max(largest, s);
            }
        }
        return largest;
    }

    // Infinity norm: largest sum of absolute values of a row (0 for an empty matrix). The elements
    // are visited in storage order: a column-major matrix keeps one running sum per row.
    template<MatrixAlgebraObject MAO>
    detail::dtype_t<MAO> norm_inf(MAO const& mao)
    {
        #ifdef MATL_DIMENSION_CHECK
        detail::check_dimensions(mao);
        #endif
        typedef detail::dtype_t<MAO> dtype;
        const detail::evaluator<MAO> source(mao);
        const index_t M = mao.rows();
        const index_t N = mao.cols();
        dtype largest = 0;
        if constexpr (detail::traits<MAO>::storage_order == StorageOrder::ColMajor)
        {
            std::vector<dtype> sums(M, dtype(0));
            for (index_t j = 1; j <= N; ++j)
                for (index_t i = 1; i <= M; ++i) sums[i - 1] += detail::magnitude(source(i, j));
            for (dtype s : sums) largest = std::max(largest, s);
        }
        else
        {
            for (index_t i = 1; i <= M; ++i)
            {
                dtype s = 0;
                for (index_t j = 1; j <= N; ++j) s += detail::magnitude(source(i, j));
                largest = std::max(largest, s);
            }
        }
        return largest;
    }
}
#endif
