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
 * (products inside them are computed into temporaries, as in any evaluation). Large reductions run
 * on several threads with OpenMP, and give the same result for any number of threads (see
 * detail::reduce_elements).
 **************************************************************************************************/
#ifndef MATL_REDUCTIONS_HPP
#define MATL_REDUCTIONS_HPP
namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * T blocked_reduce(n, per_block, partial, combine, identity, work)
        *
        * A reduction over n items in fixed blocks: partial(first, last) reduces the items
        * [first, last) of one block (from the identity), and the blocks' results are combined in
        * order, combine(...combine(combine(identity, b1), b2)..., bk). The blocks are reduced in
        * parallel when the work is large enough (see parallel_ranges); the result does not
        * depend on the number of threads. The blocks' results are kept in the temporary arena,
        * so a reduction makes no heap allocations.
        *******************************************************************************************/
        template<typename T, typename PARTIAL, typename COMBINE>
        T blocked_reduce(index_t n, index_t per_block, PARTIAL const& partial, COMBINE const& combine,
                         T identity, index_t work)
        {
            static_assert(std::is_trivially_copyable_v<T> && alignof(T) <= TemporaryArena::ALIGNMENT,
                          "blocked_reduce: partial results are kept in the temporary arena");
            if (n == 0) return identity;
            per_block = std::max<index_t>(per_block, 1);
            const index_t blocks = (n + per_block - 1) / per_block;
            if (blocks == 1) return combine(identity, partial(index_t(0), n));

            TemporaryBlock storage(blocks * sizeof(T));
            T* const results = reinterpret_cast<T*>(storage.data());
            parallel_ranges(blocks, work, [&](index_t first, index_t last)
            {
                for (index_t b = first; b < last; ++b)
                    results[b] = partial(b * per_block, std::min(n, (b + 1) * per_block));
            });
            T result = identity;
            for (index_t b = 0; b < blocks; ++b) result = combine(result, results[b]);
            return result;
        }

       /*******************************************************************************************
        * T reduce_elements(mao, identity, accumulate, combine)
        *
        * Reduces the elements of a matrix or expression: accumulate(a, x) adds element x to a
        * partial result, and combine(a, b) joins two partial results. The elements are visited in
        * storage order (with a single index if the expression allows one, see linear_index_v;
        * otherwise by (row, column), row by row when the operands have different storage orders),
        * in blocks of about REDUCTION_BLOCK elements (whole rows or columns for the (row, column)
        * order). The blocks are reduced from the identity, in parallel with OpenMP, and combined
        * in order, so the result does not depend on the number of threads (see
        * _matl_parallel_.hpp).
        *******************************************************************************************/
        template<MatrixAlgebraObject MAO, typename T, typename ACCUMULATE, typename COMBINE>
        T reduce_elements(MAO const& mao, T identity, ACCUMULATE const& accumulate, COMBINE const& combine)
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
                return blocked_reduce<T>(M * N, REDUCTION_BLOCK, [&](index_t first, index_t last)
                {
                    T a = identity;
                    for (index_t k = first + 1; k <= last; ++k) a = accumulate(a, source(k));
                    return a;
                }, combine, identity, M * N);
            }
            else if constexpr (order == StorageOrder::ColMajor)
            {
                return blocked_reduce<T>(N, REDUCTION_BLOCK / std::max<index_t>(M, 1), [&](index_t first, index_t last)
                {
                    T a = identity;
                    for (index_t j = first + 1; j <= last; ++j)
                        for (index_t i = 1; i <= M; ++i) a = accumulate(a, source(i, j));
                    return a;
                }, combine, identity, M * N);
            }
            else
            {
                return blocked_reduce<T>(M, REDUCTION_BLOCK / std::max<index_t>(N, 1), [&](index_t first, index_t last)
                {
                    T a = identity;
                    for (index_t i = first + 1; i <= last; ++i)
                        for (index_t j = 1; j <= N; ++j) a = accumulate(a, source(i, j));
                    return a;
                }, combine, identity, M * N);
            }
        }

        // Partial result of min and max: the extreme value so far, if any element has been seen
        // (the first element starts the search, as in a serial loop)
        template<typename T>
        struct extreme
        {
            bool any;
            T    value;
        };

        // min (LESS = true) or max of the elements
        template<bool LESS, MatrixAlgebraObject MAO>
        dtype_t<MAO> extreme_element(MAO const& mao)
        {
            typedef dtype_t<MAO> T;
            const auto better = [](T x, T than) {if constexpr (LESS) return x < than; else return than < x;};
            const auto accumulate = [better](extreme<T> a, T x) {return (!a.any || better(x, a.value)) ? extreme<T>{true, x} : a;};
            const auto combine = [better](extreme<T> a, extreme<T> b)
            {
                if (!b.any) return a;
                if (!a.any) return b;
                return better(b.value, a.value) ? b : a;
            };
            return reduce_elements(mao, extreme<T>{false, T{}}, accumulate, combine).value;
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
        typedef detail::dtype_t<MAO> T;
        return detail::reduce_elements(mao, T(0), [](T a, T x) {return T(a + x);}, [](T a, T b) {return T(a + b);});
    }

    // Product of all the elements (1 for an empty matrix)
    template<MatrixAlgebraObject MAO>
    detail::dtype_t<MAO> prod(MAO const& mao)
    {
        typedef detail::dtype_t<MAO> T;
        return detail::reduce_elements(mao, T(1), [](T a, T x) {return T(a * x);}, [](T a, T b) {return T(a * b);});
    }

    // Mean of all the elements (double for integer matrices)
    template<MatrixAlgebraObject MAO>
    detail::real_t<detail::dtype_t<MAO>> mean(MAO const& mao)
    {
        typedef detail::real_t<detail::dtype_t<MAO>> real;
        if (mao.size() == 0) detail::empty_reduction("mean");
        const real s = detail::reduce_elements(mao, real(0), [](real a, auto x) {return a + static_cast<real>(x);},
                                               [](real a, real b) {return a + b;});
        return s / static_cast<real>(mao.size());
    }

    // Smallest element
    template<MatrixAlgebraObject MAO>
    detail::dtype_t<MAO> min(MAO const& mao)
    {
        if (mao.size() == 0) detail::empty_reduction("min");
        return detail::extreme_element<true>(mao);
    }

    // Largest element
    template<MatrixAlgebraObject MAO>
    detail::dtype_t<MAO> max(MAO const& mao)
    {
        if (mao.size() == 0) detail::empty_reduction("max");
        return detail::extreme_element<false>(mao);
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
        const real s = detail::reduce_elements(mao, real(0),
                                               [](real a, auto x) {return a + static_cast<real>(x) * static_cast<real>(x);},
                                               [](real a, real b) {return a + b;});
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
            // The running sums are kept in the temporary arena (no heap allocations)
            detail::TemporaryBlock storage(N * sizeof(dtype));
            dtype* const sums = reinterpret_cast<dtype*>(storage.data());
            std::fill(sums, sums + N, dtype(0));
            for (index_t i = 1; i <= M; ++i)
                for (index_t j = 1; j <= N; ++j) sums[j - 1] += detail::magnitude(source(i, j));
            for (index_t j = 0; j < N; ++j) largest = std::max(largest, sums[j]);
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
            // The running sums are kept in the temporary arena (no heap allocations)
            detail::TemporaryBlock storage(M * sizeof(dtype));
            dtype* const sums = reinterpret_cast<dtype*>(storage.data());
            std::fill(sums, sums + M, dtype(0));
            for (index_t j = 1; j <= N; ++j)
                for (index_t i = 1; i <= M; ++i) sums[i - 1] += detail::magnitude(source(i, j));
            for (index_t i = 0; i < M; ++i) largest = std::max(largest, sums[i]);
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
