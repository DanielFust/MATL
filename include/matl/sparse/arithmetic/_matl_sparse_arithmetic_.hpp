/**************************************************************************************************
 * File: _matl_sparse_arithmetic_.hpp
 *
 * Description: Arithmetic of sparse matrices, with sparse results:
 *
 *     S + T, S - T                sum and difference (the union of the two patterns)
 *     S * T                       sparse matrix product (Gustavson's row-by-row algorithm)
 *     CSRMatrix<T> St = transpose(S);   the transpose, formed explicitly
 *     diagonal(S)                 the diagonal, as a dense n x 1 matrix
 *
 * Unlike dense expressions, these compute their result when they are called: the pattern of the
 * result must be known before its values can be stored. Entries that cancel to 0 are kept (the
 * pattern of a sum is the union of the patterns, and of a product all the (i,j) reachable through
 * a stored S(i,k) and T(k,j)). The dimensions are always checked (matl::dimension_error).
 **************************************************************************************************/
#ifndef MATL_SPARSE_ARITHMETIC_HPP
#define MATL_SPARSE_ARITHMETIC_HPP
namespace matl
{
    namespace detail
    {
        // R = S + sign * T, row by row: the two rows' sorted column lists are merged
        template<class T1, class T2, class I>
        CSRMatrix<promote_t<T1, T2>, I> sparse_sum(CSRMatrix<T1, I> const& s, CSRMatrix<T2, I> const& t, int sign,
                                                   char const* operation)
        {
            typedef promote_t<T1, T2> R;
            if ((s.rows() != t.rows()) || (s.cols() != t.cols()))
                incompatible_dimensions(operation, s.rows(), s.cols(), t.rows(), t.cols());

            const auto so = s.row_offsets(), sc = s.column_indices();
            const auto sv = s.values();
            const auto to = t.row_offsets(), tc = t.column_indices();
            const auto tv = t.values();
            std::vector<I> offsets(s.rows() + 1, I(0)), columns;
            std::vector<R> values;
            columns.reserve(s.nnz() + t.nnz());
            values.reserve(s.nnz() + t.nnz());
            for (index_t i = 0; i < s.rows(); ++i)
            {
                I p = so[i], q = to[i];
                while (p < so[i + 1] || q < to[i + 1])
                {
                    if (q == to[i + 1] || (p < so[i + 1] && sc[p] < tc[q]))
                        {columns.push_back(sc[p]); values.push_back(R(sv[p])); ++p;}
                    else if (p == so[i + 1] || tc[q] < sc[p])
                        {columns.push_back(tc[q]); values.push_back(R(sign) * R(tv[q])); ++q;}
                    else
                        {columns.push_back(sc[p]); values.push_back(R(sv[p]) + R(sign) * R(tv[q])); ++p; ++q;}
                }
                check_fits_index<I>(values.size(), "the number of entries");
                offsets[i + 1] = static_cast<I>(values.size());
            }
            return CSRMatrix<R, I>(s.rows(), s.cols(), std::move(offsets), std::move(columns), std::move(values));
        }
    }

    // Sum of two sparse matrices (the union of their patterns)
    template<class T1, class T2, class I>
    CSRMatrix<detail::promote_t<T1, T2>, I> operator+ (CSRMatrix<T1, I> const& s, CSRMatrix<T2, I> const& t)
    {return detail::sparse_sum(s, t, 1, "sparse matrix sum");}

    // Difference of two sparse matrices (the union of their patterns)
    template<class T1, class T2, class I>
    CSRMatrix<detail::promote_t<T1, T2>, I> operator- (CSRMatrix<T1, I> const& s, CSRMatrix<T2, I> const& t)
    {return detail::sparse_sum(s, t, -1, "sparse matrix difference");}

   /***********************************************************************************************
    * CSRMatrix operator* (CSRMatrix const& s, CSRMatrix const& t)
    *
    * Sparse matrix product by Gustavson's algorithm: row i of S T is the sum of the rows k of T
    * scaled by the stored S(i,k). The row is accumulated in a dense array of T's columns, with a
    * marker of which columns it has, and then gathered in increasing column order. Cost: one
    * multiply-add per pair of stored S(i,k) and T(k,j), plus sorting each result row's columns.
    ***********************************************************************************************/
    template<class T1, class T2, class I>
    CSRMatrix<detail::promote_t<T1, T2>, I> operator* (CSRMatrix<T1, I> const& s, CSRMatrix<T2, I> const& t)
    {
        typedef detail::promote_t<T1, T2> R;
        if (s.cols() != t.rows())
            detail::incompatible_dimensions("sparse matrix product (columns of the left operand must equal rows of the right operand)",
                                            s.rows(), s.cols(), t.rows(), t.cols());

        const auto so = s.row_offsets(), sc = s.column_indices();
        const auto sv = s.values();
        const auto to = t.row_offsets(), tc = t.column_indices();
        const auto tv = t.values();
        const index_t n = t.cols();
        std::vector<R>       accumulator(n, R(0));
        std::vector<index_t> marker(n, static_cast<index_t>(-1));   // the row that last touched column j
        std::vector<I>       row_columns;                            // columns of the current result row

        std::vector<I> offsets(s.rows() + 1, I(0)), columns;
        std::vector<R> values;
        for (index_t i = 0; i < s.rows(); ++i)
        {
            row_columns.clear();
            for (I p = so[i]; p < so[i + 1]; ++p)
            {
                const index_t k = static_cast<index_t>(sc[p]);
                const R a = R(sv[p]);
                for (I q = to[k]; q < to[k + 1]; ++q)
                {
                    const index_t j = static_cast<index_t>(tc[q]);
                    if (marker[j] != i) {marker[j] = i; accumulator[j] = R(0); row_columns.push_back(tc[q]);}
                    accumulator[j] += a * R(tv[q]);
                }
            }
            std::sort(row_columns.begin(), row_columns.end());
            for (I j : row_columns)
            {
                columns.push_back(j);
                values.push_back(accumulator[static_cast<index_t>(j)]);
            }
            detail::check_fits_index<I>(values.size(), "the number of entries");
            offsets[i + 1] = static_cast<I>(values.size());
        }
        return CSRMatrix<R, I>(s.rows(), n, std::move(offsets), std::move(columns), std::move(values));
    }

    // The diagonal of a sparse matrix, as a dense min(m, n) x 1 matrix
    template<class T, class I>
    DynamicMatrix<T, StorageOrder::ColMajor> diagonal(CSRMatrix<T, I> const& s)
    {
        const index_t k = std::min(s.rows(), s.cols());
        DynamicMatrix<T, StorageOrder::ColMajor> d(k, 1);
        for (index_t i = 1; i <= k; ++i) d(i) = s(i, i);
        return d;
    }

    namespace detail
    {
        // The transpose of a sparse matrix, formed explicitly (a counting sort by column)
        template<class T, class I>
        CSRMatrix<T, I> transpose_csr(CSRMatrix<T, I> const& s)
        {
            const auto so = s.row_offsets(), sc = s.column_indices();
            const auto sv = s.values();
            const index_t m = s.rows(), n = s.cols();
            std::vector<I> offsets(n + 1, I(0));
            for (index_t p = 0; p < s.nnz(); ++p) ++offsets[static_cast<index_t>(sc[p]) + 1];
            for (index_t j = 0; j < n; ++j) offsets[j + 1] += offsets[j];
            std::vector<I> columns(s.nnz());
            std::vector<T> values(s.nnz());
            std::vector<I> next(offsets.begin(), offsets.end() - 1);
            for (index_t i = 0; i < m; ++i)                            // rows in order: columns end up sorted
                for (I p = so[i]; p < so[i + 1]; ++p)
                {
                    const I q = next[static_cast<index_t>(sc[p])]++;
                    columns[q] = static_cast<I>(i);
                    values[q]  = sv[p];
                }
            return CSRMatrix<T, I>(n, m, std::move(offsets), std::move(columns), std::move(values));
        }
    }

    // CSRMatrix<double> St = transpose(S): the transpose, formed explicitly
    template<class T, class I>
    CSRMatrix<T, I>::CSRMatrix(expr::SparseTransposed<T, I> const& transposed) :
        CSRMatrix(detail::transpose_csr(transposed.matrix()))
    {}

    // The dense matrix of transpose(S)
    template<class T, class I>
    DynamicMatrix<T, StorageOrder::RowMajor> dense(expr::SparseTransposed<T, I> const& transposed)
    {return dense(detail::transpose_csr(transposed.matrix()));}
}
#endif
