/**************************************************************************************************
 * File: _matl_sparse_functions_.hpp
 *
 * Description: Functions of sparse matrices:
 *
 *     dense(S)                    the dense matrix (row-major) with the elements of S
 *     sparse(D), sparse<I>(D)     the sparse matrix of the nonzero elements of a dense matrix,
 *                                 view or expression (index type I, default index_t)
 *     S * s, s * S, S / s         a sparse matrix with every stored value scaled (same pattern)
 *     -S                          the negated sparse matrix
 *     os << S                     the stored entries, one "(i,j)  value" per line
 *
 * Dense factorizations (lu, cholesky, qr and their _solve forms) do not accept sparse matrices:
 * they would fill in the matrix and cost O(n^3). Calling them on a sparse matrix is a compile
 * error that names the alternatives: an iterative solver, or an explicit conversion with dense(S)
 * for a small matrix.
 **************************************************************************************************/
#ifndef MATL_SPARSE_FUNCTIONS_HPP
#define MATL_SPARSE_FUNCTIONS_HPP
namespace matl
{
    /********************************* Conversions ***********************************************/

    // Dense (row-major) matrix with the elements of a sparse matrix
    template<class T, class I>
    DynamicMatrix<T, StorageOrder::RowMajor> dense(CSRMatrix<T, I> const& s)
    {
        DynamicMatrix<T, StorageOrder::RowMajor> d(s.rows(), s.cols());
        d = T(0);
        const auto offsets = s.row_offsets();
        const auto columns = s.column_indices();
        const auto values  = s.values();
        for (index_t i = 1; i <= s.rows(); ++i)
            for (I p = offsets[i - 1]; p < offsets[i]; ++p)
                d(i, static_cast<index_t>(columns[p]) + 1) = values[p];
        return d;
    }

    // Sparse matrix of the nonzero elements of a dense matrix, view or expression (elements equal
    // to 0 are not stored)
    template<class I = index_t, MatrixAlgebraObject X>
    CSRMatrix<detail::dtype_t<X>, I> sparse(X const& x)
    {
        typedef detail::dtype_t<X> T;
        const DynamicMatrix<T, StorageOrder::RowMajor> d(x);
        detail::check_fits_index<I>(d.cols(), "the number of columns");
        std::vector<I> row_offsets(d.rows() + 1, I(0));
        std::vector<I> column_indices;
        std::vector<T> values;
        for (index_t i = 1; i <= d.rows(); ++i)
        {
            for (index_t j = 1; j <= d.cols(); ++j)
                if (d(i, j) != T(0))
                {
                    column_indices.push_back(static_cast<I>(j - 1));
                    values.push_back(d(i, j));
                }
            detail::check_fits_index<I>(values.size(), "the number of entries");
            row_offsets[i] = static_cast<I>(values.size());
        }
        return CSRMatrix<T, I>(d.rows(), d.cols(), std::move(row_offsets), std::move(column_indices), std::move(values));
    }

    /*********************************** Scaling *************************************************/

    // Sparse matrix with every stored value multiplied by a scalar (the same pattern)
    template<class T, class I, Scalar S>
    CSRMatrix<detail::promote_t<T, S>, I> operator* (CSRMatrix<T, I> const& s, S scale)
    {
        typedef detail::promote_t<T, S> R;
        std::vector<R> values(s.values().begin(), s.values().end());
        for (R& v : values) v *= scale;
        return CSRMatrix<R, I>(s.rows(), s.cols(),
                               std::vector<I>(s.row_offsets().begin(), s.row_offsets().end()),
                               std::vector<I>(s.column_indices().begin(), s.column_indices().end()),
                               std::move(values));
    }

    template<class T, class I, Scalar S>
    CSRMatrix<detail::promote_t<T, S>, I> operator* (S scale, CSRMatrix<T, I> const& s)
    {return s * scale;}

    template<class T, class I, Scalar S>
    CSRMatrix<detail::promote_t<T, S>, I> operator/ (CSRMatrix<T, I> const& s, S scale)
    {
        CSRMatrix<detail::promote_t<T, S>, I> r = s * 1;
        r /= scale;
        return r;
    }

    template<class T, class I>
    CSRMatrix<T, I> operator- (CSRMatrix<T, I> const& s)
    {
        CSRMatrix<T, I> r = s;
        for (T& v : r.values()) v = -v;
        return r;
    }

    /*********************************** Output **************************************************/

    // The dimensions, then the stored entries, one "(i,j)  value" per line (1-based indices)
    template<class T, class I>
    std::ostream& operator<< (std::ostream& os, CSRMatrix<T, I> const& s)
    {
        os << s.rows() << "x" << s.cols() << " sparse matrix, " << s.nnz() << " stored entries";
        const auto offsets = s.row_offsets();
        const auto columns = s.column_indices();
        const auto values  = s.values();
        for (index_t i = 1; i <= s.rows(); ++i)
            for (I p = offsets[i - 1]; p < offsets[i]; ++p)
                os << "\n(" << i << "," << static_cast<index_t>(columns[p]) + 1 << ")  " << values[p];
        return os;
    }

    /************************** Dense factorizations of sparse matrices *************************/
    // A dense factorization of a sparse matrix would fill it in and cost O(n^3), so these
    // overloads turn the call into a compile error with a message. (They return a placeholder so
    // that the message is the only error, e.g. for 'auto F = lu(S);'.)

    #define MATL_SPARSE_REJECTION_MESSAGE \
        "a dense factorization of a sparse matrix fills it in and costs O(n^3): use an iterative " \
        "solver, or convert a small matrix explicitly with dense(S)"

    template<SparseMatrixType S>
    auto lu(S const&) {static_assert(detail::dependent_false_v<S>, "lu: " MATL_SPARSE_REJECTION_MESSAGE); return 0;}

    template<SparseMatrixType S, typename B>
    auto lu_solve(S const&, B const&) {static_assert(detail::dependent_false_v<S>, "lu_solve: " MATL_SPARSE_REJECTION_MESSAGE); return 0;}

    template<SparseMatrixType S>
    auto cholesky(S const&) {static_assert(detail::dependent_false_v<S>, "cholesky: " MATL_SPARSE_REJECTION_MESSAGE); return 0;}

    template<SparseMatrixType S, typename B>
    auto cholesky_solve(S const&, B const&) {static_assert(detail::dependent_false_v<S>, "cholesky_solve: " MATL_SPARSE_REJECTION_MESSAGE); return 0;}

    template<SparseMatrixType S>
    auto qr(S const&) {static_assert(detail::dependent_false_v<S>, "qr: " MATL_SPARSE_REJECTION_MESSAGE); return 0;}

    template<SparseMatrixType S, typename B>
    auto qr_solve(S const&, B const&) {static_assert(detail::dependent_false_v<S>, "qr_solve: " MATL_SPARSE_REJECTION_MESSAGE); return 0;}

    #undef MATL_SPARSE_REJECTION_MESSAGE
}
#endif
