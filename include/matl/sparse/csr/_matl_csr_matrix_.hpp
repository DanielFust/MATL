/**************************************************************************************************
 * File: _matl_csr_matrix_.hpp
 *
 * Description: Sparse matrices in compressed sparse row (CSR) format.
 *
 * A CSR matrix stores only its nonzero entries (the "stored entries"), row by row, in three
 * arrays:
 *
 *     row_offsets      rows + 1 positions: the entries of row i are at positions
 *                      row_offsets[i-1] .. row_offsets[i]-1 of the two arrays below
 *     column_indices   the column of each entry, increasing within each row
 *     values           the value of each entry
 *
 * For example, the 3x4 matrix
 *
 *     [ 1  0  2  0 ]
 *     [ 0  0  0  0 ]          row_offsets    = {0, 2, 2, 4}
 *     [ 0  3  0  4 ]          column_indices = {0, 2, 1, 3}
 *                             values         = {1, 2, 3, 4}
 *
 * The arrays are 0-based, as in SciPy, MKL, PETSc and SuiteSparse, so that they can be exchanged
 * with those libraries without conversion. Element access is 1-based, like the rest of MATL:
 * S(1,3) is 2 in the example.
 *
 * Sparse matrices are built with a SparseBuilder (from (row, column, value) entries, e.g. during
 * finite-element assembly), converted from a dense matrix with sparse(D), or constructed from the
 * three arrays. The sparsity pattern (which entries are stored) is fixed once a matrix is built;
 * the values can be changed (values(), entry(i,j)). A new pattern means building a new matrix.
 *
 * Sparse matrices are not dense matrices or expressions (MatrixType, ExpressionType): dense
 * algorithms do not accept them, so that a sparse matrix is never silently treated as a dense
 * one. Products with dense matrices and vectors (S * x, transpose(S) * x) are expressions that
 * evaluate into dense matrices (see _matl_sparse_product_.hpp).
 *
 * Template Parameters:
 *     T         Element type (an arithmetic type)
 *     I         Index type of the column indices and row offsets (default index_t; a 32-bit type
 *               such as std::int32_t uses less memory and makes products faster, for matrices
 *               with fewer than 2^31 stored entries)
 **************************************************************************************************/
#ifndef MATL_CSR_MATRIX_HPP
#define MATL_CSR_MATRIX_HPP
namespace matl
{
    template<class T, class I = index_t> class CSRMatrix;
    template<class T, class I = index_t> class SparseBuilder;

    namespace expr
    {
        // transpose(S) of a sparse matrix (see _matl_sparse_product_.hpp)
        template<class T, class I> class SparseTransposed;
    }

    namespace detail
    {
        // Sparse matrix types
        template<typename S> inline constexpr bool is_sparse_v = false;
        template<class T, class I> inline constexpr bool is_sparse_v< CSRMatrix<T, I> > = true;

        // 'false' for any type (for static_assert in templates that must not be instantiated)
        template<typename> inline constexpr bool dependent_false_v = false;

        // Reports an index that does not fit the index type of a sparse matrix
        template<class I>
        void check_fits_index(index_t value, char const* what)
        {
            if (value > static_cast<index_t>(std::numeric_limits<I>::max()))
                detail::raise<dimension_error>(std::string("sparse matrix: ") + what + " (" + std::to_string(value) +
                                       ") does not fit the index type (at most " +
                                       std::to_string(std::numeric_limits<I>::max()) + ")");
        }
    }

    // SparseMatrixType: a sparse matrix (CSRMatrix), ignoring references and cv-qualifiers
    template<typename S>
    concept SparseMatrixType = detail::is_sparse_v< std::remove_cvref_t<S> >;


   /***********************************************************************************************
    * class CSRMatrix<T, I>
    *
    * Sparse matrix in compressed sparse row format (see the description above).
    ***********************************************************************************************/
    template<class T, class I>
    class CSRMatrix
    {
        static_assert(std::is_arithmetic_v<T>, "CSRMatrix: the element type must be an arithmetic type");
        static_assert(std::is_integral_v<I> && !std::is_same_v<I, bool>, "CSRMatrix: the index type must be an integer type");

      public:
        typedef T dtype;
        typedef I index_type;

      private:
        index_t        rows_;
        index_t        cols_;
        std::vector<I> row_offsets_;       // rows_ + 1 positions (0-based)
        std::vector<I> column_indices_;    // 0-based, increasing within each row
        std::vector<T> values_;

        // Position of the stored entry (i,j) (1-based indices) in values_, or nnz() if (i,j) is
        // not stored (a binary search within row i)
        index_t position(index_t i, index_t j) const
        {
            const auto first = this->column_indices_.begin() + static_cast<std::ptrdiff_t>(this->row_offsets_[i - 1]);
            const auto end   = this->column_indices_.begin() + static_cast<std::ptrdiff_t>(this->row_offsets_[i]);
            const auto it    = std::lower_bound(first, end, static_cast<I>(j - 1));
            return (it != end && *it == static_cast<I>(j - 1)) ? static_cast<index_t>(it - this->column_indices_.begin())
                                                               : this->nnz();
        }

        void check_index(index_t i, index_t j) const
        {
            if ((i < 1) || (i > this->rows_) || (j < 1) || (j > this->cols_))
                detail::index_out_of_range(i, j, this->rows_, this->cols_);
        }

        // Checks the three arrays (always: an invalid structure would make every operation read
        // out of bounds)
        void validate() const
        {
            const index_t nnz = this->values_.size();
            if (this->row_offsets_.size() != this->rows_ + 1)
                detail::raise<dimension_error>("CSRMatrix: row_offsets has " + std::to_string(this->row_offsets_.size()) +
                                       " elements, but a matrix with " + std::to_string(this->rows_) +
                                       " rows needs " + std::to_string(this->rows_ + 1));
            if (this->column_indices_.size() != nnz)
                detail::raise<dimension_error>("CSRMatrix: column_indices has " + std::to_string(this->column_indices_.size()) +
                                       " elements, values " + std::to_string(nnz));
            if (this->row_offsets_[0] != 0 || static_cast<index_t>(this->row_offsets_[this->rows_]) != nnz)
                detail::raise<index_error>("CSRMatrix: row_offsets must start at 0 and end at the number of entries (" +
                                   std::to_string(nnz) + ")");
            for (index_t i = 1; i <= this->rows_; ++i)
            {
                if (this->row_offsets_[i] < this->row_offsets_[i - 1])
                    detail::raise<index_error>("CSRMatrix: row_offsets decrease at row " + std::to_string(i));
                for (I p = this->row_offsets_[i - 1]; p < this->row_offsets_[i]; ++p)
                {
                    const I c = this->column_indices_[p];
                    bool negative = false;
                    if constexpr (std::is_signed_v<I>) negative = (c < 0);
                    if (negative || static_cast<index_t>(c) >= this->cols_)
                        detail::raise<index_error>("CSRMatrix: column index " + std::to_string(c) + " in row " + std::to_string(i) +
                                           " is out of range for " + std::to_string(this->cols_) +
                                           " columns (the arrays are 0-based)");
                    if (p > this->row_offsets_[i - 1] && c <= this->column_indices_[p - 1])
                        detail::raise<index_error>("CSRMatrix: the column indices of row " + std::to_string(i) +
                                           " must be strictly increasing");
                }
            }
        }

      public:

        /********************* Constructors *********************/

        // Empty 0 x 0 matrix
        CSRMatrix() : rows_(0), cols_(0), row_offsets_(1, I(0)), column_indices_(), values_() {}

        // rows x cols matrix with no stored entries (all zero)
        CSRMatrix(index_t rows, index_t cols) : rows_(rows), cols_(cols), row_offsets_(rows + 1, I(0)),
                                                column_indices_(), values_()
        {}

        // From the three CSR arrays (0-based, see the description above), which are checked and
        // moved into the matrix
        CSRMatrix(index_t rows, index_t cols, std::vector<I> row_offsets, std::vector<I> column_indices,
                  std::vector<T> values) :
            rows_(rows), cols_(cols), row_offsets_(std::move(row_offsets)),
            column_indices_(std::move(column_indices)), values_(std::move(values))
        {
            detail::check_fits_index<I>(cols, "the number of columns");
            this->validate();
        }

        // The transpose, formed explicitly: CSRMatrix<double> St = transpose(S);
        // (defined in _matl_sparse_arithmetic_.hpp)
        CSRMatrix(expr::SparseTransposed<T, I> const& transposed);

        /********************** Dimensions **********************/

        index_t rows() const {return this->rows_;}
        index_t cols() const {return this->cols_;}

        // Number of stored entries
        index_t nnz() const {return this->values_.size();}

        /*********************** Elements ***********************/

        // Element (i,j) (1-based): the stored value, or 0 if (i,j) is not stored
        T operator() (index_t i, index_t j) const
        {
            #ifdef MATL_BOUNDS_CHECK
            this->check_index(i, j);
            #endif
            const index_t p = this->position(i, j);
            return (p < this->nnz()) ? this->values_[p] : T(0);
        }

        // 'true' if the entry (i,j) is stored
        bool is_stored(index_t i, index_t j) const
        {
            this->check_index(i, j);
            return this->position(i, j) < this->nnz();
        }

        // The stored entry (i,j), to read or change its value; matl::index_error if (i,j) is not
        // stored (the pattern cannot change)
        T& entry(index_t i, index_t j)
        {
            this->check_index(i, j);
            const index_t p = this->position(i, j);
            if (p == this->nnz())
                detail::raise<index_error>("CSRMatrix: entry (" + std::to_string(i) + "," + std::to_string(j) +
                                   ") is not stored (the sparsity pattern is fixed)");
            return this->values_[p];
        }

        /************************* Arrays ***********************/

        // The CSR arrays (0-based). The structure is read-only; the values can be changed.
        std::span<const I> row_offsets()    const {return this->row_offsets_;}
        std::span<const I> column_indices() const {return this->column_indices_;}
        std::span<const T> values()         const {return this->values_;}
        std::span<T>       values()               {return this->values_;}

        /************************ Scaling ***********************/

        // Multiplies (divides) every stored value by a scalar
        template<Scalar S>
        CSRMatrix& operator*= (S s)
        {
            for (T& v : this->values_) v *= s;
            return *this;
        }

        template<Scalar S>
        CSRMatrix& operator/= (S s)
        {
            for (T& v : this->values_) v /= s;
            return *this;
        }
    };
}
#endif
