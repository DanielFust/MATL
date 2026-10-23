/**************************************************************************************************
 * File: _matl_sparse_builder_.hpp
 *
 * Description: SparseBuilder: assembles a sparse matrix from (row, column, value) entries.
 *
 *     SparseBuilder<double> B(n, n);
 *     for (each element e)                      // e.g. finite-element assembly
 *         for (each local entry (a, b))
 *             B.add(row(e, a), col(e, b), k_e(a, b));
 *     CSRMatrix<double> K = B.build();
 *
 * Entries may be added in any order, and the same (row, column) may be added more than once:
 * build() sums the values of duplicate entries (in the order they were added), which is what
 * assembly needs. Entries whose value is 0 are kept, so that a matrix assembled once keeps its
 * pattern even where values happen to cancel.
 *
 * Indices are 1-based and always checked (matl::index_error): an index outside the matrix would
 * otherwise corrupt the structure.
 **************************************************************************************************/
#ifndef MATL_SPARSE_BUILDER_HPP
#define MATL_SPARSE_BUILDER_HPP
namespace matl
{
   /***********************************************************************************************
    * class SparseBuilder<T, I>
    *
    * Collects entries for a rows x cols CSRMatrix<T, I> (see the description above).
    ***********************************************************************************************/
    template<class T, class I>
    class SparseBuilder
    {
      private:
        index_t        rows_;
        index_t        cols_;
        std::vector<I> row_;        // 0-based rows of the entries, in the order added
        std::vector<I> col_;        // 0-based columns
        std::vector<T> value_;

      public:
        // Builder for a rows x cols matrix
        SparseBuilder(index_t rows, index_t cols) : rows_(rows), cols_(cols)
        {
            detail::check_fits_index<I>(rows, "the number of rows");
            detail::check_fits_index<I>(cols, "the number of columns");
        }

        index_t rows() const {return this->rows_;}
        index_t cols() const {return this->cols_;}

        // Number of entries added (duplicates included)
        index_t entries() const {return this->value_.size();}

        // Reserves memory for the given number of entries
        void reserve(index_t entries)
        {
            this->row_.reserve(entries);
            this->col_.reserve(entries);
            this->value_.reserve(entries);
        }

        // Adds value at (i,j) (1-based); adding (i,j) again adds to its value
        void add(index_t i, index_t j, T value)
        {
            if ((i < 1) || (i > this->rows_) || (j < 1) || (j > this->cols_))
                detail::index_out_of_range(i, j, this->rows_, this->cols_);
            this->row_.push_back(static_cast<I>(i - 1));
            this->col_.push_back(static_cast<I>(j - 1));
            this->value_.push_back(value);
        }

        // Removes all entries (keeps the dimensions and the reserved memory)
        void clear()
        {
            this->row_.clear();
            this->col_.clear();
            this->value_.clear();
        }

        // The sparse matrix of the entries: sorted by row and column, duplicates summed. Takes
        // O(entries log(entries per row)) time; the builder is unchanged, so it can build again.
        CSRMatrix<T, I> build() const
        {
            const index_t m = this->rows_;
            const index_t count = this->value_.size();
            detail::check_fits_index<I>(count, "the number of entries");

            // Counting sort by row: start[i] is the first position of row i
            std::vector<index_t> start(m + 1, 0);
            for (index_t k = 0; k < count; ++k) ++start[static_cast<index_t>(this->row_[k]) + 1];
            for (index_t i = 0; i < m; ++i) start[i + 1] += start[i];

            std::vector<std::pair<I, T>> by_row(count);
            {
                std::vector<index_t> next(start.begin(), start.end() - 1);
                for (index_t k = 0; k < count; ++k)
                    by_row[next[static_cast<index_t>(this->row_[k])]++] = {this->col_[k], this->value_[k]};
            }

            // Within each row: sort by column (stable, so duplicates are summed in the order they
            // were added) and sum duplicates
            std::vector<I> row_offsets(m + 1, I(0));
            std::vector<I> column_indices;
            std::vector<T> values;
            column_indices.reserve(count);
            values.reserve(count);
            for (index_t i = 0; i < m; ++i)
            {
                const auto first = by_row.begin() + static_cast<std::ptrdiff_t>(start[i]);
                const auto end   = by_row.begin() + static_cast<std::ptrdiff_t>(start[i + 1]);
                if (end - first <= 32)
                {
                    // Insertion sort: rows are usually short (e.g. a few dozen entries in
                    // finite-element matrices), and std::stable_sort would allocate per row
                    for (auto it = first; it != end; ++it)
                    {
                        auto entry = *it;
                        auto hole = it;
                        for (; hole > first && (hole - 1)->first > entry.first; --hole) *hole = *(hole - 1);
                        *hole = entry;
                    }
                }
                else
                    std::stable_sort(first, end, [](auto const& a, auto const& b) {return a.first < b.first;});
                for (auto it = first; it != end; ++it)
                {
                    if (it != first && it->first == column_indices.back()) values.back() += it->second;
                    else {column_indices.push_back(it->first); values.push_back(it->second);}
                }
                row_offsets[i + 1] = static_cast<I>(values.size());
            }

            return CSRMatrix<T, I>(m, this->cols_, std::move(row_offsets), std::move(column_indices), std::move(values));
        }
    };
}
#endif
