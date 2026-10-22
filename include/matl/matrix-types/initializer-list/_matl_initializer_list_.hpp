/**************************************************************************************************
 * File: _matl_initializer_list_.hpp
 *
 * Description: Helpers for constructing and assigning matrices from nested initializer lists,
 *              one inner list per row:
 *
 *                  Matrix<double> A = {{1, 2, 3},
 *                                      {4, 5, 6}};      // 2x3
 **************************************************************************************************/
#ifndef MATL_INITIALIZER_LIST_HPP
#define MATL_INITIALIZER_LIST_HPP

namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * void initializer_list_dimensions(rows, M, N)
        *
        * Dimensions of a nested initializer list: M rows of N elements. Reports a
        * matl::dimension_error if the rows have different lengths (always checked: an initializer
        * list is written by hand, and the check is not in an inner loop).
        *******************************************************************************************/
        template<class T>
        void initializer_list_dimensions(matrix_rows<T> rows, index_t& M, index_t& N)
        {
            M = rows.size();
            N = (M > 0) ? rows.begin()->size() : 0;
            index_t i = 1;
            for (std::initializer_list<T> const& row : rows)
            {
                if (row.size() != N)
                    raise<dimension_error>("the rows of a matrix initializer list must have the same length: row " +
                                           std::to_string(i) + " has " + std::to_string(row.size()) +
                                           " elements, row 1 has " + std::to_string(N));
                ++i;
            }
        }

       /*******************************************************************************************
        * void copy_initializer_list(target, rows)
        *
        * Copies the elements of a nested initializer list into a matrix of its dimensions.
        *******************************************************************************************/
        template<MatrixType MAT, class T>
        void copy_initializer_list(MAT& target, matrix_rows<T> rows)
        {
            index_t i = 1;
            for (std::initializer_list<T> const& row : rows)
            {
                index_t j = 1;
                for (T const& value : row) target(i, j++) = value;
                ++i;
            }
        }
    }
}
#endif
