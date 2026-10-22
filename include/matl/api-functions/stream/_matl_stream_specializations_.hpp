/**************************************************************************************************
 * File: _matl_stream_specializations_.hpp
 *
 * Description: Writing matrices and expressions to output streams.
 **************************************************************************************************/
#ifndef MATL_STREAM_SPECIALIZATIONS_HPP
#define MATL_STREAM_SPECIALIZATIONS_HPP

namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * void write_matrix(os, mat)
        *
        * Writes the rows of a matrix on separate lines, with the columns right-aligned and
        * separated by two spaces. Elements are formatted with the stream's settings (e.g.
        * std::setprecision, std::fixed). No newline follows the last row.
        *******************************************************************************************/
        template<MatrixType MAT>
        void write_matrix(std::ostream& os, MAT const& mat)
        {
            const index_t M = mat.rows();
            const index_t N = mat.cols();

            // Format every element with the stream's settings, and find the width of each column
            std::ostringstream field;
            field.copyfmt(os);
            field.width(0);
            std::vector<std::string> text(M * N);
            std::vector<std::size_t> width(N, 0);
            for (index_t i = 1; i <= M; ++i)
                for (index_t j = 1; j <= N; ++j)
                {
                    field.str("");
                    field << mat(i, j);
                    text[(i - 1) * N + (j - 1)] = field.str();
                    width[j - 1] = std::max(width[j - 1], field.str().size());
                }

            os.width(0);
            for (index_t i = 1; i <= M; ++i)
            {
                if (i > 1) os << '\n';
                for (index_t j = 1; j <= N; ++j)
                {
                    std::string const& t = text[(i - 1) * N + (j - 1)];
                    if (j > 1) os << "  ";
                    os << std::string(width[j - 1] - t.size(), ' ') << t;
                }
            }
        }
    }

   /***********************************************************************************************
    * std::ostream& operator<< (std::ostream&, MAO const&)
    *
    * Writes a matrix, or the value of an expression, to a stream (see detail::write_matrix):
    *
    *     std::cout << A << '\n';
    *     std::cout << std::setprecision(3) << A * B << '\n';
    *
    * An expression is evaluated into a temporary matrix first.
    ***********************************************************************************************/
    template<MatrixAlgebraObject MAO>
    std::ostream& operator<< (std::ostream& os, MAO const& mao)
    {
        if constexpr (MatrixType<MAO>)
            detail::write_matrix(os, mao);
        else
        {
            const detail::result_matrix_t<MAO> value(mao);
            detail::write_matrix(os, value);
        }
        return os;
    }
}
#endif
