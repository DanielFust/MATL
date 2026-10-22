#ifndef MATL_TRANSPOSE_HPP
#define MATL_TRANSPOSE_HPP

namespace matl
{
    // Transpose of a matrix or expression
    template<MatrixAlgebraObject MAO>
    expr::Transpose<MAO> transpose(MAO const& mat)
    {return expr::Transpose<MAO>(mat);}
}
#endif
