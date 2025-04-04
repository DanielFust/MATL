#ifndef _MATL_TRANSPOSE__HPP
#define _MATL_TRANSPOSE__HPP

namespace matl
{
    template<MatrixAlgebraObject MAO>
    _MATRIX_TRANSPOSE_<MAO> transpose(MAO const& mat)
    {return _MATRIX_TRANSPOSE_<MAO>(mat);}
}

#endif