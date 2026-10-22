#ifndef MATL_ZEROS_SPECIALIZATIONS_HPP
#define MATL_ZEROS_SPECIALIZATIONS_HPP

namespace matl
{
   /***********************************************************************************************
    * Matrix zeros<T>(M, N), zeros<T>(N), zeros<T, M, N>(), zeros<T, N>()
    *
    * Matrices with every element zero: dynamic (M x N, or N x N) or fixed-size.
    ***********************************************************************************************/
    template<class T>
    Matrix<T> zeros(index_t M, index_t N)
    {
        Matrix<T> mat(M, N);
        mat = T(0);
        return mat;
    }

    template<class T>
    Matrix<T> zeros(index_t N)
    {return zeros<T>(N, N);}

    template<class T, index_t M, index_t N>
    Matrix<T,M,N> zeros()
    {
        Matrix<T,M,N> mat;
        mat = T(0);
        return mat;
    }

    template<class T, index_t N>
    Matrix<T,N,N> zeros()
    {return zeros<T,N,N>();}
}
#endif
