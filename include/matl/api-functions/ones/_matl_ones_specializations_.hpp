#ifndef MATL_ONES_SPECIALIZATIONS_HPP
#define MATL_ONES_SPECIALIZATIONS_HPP

namespace matl
{
   /***********************************************************************************************
    * Matrix ones<T>(M, N), ones<T>(N), ones<T, M, N>(), ones<T, N>()
    *
    * Matrices with every element one: dynamic (M x N, or N x N) or fixed-size.
    ***********************************************************************************************/
    template<class T>
    Matrix<T> ones(index_t M, index_t N)
    {
        Matrix<T> mat(M, N);
        mat = T(1);
        return mat;
    }

    template<class T>
    Matrix<T> ones(index_t N)
    {return ones<T>(N, N);}

    template<class T, index_t M, index_t N>
    Matrix<T,M,N> ones()
    {
        Matrix<T,M,N> mat;
        mat = T(1);
        return mat;
    }

    template<class T, index_t N>
    Matrix<T,N,N> ones()
    {return ones<T,N,N>();}
}
#endif
