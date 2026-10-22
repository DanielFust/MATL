#ifndef MATL_API_FUNCTIONS_DECLARATIONS_HPP
#define MATL_API_FUNCTIONS_DECLARATIONS_HPP

namespace matl
{
    // eye (make identity matrix)
    template<class T>                       Matrix<T>     eye(index_t);
    template<class T>                       Matrix<T>     eye(index_t,index_t);
    template<class T, index_t N>            Matrix<T,N,N> eye();
    template<class T, index_t M, index_t N> Matrix<T,M,N> eye();

    // zeros (every element zero)
    template<class T>                       Matrix<T>     zeros(index_t);
    template<class T>                       Matrix<T>     zeros(index_t,index_t);
    template<class T, index_t N>            Matrix<T,N,N> zeros();
    template<class T, index_t M, index_t N> Matrix<T,M,N> zeros();

    // ones (every element one)
    template<class T>                       Matrix<T>     ones(index_t);
    template<class T>                       Matrix<T>     ones(index_t,index_t);
    template<class T, index_t N>            Matrix<T,N,N> ones();
    template<class T, index_t M, index_t N> Matrix<T,M,N> ones();

    // Write a matrix or the value of an expression to a stream
    template<MatrixAlgebraObject MAO> std::ostream& operator<< (std::ostream&, MAO const&);
} 
#endif
