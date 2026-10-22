#ifndef MATL_API_FUNCTIONS_DECLARATIONS_HPP
#define MATL_API_FUNCTIONS_DECLARATIONS_HPP

namespace matl
{
    // eye (make identity matrix)
    template<class T>                       Matrix<T>     eye(index_t);
    template<class T>                       Matrix<T>     eye(index_t,index_t);
    template<class T, index_t N>            Matrix<T,N,N> eye();
    template<class T, index_t M, index_t N> Matrix<T,M,N> eye();
} 
#endif
