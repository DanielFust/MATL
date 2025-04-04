#ifndef _MATL_API_FUNCTIONS_DECLARATIONS__HPP
#define _MATL_API_FUNCTIONS_DECLARATIONS__HPP

namespace matl
{
    // eye (make identity matrix)
    template<class T>                       Matrix<T>     eye(index_t);
    template<class T>                       Matrix<T>     eye(index_t,index_t);
    template<class T, index_t N>            Matrix<T,N,N> eye();
    template<class T, index_t N, index_t M> Matrix<T,M,N> eye();
} 
#endif
