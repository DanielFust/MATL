#ifndef MATL_EYE_SPECIALIZATIONS_HPP
#define MATL_EYE_SPECIALIZATIONS_HPP

namespace matl 
{
   /***********************************************************************************************
    * Matrix eye(index_t, index_t)
    *
    * Creates an identity matrix of size MxN.
    * 
    * Parameters:
    *     M           Number of rows
    *     N           Number of columns
    * Returns:
    *     Matrix      Identity matrix of size MxN
    ***********************************************************************************************/
   template<typename T>
   Matrix<T> eye(index_t M, index_t N)
   {
       Matrix<T> mat(M,N);
       for (index_t i = 1; i <= M; ++i)
           for (index_t j = 1; j <= N; ++j)
               mat(i,j) = (i==j);
       return mat;
   }
   
  /***********************************************************************************************
   * Matrix eye(index_t)
   *
   * Creates an identity matrix of size NxN.
   *
   * Parameters:
   *     N           Number of rows and columns
   * Returns:
   *     Matrix      Identity matrix of size NxN
   ***********************************************************************************************/
   template<typename T>
   Matrix<T> eye(index_t N)
   {return eye<T>(N,N);} 

  /***********************************************************************************************
   * Matrix eye<index_t,index_t>
   * 
   * Fixed size identity matrix
   *
   * (Template) Parameters:
   *     M           Number of rows
   *     N           Number of columns
   * Returns:
   *     Matrix      Fixed size identity matrix of size MxN
   ***********************************************************************************************/
   template<typename T, index_t M, index_t N>  
    Matrix<T,M,N> eye()
    {
        Matrix<T,M,N> mat;
        for (index_t i = 1; i <= M; ++i)
            for (index_t j = 1; j <= N; ++j)
                mat(i,j) = (i==j);
        return mat;
    }

   /***********************************************************************************************
    * Matrix eye<index_t,index_t>
    * 
    * Fixed size identity matrix
    *
    * (Template) Parameters:
    *     N           Number of rows and columns
    * Returns:
    *     Matrix      Fixed size identity matrix of size NxN
    ***********************************************************************************************/
    template<typename T, index_t N>
    Matrix<T,N,N> eye()
    {
        Matrix<T,N,N> mat;
        for (index_t i = 1; i <= N; ++i)
            for (index_t j = 1; j <= N; ++j)
                mat(i,j) = (i==j);
                
        return mat;
    }
}

#endif