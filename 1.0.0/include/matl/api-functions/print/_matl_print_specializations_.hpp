#ifndef _MATL_PRINT_SPECIALIZATIONS__HPP
#define _MATL_PRINT_SPECIALIZATIONS__HPP

namespace matl 
{
   /***********************************************************************************************
    * void print(_MATRIX_) 
    * 
    * Prints the matrix to the standard output stream.
    *
    * Parameters:
    *     mat         Matrix to print
    ***********************************************************************************************/
   template<typename MAT>
   void print(_MATRIX_<MAT> const& mat)
   {
       const index_t M = mat.rows();
       const index_t N = mat.cols();
       for (index_t i = 1; i <= M; ++i)
       {
           for (index_t j = 1; j <= N; ++j)
               std::cout << mat(i,j) << " ";
           std::cout << std::endl;
       }
   } 

}

#endif