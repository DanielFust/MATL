#ifndef MATL_PRINT_SPECIALIZATIONS_HPP
#define MATL_PRINT_SPECIALIZATIONS_HPP

namespace matl 
{
   /***********************************************************************************************
    * void print(MAT const&)
    * 
    * Prints the matrix to the standard output stream.
    *
    * Parameters:
    *     mat         Matrix to print
    ***********************************************************************************************/
   template<MatrixType MAT>
   void print(MAT const& mat)
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