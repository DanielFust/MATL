#ifndef _MATL_MATRIX_ALGEBRA_OBJECT_IMPLEMENTATION__HPP
#define _MATL_MATRIX_ALGEBRA_OBJECT_IMPLEMENTATION__HPP

namespace matl
{
   /***********************************************************************************************
    * index_t rows() const
    * 
    * Returns:
    *     Number of rows in the matrix algebra object
    ***********************************************************************************************/
    template<typename MAO>
    index_t _MATRIX_ALGEBRA_OBJECT_<MAO>::rows() const
    {return static_cast<MAO const&>(*this).rows();}

   /***********************************************************************************************
    * index_t cols() const
    * 
    * Returns:
    *     Number of columns in the matrix algebra object
    ***********************************************************************************************/
    template<typename MAO>
    index_t _MATRIX_ALGEBRA_OBJECT_<MAO>::cols() const
    {return static_cast<MAO const&>(*this).cols();}  

   /***********************************************************************************************
    * index_t size() const
    * 
    * Returns:
    *     Number of elements in the matrix algebra object
    ***********************************************************************************************/
    template<typename MAO>
    index_t _MATRIX_ALGEBRA_OBJECT_<MAO>::size() const
    {return static_cast<MAO const&>(*this).size();}  
    
}

#endif