#ifndef _MATL_ELEMENTWISE_IMPLEMENTATION__HPP
#define _MATL_ELEMENTWISE_IMPLEMENTATION__HPP
namespace matl
{
   /***********************************************************************************************
    * operator() (index_t, index_t) const
    * 
    * Get indexed result of matrix or elementwise-expression
    * 
    * Parameters:
    *     row        Row index
    *     col        Column index
    * Returns:
    *     dtype      Element at row, col
    ***********************************************************************************************/
    template<typename EW>
    typename _ELEMENTWISE_<EW>::dtype _ELEMENTWISE_<EW>::operator() (index_t row, index_t col) const 
    {return static_cast<EW const&>(*this)(row, col);}

   /***********************************************************************************************
    * operator() (index_t) const
    * 
    * Get single-indexed result of matrix or elementwise-expression
    * 
    * Parameters:
    *     index_t i        Index of the element to get
    ***********************************************************************************************/
    template<typename EW>
    typename _ELEMENTWISE_<EW>::dtype _ELEMENTWISE_<EW>::operator() (index_t i) const 
    {return static_cast<EW const&>(*this)(i);} 
    
   /***********************************************************************************************
    * rows() const
    * 
    * Get number of rows of the matrix or expression
    * 
    * Returns:
    *     index_t          Number of rows in resultant matrix
    ***********************************************************************************************/  
    template<typename EW>        
    index_t _ELEMENTWISE_<EW>::rows() const 
    {return static_cast<EW const&>(*this).rows();}
    
   /***********************************************************************************************
    * cols() const
    * 
    * Get number of columns of matrix or expression
    * 
    * Returns:
    *     index_t  Number of columns in resultant matrix
    ***********************************************************************************************/
    template<typename EW>
    index_t _ELEMENTWISE_<EW>::cols() const 
    {return static_cast<EW const&>(*this).cols();}
    
    /***********************************************************************************************
    * size() const
    * 
    * Get the number of elements of the matrix or expression
    * 
    * Returns:
    *    index_t          Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename EW>
    index_t _ELEMENTWISE_<EW>::size() const
    {return static_cast<EW const&>(*this).size();} 
}
#endif