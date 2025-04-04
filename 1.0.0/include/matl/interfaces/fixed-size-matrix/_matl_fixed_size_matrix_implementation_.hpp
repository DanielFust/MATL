#ifndef _MATL_FIXED_SIZE_MATRIX_IMPLEMENTATION__HPP
#define _MATL_FIXED_SIZE_MATRIX_IMPLEMENTATION__HPP
namespace matl
{

   /***********************************************************************************************
    * FSM& operator= (_MATRIX_&)
    * 
    * Matrix copy.
    * 
    * Parameters:
    *     _MATRIX_ mat        Matrix to copy
    ***********************************************************************************************/
    template<typename FSM>
    template<typename MAT>
    FSM& _FIXED_SIZE_MATRIX_<FSM>::operator= (_MATRIX_<MAT>& mat)
    {
        static_cast<FSM&>(*this)=mat;
        return *this;
    }
   
   /***********************************************************************************************
    * FSM& operator= (_EXPRESSION_&)
    * 
    * Trigger expression evaluation to this matrix. 
    * 
    * Parameters:
    *     _EXPRESSION_ expr      Matrix expression  
    ***********************************************************************************************/
    template<typename FSM>
    template<typename EXPR>
    FSM& _FIXED_SIZE_MATRIX_<FSM>::operator=(_EXPRESSION_<EXPR> const& expr)
    {
        expr.eval(*this);
        return *this;
    }

   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Single index getter. Downcasts to derived type.
    * 
    * Parameters:
    *     index_t i         Index
    * Returns:
    *     dtype            Element at index
    ***********************************************************************************************/
    template<typename FSM>
    typename _FIXED_SIZE_MATRIX_<FSM>::dtype _FIXED_SIZE_MATRIX_<FSM>::operator() (index_t i) const
    {return static_cast<FSM const&>(*this)(i);}
    
   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Double index getter. Downcasts to derived type.
    * 
    * Parameters:
    *     index_t i         Row index  
    *     index_t j         Column index  
    * Returns:
    *     dtype            Element at row, col
    ***********************************************************************************************/
    template<typename FSM>
    typename _FIXED_SIZE_MATRIX_<FSM>::dtype _FIXED_SIZE_MATRIX_<FSM>::operator() (index_t i, index_t j) const
    {return static_cast<FSM const&>(*this)(i,j);}
    
   /***********************************************************************************************
    * dtype& operator() (index_t i)
    * 
    * Single index accessor. Downcasts to derived type.
    * 
    * Parameters:
    *     index_t i         Index  
    * Returns:
    *     dtype&           Element at index
    ***********************************************************************************************/
    template<typename FSM>
    typename _FIXED_SIZE_MATRIX_<FSM>::dtype& _FIXED_SIZE_MATRIX_<FSM>::operator() (index_t i)
    {return static_cast<FSM&>(*this)(i);}
    
   /***********************************************************************************************
    * dtype& operator() (index_t i, index_t j)
    * 
    * Double index accessor. Downcasts to derived type.
    * 
    * Parameters:
    *     index_t i         Row index  
    *     index_t j         Column index  
    * Returns:
    *     dtype&           Element at row, col
    ***********************************************************************************************/ 
    template<typename FSM>
    typename _FIXED_SIZE_MATRIX_<FSM>::dtype& _FIXED_SIZE_MATRIX_<FSM>::operator() (index_t i, index_t j)
    {return static_cast<FSM&>(*this)(i,j);}
    
   /***********************************************************************************************
    * index_t rows() const
    * 
    * Get number of rows of the matrix
    * 
    * Returns:
    *     index_t          Number of rows in matrix
    ***********************************************************************************************/  
    template<typename FSM>        
    index_t _FIXED_SIZE_MATRIX_<FSM>::rows() const
    {return FSM::ROWS;}
    
   /***********************************************************************************************
    * index_t cols() const
    * 
    * Get number of columns of the matrix
    * 
    * Returns:
    *     index_t          Number of columns in matrix
    ***********************************************************************************************/
    template<typename FSM>
    index_t _FIXED_SIZE_MATRIX_<FSM>::cols() const
    {return FSM::COLS;}
    
    
   /***********************************************************************************************
    * index_t size() const
    * 
    * Get the number of elements of the matrix
    * 
    * Returns:
    *    index_t          Number of elements in matrix
    ***********************************************************************************************/
    template<typename FSM>
    index_t _FIXED_SIZE_MATRIX_<FSM>::size() const
    {return FSM::SIZE;} 
}
#endif   