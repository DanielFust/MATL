#ifndef _MATL_MATRIX_IMPLEMENTATION__HPP
#define _MATL_MATRIX_IMPLEMENTATION__HPP
namespace matl 
{

   /***********************************************************************************************
    * void operator= (T val)
    * 
    * Scalar data writing. Downcasts to derived type.
    * 
    * Parameters
    *      T val            Scalar value
    ***********************************************************************************************/
    /*
    template<typename MAT>
    template<class T, enable_if_scalar<T> >
    void _MATRIX_<MAT>::operator= (T val)
    {static_cast<MAT&>(*this)=val;}
    */


   /***********************************************************************************************
    * MAT1& eval(_MATRIX_& target) const
    * 
    * Generic matrix copy. Downcasts to derived type.
    * 
    * Parameters
    *     _MATRIX_ mat        Template matrix  
    ***********************************************************************************************/
    template<typename MAT1>
    template<typename MAT2>
    MAT1& _MATRIX_<MAT1>::operator= (_MATRIX_<MAT2>& mat)
    {
        static_cast<MAT1&>(*this)=mat;
        return *this;
    }



   /***********************************************************************************************
    *  MAT& operator= (_EXPRESSION_& expr)
    * 
    *  Trigger expression evaluation to this matrix. 
    * 
    * Parameters
    *     _EXPRESSION_ expr      Matrix expression  
    ***********************************************************************************************/
    template<typename MAT>
    template<typename EXPR>
    MAT& _MATRIX_<MAT>::operator=(_EXPRESSION_<EXPR> const& expr)
    {
        expr.eval(*this);
        return *this;
    }


   
   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Single index getter. Downcasts to derived type.
    * 
    * Parameters
    *     index_t i         index  
    * Returns:
    *     dtype            Element at index
    ***********************************************************************************************/
    template<typename MAT>
    typename _MATRIX_<MAT>::dtype _MATRIX_<MAT>::operator() (index_t i) const 
    {return static_cast<MAT const&>(*this)(i);} 


   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Double index getter. Downcasts to derived type.
    * 
    * Parameters
    *     index_t i         Row index  
    *     index_t j         Column index  
    * Returns:
    *     dtype            Element at row, col
    ***********************************************************************************************/ 
    template<typename MAT>
    typename _MATRIX_<MAT>::dtype _MATRIX_<MAT>::operator() (index_t i, index_t j) const 
    {return static_cast<MAT const&>(*this)(i,j);}



   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Single index accessor. Downcasts to derived type.
    * 
    * Parameters
    *     index_t i         index  
    * Returns:
    *     dtype            Element at index
    ***********************************************************************************************/
    template<typename MAT>
    typename _MATRIX_<MAT>::dtype& _MATRIX_<MAT>::operator() (index_t i)
    {return static_cast<MAT&>(*this)(i);}

    

   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Double index accessor. Downcasts to derived type.
    * 
    * Parameters
    *     index_t i         Row index  
    *     index_t j         Column index  
    * Returns:
    *     dtype&            Element at row, col
    ***********************************************************************************************/ 
    template<typename MAT>
    typename _MATRIX_<MAT>::dtype& _MATRIX_<MAT>::operator() (index_t i, index_t j)
    {return static_cast<MAT&>(*this)(i,j);}




   /***********************************************************************************************
    * index_t rows() const
    * 
    * Number of rows in matrix. Downcasts to derived type.
    * 
    *  Returns:
    *      index_t          Number of rows in matrix
    ***********************************************************************************************/
    template<typename MAT>
    index_t _MATRIX_<MAT>::rows() const
    {return static_cast<MAT const&>(*this).rows();}



   /***********************************************************************************************
    * index_t cols() const
    * 
    * Number of columns in matrix. Downcasts to derived type.
    * 
    *  Returns:
    *     index_t          Number of columns in matrix
    ***********************************************************************************************/
    template<typename MAT>
    index_t _MATRIX_<MAT>::cols() const
    {return static_cast<MAT const&>(*this).cols();}



   /***********************************************************************************************
    * index_t size() const
    * 
    * Number of elements in matrix. Downcasts to derived type.
    * 
    *  Returns:
    *      index_t          Number of elements in matrix
    ***********************************************************************************************/
    template<typename MAT>
    index_t _MATRIX_<MAT>::size() const
    {return static_cast<MAT const&>(*this).size();}

}
#endif