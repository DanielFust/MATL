//#ifndef _MATL_MATRIX_IMPLEMENTATION__HPP
//#define _MATL_MATRIX_IMPLEMENTATION__HPP

#include <iostream>
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"


namespace matl 
{

   /***********************************************************************************************
    * void operator= (T val)
    * 
    * Scalar data writing. Downcasts to derived type.
    * 
    * Template Parameters
    *      MAT_T            Matrix type 
    *      T                Scalar type
    * 
    * Parameters
    *      T val            Scalar value
    ***********************************************************************************************/
    template<typename MAT_T>
    template<class T, enable_if_scalar<T> >
    void Matrix<MAT_T>::operator= (T val)
    {static_cast<MAT_T&>(*this)=val;}



   /***********************************************************************************************
    * void operator= (_MATEXPRESSION_& expr)
    * 
    * Generic matrix copy. Downcasts to derived type.
    * 
    *  Template Parameters
    *      MAT_T1           Matrix type (this)
    *      MAT_T2           Matrix type (template)
    * 
    * Parameters
    *     Matrix mat        Template matrix  
    ***********************************************************************************************/
    template<typename MAT_T1>
    template<typename MAT_T2>
    void Matrix<MAT_T1>::operator= (Matrix<MAT_T2>& mat)
    {static_cast<MAT_T1&>(*this)=mat;}



   /***********************************************************************************************
    * void eval(Matrix& target) const
    * 
    *  Trigger expression evaluation to this matrix. 
    * 
    *  Template Parameters
    *      MAT_T                    Matrix type (this)
    *      EXPR_T                   Expression type
    * 
    * Parameters
    *     _MATEXPRESSION_ expr      Matrix expression  
    ***********************************************************************************************/
    template<typename MAT_T>
    template<typename EXPR_T>
    void Matrix<MAT_T>::operator= (_MATEXPRESSION_<EXPR_T>& expr)
    {expr.eval(*this);}



   /***********************************************************************************************
    * Matrix blank_copy() const
    * 
    *  Create a blank Matrix of same type, data type, and size.
    * 
    * Template Parameters:
    *     MAT_T                    Matrix type
    * 
    * Parameters:
    * Returns:
    *     Matrix                   Blank copy of matrix (this)            
    * 
    ***********************************************************************************************/
    template<typename MAT_T>
    Matrix<MAT_T> Matrix<MAT_T>::blank_copy() const
    {return static_cast<MAT_T const&>(*this).blank_copy();}  


   
   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Single index getter. Downcasts to derived type.
    * 
    *  Template Parameters
    *      MAT_T            Matrix type
    * 
    * Parameters
    *     index_t i         index  
    ***********************************************************************************************/
    template<typename MAT_T>
    typename Matrix<MAT_T>::dtype Matrix<MAT_T>::operator() (index_t i) const 
    {return static_cast<MAT_T const&>(*this)(i);} 


   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Double index getter. Downcasts to derived type.
    * 
    *  Template Parameters
    *     MAT_T            Matrix type
    * 
    * Parameters
    *     index_t i         Row index  
    *     index_t j         Column index  
    ***********************************************************************************************/ 
    template<typename MAT_T>
    typename Matrix<MAT_T>::dtype Matrix<MAT_T>::operator() (index_t i, index_t j) const 
    {return static_cast<MAT_T const&>(*this)(i,j);}



   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Single index accessor. Downcasts to derived type.
    * 
    *  Template Parameters
    *      MAT_T            Matrix type
    * 
    * Parameters
    *     index_t i         index  
    ***********************************************************************************************/
    template<typename MAT_T>
    typename Matrix<MAT_T>::dtype& Matrix<MAT_T>::operator() (index_t i)
    {return static_cast<MAT_T&>(*this)(i);}

    

   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Double index accessor. Downcasts to derived type.
    * 
    *  Template Parameters
    *     MAT_T            Matrix type
    * 
    * Parameters
    *     index_t i         Row index  
    *     index_t j         Column index  
    ***********************************************************************************************/ 
    template<typename MAT_T>
    typename Matrix<MAT_T>::dtype& Matrix<MAT_T>::operator() (index_t i, index_t j)
    {return static_cast<MAT_T&>(*this)(i,j);}




   /***********************************************************************************************
    * index_t rows() const
    * 
    * Number of rows in matrix. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      MAT_T            Matrix type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename MAT_T>
    index_t Matrix<MAT_T>::rows() const
    {return static_cast<MAT_T const&>(*this).rows();}



   /***********************************************************************************************
    * index_t cols() const
    * 
    * Number of columns in matrix. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      MAT_T            Matrix type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename MAT_T>
    index_t Matrix<MAT_T>::cols() const
    {return static_cast<MAT_T const&>(*this).cols();}



   /***********************************************************************************************
    * index_t size() const
    * 
    * Number of elements in matrix. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      MAT_T            Matrix type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename MAT_T>
    index_t Matrix<MAT_T>::size() const
    {return static_cast<MAT_T const&>(*this).size();}

    
   /***********************************************************************************************
    * index_t size() const
    * 
    * Number of elements in buffer. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      MAT_T            Matrix type
    * 
    * Parameters:
    ***********************************************************************************************/ 
    template<typename MAT_T>
    index_t Matrix<MAT_T>::buffer() const
    {return static_cast<MAT_T const&>(*this).buffer();} 




   /***********************************************************************************************
    * void resize(index_t, index_t) 
    * 
    * Resize matrix to the specified shape. Attempts to use buffer if size increases
    * 
    * Template Parameters:
    *     MAT_T            Matrix type
    * 
    * Parameters:
    ***********************************************************************************************/ 
    template<typename MAT_T>
    void Matrix<MAT_T>::resize(index_t rows, index_t cols)
    {static_cast<MAT_T&>(*this).resize(rows, cols);}


    template<typename MAT_T>
    void Matrix<MAT_T>::println() const
    {static_cast<MAT_T const&>(*this).println();}

}


//#endif