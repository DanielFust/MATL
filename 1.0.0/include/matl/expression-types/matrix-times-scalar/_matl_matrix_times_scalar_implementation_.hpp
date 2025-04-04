/**************************************************************************************************
 * File: _matl_matrix_times_scalar_implementation_.hpp
 * 
 * Description: Implementations for the matrix times scalar expression class.
 **************************************************************************************************/
#ifndef _MATL_MATRIX_TIMES_SCALAR_IMPLEMENTATION__HPP
#define _MATL_MATRIX_TIMES_SCALAR_IMPLEMENTATION__HPP

namespace matl
{
   /***********************************************************************************************
    * Constructor (elementwise * scalar) 
    * 
    * Initialize matrix times scalar expression from a matrix or elementwise expression and a 
    * scalar.
    * 
    * Parameters:
    *     mat          Matrix or elementwise expression
    *    scalar        Scalar
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAO_, 
             enable_if_elementwise<MAO_>,
             enable_if_scalar<T> >
    _MATRIX_TIMES_SCALAR_<MAO,T>::_MATRIX_TIMES_SCALAR_(MAO_ const& mat, T scalar) : 
        mat_(unqualified_cast(mat)), 
        scalar_(scalar)
    {}

   /***********************************************************************************************
    * Constructor (non-elementwise * scalar)
    * 
    * Initialize matrix times scalar expression from a non-elementwise expression and a scalar.
    * 
    * Parameters:
    *     expr         Non-elementwise expression
    *     scalar       Scalar
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAO_, 
             disable_if_elementwise<MAO_>,
             enable_if_scalar<T> >
    _MATRIX_TIMES_SCALAR_<MAO,T>::_MATRIX_TIMES_SCALAR_(MAO_ const& expr, T scalar) : scalar_(scalar)
    {
        MemoryManager::allocate(this->mat_, expr.rows(), expr.cols());
        expr.eval(this->mat_);
    }
    
    
   /***********************************************************************************************
    * index_t rows() const
    * 
    * Returns:
    *      index_t      Number of rows in resultant matrix
    **********************************************************************************************/ 
    template<typename MAO, class T>
    index_t _MATRIX_TIMES_SCALAR_<MAO,T>::rows() const
    {return this->mat_.rows();}
    
   /***********************************************************************************************
    * index_t cols() const
    *
    * Returns: 
    *     index_t      Number of columns in resultant matrix
    ***********************************************************************************************/
    template<typename MAO, class T>
    index_t _MATRIX_TIMES_SCALAR_<MAO,T>::cols() const
    {return this->mat_.cols();}
    
   /***********************************************************************************************
    * index_t size() const
    *
    * Returns: 
    *     index_t      Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename MAO, class T>
    index_t _MATRIX_TIMES_SCALAR_<MAO,T>::size() const
    {return this->mat_.size();}
    
    
   /***********************************************************************************************
    * operator() (index_t i, index_t j) const
    * 
    * Parameters:
    *     i            Row index
    *     j            Column index
    * Returns:
    *     dtype        Value of resultant matrix at index pair
    ***********************************************************************************************/
    template<typename MAO, class T>
    typename _MATRIX_TIMES_SCALAR_<MAO,T>::dtype
    _MATRIX_TIMES_SCALAR_<MAO,T>::operator() (index_t i, index_t j) const
    {return this->mat_(i,j) * this->scalar_;}

   /***********************************************************************************************
    * operator() (index_t i) const
    * 
    * Parameters:
    *     i            Index of the element to get
    * Returns:
    *     dtype        Value of resultant matrix at index 
    *************************************************************************************************/
    template<typename MAO, class T>
    typename _MATRIX_TIMES_SCALAR_<MAO,T>::dtype
    _MATRIX_TIMES_SCALAR_<MAO,T>::operator() (index_t i) const
    {return this->mat_(i) * this->scalar_;} 

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const*)
    * 
    * Check if evaluation of the expression to the target matrix would corrupt memory.
    * 
    * Parameters:
    *     mat          Matrix to check
    * Returns:
    *     'true' if the evaluation of the expression to target matrix would corrupt memory 
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    constexpr bool _MATRIX_TIMES_SCALAR_<MAO,T>::unsafe_evaluation_target(_MATRIX_<MAT> const& mat) const 
    {
        return (matl::references_address(&mat,this->mat_) && transposed_matrix(this->mat_))
               || matl::unsafe_evaluation_target(mat,this->mat_);
    }

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const*)
    * 
    * Check if the expression references the address of the target matrix.
    * 
    * Parameters:
    *     mat          Matrix to check
    * Returns:
    *     'true' if the expression references the address of the target matrix
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    constexpr bool _MATRIX_TIMES_SCALAR_<MAO,T>::references_address(_MATRIX_<MAT> const* mat) const
    {return matl::references_address(mat,this->mat_);} 

   /****************************************************************************************
    * void eval(_MATRIX_<MAT_T>& mat) const
    * 
    * Evaluate the expression to the target matrix.
    * 
    * Parameters:
    *     mat          Target matrix
    ****************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    void _MATRIX_TIMES_SCALAR_<MAO,T>::eval(_MATRIX_<MAT>& mat) const
    {
        if (is_same_storage_order<MAT,MAO>::value)
        {
            const index_t N = this->size();
            for (index_t i = 1; i <= N; i++)
                mat(i) = (*this)(i);
        }
        else
        {
            const index_t M = this->rows();
            const index_t N = this->cols();
            for (index_t i = 1; i <= M; i++)
                for (index_t j = 1; j <= N; j++)
                    mat(i,j) = (*this)(i,j);
        }
    } 
}

#endif