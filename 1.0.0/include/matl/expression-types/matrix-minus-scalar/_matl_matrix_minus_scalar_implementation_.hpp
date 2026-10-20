/**************************************************************************************************
 * File: _matl_matrix_minus_scalar_implementation_.hpp
 * 
 * Description: Implementations for the matrix minus scalar expression.
 **************************************************************************************************/
namespace matl
{
   /***********************************************************************************************
    * Constructor (elementwise)
    * 
    * Initialize matrix minus scalar expression
    * Parameters:
    *    mat              Matrix or elementwise expression
    *    scalar           Scalar
    ***********************************************************************************************/ 
    template<typename MAO, class T>
    template<typename MAO_, 
             class T_,
             enable_if_elementwise<MAO_>,
             enable_if_scalar<T_> >
    _MATRIX_MINUS_SCALAR_<MAO,T>::_MATRIX_MINUS_SCALAR_(MAO_ const& mat, T_ scalar) : mat_(unqualified_cast(mat)),
                                                                                      scalar_(scalar)
    {}


   /***********************************************************************************************
    * Constructor (non-elementwise)
    * 
    * Initialize matrix minus scalar expression
    * Parameters:
    *    expr           Non-elementwise expression
    *    scalar         Scalar
    ***********************************************************************************************/ 
    template<typename MAO, class T>
    template<typename MAO_, 
             class T_,
             disable_if_elementwise<MAO_>,
             enable_if_scalar<T_> >
    _MATRIX_MINUS_SCALAR_<MAO,T>::_MATRIX_MINUS_SCALAR_(MAO_ const& expr, T_ scalar) : scalar_(scalar)
    {
        MemoryManager::allocate(this->mat_, expr.rows(), expr.cols());
        expr.eval(this->mat_);
    } 
    
   /***********************************************************************************************
    * index_t rows() const
    *
    * Get number of rows of resultant matrix
    * Number of rows in expression. Downcasts to derived type.
    * 
    * Returns:
    *    index_t          Number of rows in resultant matrix
    **********************************************************************************************/ 
    template<typename MAO, class T>
    index_t _MATRIX_MINUS_SCALAR_<MAO,T>::rows() const
    {return this->mat_.rows();}
     
   /***********************************************************************************************
    * index_t cols() const
    * 
    * Get number of columns of resultant matrix
    * 
    * Number of columns in expression. Downcasts to derived type.
    * Returns:
    *   index_t          Number of columns in resultant matrix
    ***********************************************************************************************/
    template<typename MAO, class T>
    index_t _MATRIX_MINUS_SCALAR_<MAO,T>::cols() const
    {return this->mat_.cols();}   
     
   /***********************************************************************************************
    * index_t size() const
    * 
    * Get number of elements in of resultant matrix
    * 
    * Number of elements in expression. Downcasts to derived type.
    * Returns:
    *   index_t          Number of elements in resultant matrix
    ***********************************************************************************************/ 
    template<typename MAO, class T>
    index_t _MATRIX_MINUS_SCALAR_<MAO,T>::size() const
    {return this->mat_.size();}
    
   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Retrieve the resultant value at specified index pair
    * 
    * Template Parameters:
    *     MAO           Atomic expression type
    * 
    * Parameters:
    *     i                Row index
    *     j                Column index 
    * Returns:
    *     dtype            Value of resultant matrix at index pair
    ************************************************************************************************/
    template<typename MAO, class T>
    typename _MATRIX_MINUS_SCALAR_<MAO,T>::dtype
    _MATRIX_MINUS_SCALAR_<MAO,T>::operator() (index_t i, index_t j) const
    {return this->mat_(i,j) - this->scalar_;}

   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Retrieve the resultant value at specified index
    * 
    * Template Parameters:
    *     MAO           Atomic expression type
    * 
    * Parameters:
    *     i                Index of the element to get
    * Returns:
    *     dtype            Value of resultant matrix at index
    ************************************************************************************************/
    template<typename MAO, class T>
    typename _MATRIX_MINUS_SCALAR_<MAO,T>::dtype
    _MATRIX_MINUS_SCALAR_<MAO,T>::operator() (index_t i) const
    {return this->mat_(i) - this->scalar_;} 

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const*)
    * 
    * Check if evaluation of the expression to the target matrix would corrupt memory.
    * 
    * Parameters:
    *     mat            Matrix to check
    * Returns:
    *     'true' if the evaluation of the expression to target matrix would corrupt memory 
    *     'false' otherwise
    ***********************************************************************************************/ 
    template<typename MAO, class T>
    template<typename MAT>
    constexpr bool _MATRIX_MINUS_SCALAR_<MAO,T>::unsafe_evaluation_target(_MATRIX_<MAT> const& mat) const 
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
    *     mat            Matrix to check
    * Returns:
    *     'true' if the expression references the address of the target matrix
    *     'false' otherwise
    ***********************************************************************************************/ 
    template<typename MAO, class T>
    template<typename MAT>
    constexpr bool _MATRIX_MINUS_SCALAR_<MAO,T>::references_address(_MATRIX_<MAT> const* mat) const
    {return matl::references_address(mat,this->mat_);}
    


   /***********************************************************************************************
    *  void eval(MAT& target) const
    * 
    * Evaluate the expression to the target matrix
    * Parameters:
    *     target           Target matrix
    ***********************************************************************************************/ 
    template<typename MAO, class T>
    template<typename MAT>
    void _MATRIX_MINUS_SCALAR_<MAO,T>::eval(_MATRIX_<MAT>& target) const
    {
        if (is_same_storage_order<MAT,MAO>::value)
        {
            const index_t N = this->size();
            for (index_t i = 1; i <= N; ++i)
                target(i) = (*this)(i);
        }
        else
        {
            const index_t M = this->rows();
            const index_t N = this->cols();
            for (index_t i = 1; i <= M; ++i)
                for (index_t j = 1; j <= N; ++j)
                    target(i,j) = (*this)(i,j);
        }
    }
}
