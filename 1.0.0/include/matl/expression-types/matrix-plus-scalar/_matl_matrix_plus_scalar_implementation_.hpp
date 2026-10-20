/**************************************************************************************************
 * File: _matl_matrix_plus_scalar_implementation_.hpp
 * 
 * Description: Implementations for the matrix plus scalar expression.
 **************************************************************************************************/
namespace matl
{
   /***********************************************************************************************
    * Constructor (elementwise + scalar)
    * 
    * Initialize matrix plus scalar expression
    * Parameters:
    *    mat            Matrix or elementwise expression
    *    scalar         Scalar
    ***********************************************************************************************/ 
    template<typename MAO, class T>
    template<typename MAO_, 
             class T_,
             enable_if_elementwise<MAO_>,
             enable_if_scalar<T_> >
    _MATRIX_PLUS_SCALAR_<MAO,T>::_MATRIX_PLUS_SCALAR_(MAO_ const& mat, T_ scalar) : 
        mat_(unqualified_cast(mat)), 
        scalar_(scalar)
    {}

   /***********************************************************************************************
    * Constructor (non-elementwise + scalar)
    * 
    * Initialize matrix plus scalar expression from a non-elementwise expression and a scalar
    * 
    * Parameters:
    *     expr           Non-elementwise expression
    *     scalar         Scalar 
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAO_, 
             class T_,
             disable_if_elementwise<MAO_>,
             enable_if_scalar<T_> >
    _MATRIX_PLUS_SCALAR_<MAO,T>::_MATRIX_PLUS_SCALAR_(MAO_ const& expr, T_ scalar) : scalar_(scalar)
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
    index_t _MATRIX_PLUS_SCALAR_<MAO,T>::rows() const
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
    index_t _MATRIX_PLUS_SCALAR_<MAO,T>::cols() const
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
    index_t _MATRIX_PLUS_SCALAR_<MAO,T>::size() const
    {return this->mat_.size();}

   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Get the element at the i-th row and j-th column of the resultant matrix.
    * Parameters:
    *    i                Row index
    *    j                Column index
    * Returns:
    *   dtype            Element at the i-th row and j-th column
    ***********************************************************************************************/
    template<typename MAO, class T>
    typename _MATRIX_PLUS_SCALAR_<MAO,T>::dtype 
    _MATRIX_PLUS_SCALAR_<MAO,T>::operator() (index_t i, index_t j) const
    {return this->mat_(i,j) + this->scalar_;}

   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Get the element at the i-th index of the resultant matrix.
    * 
    * Parameters:
    *     i                Index of the element to get 
    * Returns:
    *     dtype            Element at the i-th index
    *************************************************************************************************/
    template<typename MAO, class T>
    typename _MATRIX_PLUS_SCALAR_<MAO,T>::dtype
    _MATRIX_PLUS_SCALAR_<MAO,T>::operator() (index_t i) const
    {return this->mat_(i) + this->scalar_;} 

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const* target)
    *
    * Parameters:
    *     target       Candidate target matrix for evaluation of expression
    * Returns:
    *     'true' if the evaluation of the expression would corrupt memory.
    *     'false' otherwise.
    ************************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    constexpr bool _MATRIX_PLUS_SCALAR_<MAO,T>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {
        return (matl::references_address(&target,this->mat_) && transposed_matrix(this->mat_)) 
               || matl::unsafe_evaluation_target(target,this->mat_);
    } 

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const* target)
    * 
    * Check if the expression references the address of the target matrix.
    * 
    * Parameters:
    *     target       Candidate target matrix for evaluation of expression
    * Returns:
    *     'true' if the expression references the address of the target matrix.
    *     'false' otherwise.
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    constexpr bool _MATRIX_PLUS_SCALAR_<MAO,T>::references_address(_MATRIX_<MAT> const* target) const
    {return matl::references_address(target,this->mat_);} 

   /***********************************************************************************************
    *  void eval(MAT_T& target) const
    * 
    * Evaluate the expression to the target matrix
    * Parameters:
    *   target           Target matrix
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    void _MATRIX_PLUS_SCALAR_<MAO,T>::eval(_MATRIX_<MAT>& target) const
    {
        if (is_same_storage_order<MAT,MAO>::value)
        {
            const index_t N = this->mat_.size();
            for (index_t i=1; i<=N; ++i)
                target(i) = (*this)(i);
        }
        else
        {
            const index_t M = this->mat_.rows();
            const index_t N = this->mat_.cols();
            for (index_t i=1; i<=M; ++i)
                for (index_t j=1; j<=N; ++j)
                    target(i,j) = (*this)(i,j);
        }
    } 
}




