#ifndef _MATL_MATRIX_DIVIDED_BY_SCALAR_IMPLEMENTATION__HPP
#define _MATL_MATRIX_DIVIDED_BY_SCALAR_IMPLEMENTATION__HPP

namespace matl
{
   /***********************************************************************************************
    *  Constructor (elementwise / scalar)
    * 
    * Constructor for matrix or elementwise expression divided by a scalar.
    * 
    * Parameters:
    *     mat       Matrix or elementwise expression
    *     scalar    Scalar
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAO_, 
             class T_,
             enable_if_elementwise<MAO_>,
             enable_if_scalar<T_> >
    _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::_MATRIX_DIVIDED_BY_SCALAR_(MAO_ const& mat, T_ scalar) : mat_(unqualified_cast(mat)),
                                                                                                scalar_(scalar)
    {}

   /***********************************************************************************************
    * Constructor (non-elementwise / scalar)
    * 
    * Constructor for non-elementwise expression divided by a scalar.
    * 
    * Parameters:
    *     expr      Non-elementwise expression
    *     scalar    Scalar
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAO_, 
             class T_,
             disable_if_elementwise<MAO_>,
             enable_if_scalar<T_> >
    _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::_MATRIX_DIVIDED_BY_SCALAR_(MAO_ const& expr, T_ scalar) : scalar_(scalar)
    {
        MemoryManager::allocate(this->mat_, expr.rows(), expr.cols());
        expr.eval(this->mat_);

        // Call destructor on temporary object
        const_cast<MAO_ &>(expr).~MAO_();
    }
    
    
   /***********************************************************************************************
    * index_t rows() const
    * 
    * Get number of rows of resultant matrix
    * 
    * Returns:
    *     index_t          Number of rows in resultant matrix
    ***********************************************************************************************/
    template<typename MAO, class T>
    index_t _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::rows() const
    {return this->mat_.rows();} 

   /***********************************************************************************************
    * index_t cols() const
    * 
    * Get number of columns of resultant matrix
    * 
    * Returns:
    *     index_t          Number of columns in resultant matrix
    * **********************************************************************************************/
    template<typename MAO, class T>
    index_t _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::cols() const
    {return this->mat_.cols();} 

   /***********************************************************************************************
    * index_t size() const
    * 
    * Get number of elements in resultant matrix
    * 
    * Returns:
    *     index_t          Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename MAO, class T>
    index_t _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::size() const
    {return this->mat_.size();} 

   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Get the result of the division operation at the i-th row and j-th column.
    * 
    * Parameters:
    *     i          Row index
    *     j          Column index
    * Returns:
    *     dtype      Result of division operation
    ***********************************************************************************************/
    template<typename MAO, class T>
    typename _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::dtype _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::operator() (index_t i, index_t j) const
    {return this->mat_(i,j) / this->scalar_;}

   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Get the result of the division operation at the i-th index.
    * 
    * Parameters:
    *     i          Index of the element to get
    * Returns:
    *     dtype      Result of division operation
    ***********************************************************************************************/
    template<typename MAO, class T>
    typename _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::dtype 
    _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::operator() (index_t i) const
    {return this->mat_(i) / this->scalar_;} 


   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const* target)
    * 
    * Check if the expression references the address of the target matrix.
    *
    * Parameters:
    *     target     Target matrix
    * Returns:
    *     'true' if evaluation of the expression would corrupt memory
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    constexpr bool _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {
        const bool referenced_in_transpose = matl::references_address(&target, this->mat_) && transposed_matrix(target);
        return referenced_in_transpose || matl::unsafe_evaluation_target(target,this->mat_);
    } 

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const* target)
    * 
    * Check if the expression references the address of the target matrix.
    *
    * Parameters:
    *     target     Target matrix
    * Returns:
    *     'true' if the expression references the address of the target matrix
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    constexpr bool _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::references_address(_MATRIX_<MAT> const* target) const
    {return matl::references_address(target, this->mat_);} 

   /***********************************************************************************************
    * void eval(_MATRIX_& target) const 
    * 
    * Evaluate the expression to the target matrix.
    *
    * Parameters:
    *     target     Target matrix
    ***********************************************************************************************/
    template<typename MAO, class T>
    template<typename MAT>
    void _MATRIX_DIVIDED_BY_SCALAR_<MAO,T>::eval(_MATRIX_<MAT>& target) const
    {
        if (is_same_storage_order<MAT,MAO>::value)
        {
            const index_t N = this->size();
            for (index_t i = 1; i <= N; i++)
                target(i) = (*this)(i);
        }
        else
        {
            const index_t M = this->rows();
            const index_t N = this->cols();
            for (index_t i = 1; i <= M; i++)
                for (index_t j = 1; j <= N; j++)
                    target(i,j) = (*this)(i,j);
        }
    }
}

#endif