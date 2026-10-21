/**************************************************************************************************
 * File: _matl_scalar_minus_matrix_implementation_.hpp
 * 
 * Description: Implementations for the scalar minus matrix expression template.
 **************************************************************************************************/
namespace matl
{
   /***********************************************************************************************
    * Constructor (elementwise object)
    * 
    * Initialize scalar-matrix subtraction expression
    * Parameters:
    *    mat              Matrix or elementwise expression
    *    scalar           Scalar
     ***********************************************************************************************/ 
    template<class T, typename MAO>
    template<typename MAO_, 
             enable_if_elementwise<MAO_>, 
             enable_if_scalar<T> >
    _SCALAR_MINUS_MATRIX_<T,MAO>::_SCALAR_MINUS_MATRIX_(T scalar, MAO_ const& mat) : 
        mat_(unqualified_cast(mat)),
        scalar_(scalar)
    {}

   /***********************************************************************************************
    * Constructor (non-elementwise object)
    * 
    * Initialize scalar-matrix subtraction expression
    * 
    * Parameters:
    *   expr             non-elementwise expression
    *   scalar           Scalar
    ***********************************************************************************************/  
    template<class T, typename MAO>
    template<typename MAO_, 
             disable_if_elementwise<MAO_>, 
             enable_if_scalar<T> >
    _SCALAR_MINUS_MATRIX_<T,MAO>::_SCALAR_MINUS_MATRIX_(T scalar, MAO_ const& expr) : scalar_(scalar)
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
    *     index_t          Number of rows in resultant matrix
    **********************************************************************************************/ 
    template<class T, typename MAO>
    index_t _SCALAR_MINUS_MATRIX_<T,MAO>::rows() const
    {return this->mat_.rows();}
      
   /***********************************************************************************************
    * index_t cols() const
    * 
    * Get number of columns of resultant matrix
    * 
    * Number of columns in expression. Downcasts to derived type.
    * Returns:
    *     index_t          Number of columns in resultant matrix
    ***********************************************************************************************/
    template<class T, typename MAO>
    index_t _SCALAR_MINUS_MATRIX_<T,MAO>::cols() const
    {return this->mat_.cols();}   
      
   /***********************************************************************************************
    * index_t size() const
    * 
    * Get number of elements in of resultant matrix
    * 
    * Number of elements in expression. Downcasts to derived type.
    * Returns:
    *     index_t          Number of elements in resultant matrix
    ***********************************************************************************************/ 
    template<class T, typename MAO>
    index_t _SCALAR_MINUS_MATRIX_<T,MAO>::size() const
    {return this->mat_.size();}
     
   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Retrieve the resultant value at specified index pair
    * 
    * Parameters:
    *     i          Row index
    *     j          Column index 
    * Returns:
    *     dtype      Value of resultant matrix at index pair  
    ***********************************************************************************************/
    template<class T, typename MAO>
    typename _SCALAR_MINUS_MATRIX_<T,MAO>::dtype 
    _SCALAR_MINUS_MATRIX_<T,MAO>::operator() (index_t i, index_t j) const
    {return this->scalar_ - this->mat_(i,j);}

   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Retrieve the resultant value at specified index
    * 
    * Parameters:
    *     i          Index of the element to get
    * Returns:
    *     dtype      Value of resultant matrix at index
    ***********************************************************************************************/
    template<class T, typename MAO>
    typename _SCALAR_MINUS_MATRIX_<T,MAO>::dtype 
    _SCALAR_MINUS_MATRIX_<T,MAO>::operator() (index_t i) const
    {return this->scalar_ - this->mat_(i);} 

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const* target)
    * 
    * Check if evaluation target is safe
    * 
    * Parameters:
    *     target       Candidate target matrix for evaluation of expression
    * Returns:
    *     'true' if the evaluation of the expression would corrupt memory.
    *     'false' otherwise.
    ***********************************************************************************************/
    template<class T, typename MAO>
    template<typename MAT>
    constexpr bool _SCALAR_MINUS_MATRIX_<T,MAO>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {
        return (matl::references_address(&target,this->mat_) && transposed_matrix(this->mat_)) 
               || matl::unsafe_evaluation_target(target,this->mat_);
    }

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const* target)
    * 
    * Check if the expression references the address of the target matrix
    * 
    * Parameters:
    *     target       Matrix to check address of
    * Returns:
    *     'true' if the expression references the address of the target matrix
    *     'false' otherwise
    ***********************************************************************************************/
    template<class T, typename MAO>
    template<typename MAT>
    constexpr bool _SCALAR_MINUS_MATRIX_<T,MAO>::references_address(_MATRIX_<MAT> const* target) const
    {return matl::references_address(target,this->mat_);} 

   /***********************************************************************************************
    * void eval(<MAT>& target) const
    * 
    * Evaluate the expression to the target matrix
    * 
    * Parameters:
    *    target      Target matrix
    ***********************************************************************************************/
    template<class T, typename MAO>
    template<typename MAT>
    void _SCALAR_MINUS_MATRIX_<T,MAO>::eval(_MATRIX_<MAT>& target) const
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