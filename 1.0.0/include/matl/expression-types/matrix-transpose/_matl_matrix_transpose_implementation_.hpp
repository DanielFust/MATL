#ifndef _MATL_MATRIX_TRANSPOSE_IMPLEMENTATION__HPP
#define _MATL_MATRIX_TRANSPOSE_IMPLEMENTATION__HPP

namespace matl 
{
   /***********************************************************************************************
    * Constructor (elementwise object)
    * 
    * Constructs a matrix transpose expression from an elementwise matrix algebra object.
    * 
    * Parameters:
    *     mat        matrix or elementwise expression
    ***********************************************************************************************/
    template<typename MAO>
    template<typename MAO_,
             enable_if_elementwise<MAO_> >
    _MATRIX_TRANSPOSE_<MAO>::_MATRIX_TRANSPOSE_(MAO_ const& mat) : mat_(unqualified_cast(mat))    
    {}
    
   /***********************************************************************************************
    * Constructor (non-elementwise object)
    * 
    * Constructs a matrix transpose expression from a non-elementwise matrix algebra object.
    * 
    * Parameters:
    *     mat        matrix or elementwise expression
    ***********************************************************************************************/
    template<typename MAO>
    template<typename MAO_,
             disable_if_elementwise<MAO_> >
    _MATRIX_TRANSPOSE_<MAO>::_MATRIX_TRANSPOSE_(MAO_ const& expr) 
    {
        // Allocate memory for the expression result
        MemoryManager::allocate(this->mat_, expr.rows(), expr.cols());

        // Evaluate the expression to the transpose
        expr.eval(this->mat_);

        // Call early destructor for temporary object
        static_cast<MAO_&>(const_cast<_EXPRESSION_<MAO_>&>(expr)).~MAO_();
    } 


   /***********************************************************************************************
    * index_t rows() const
    * 
    * Returns:
    *   index_t          Number of rows in resultant matrix
    ***********************************************************************************************/
    template<typename MAO>
    index_t _MATRIX_TRANSPOSE_<MAO>::rows() const
    {return mat_.cols();}
    
   /***********************************************************************************************
    * index_t cols() const
    *
    * Returns:
    *   index_t          Number of columns in resultant matrix
    ***********************************************************************************************/
    template<typename MAO>
    index_t _MATRIX_TRANSPOSE_<MAO>::cols() const
    {return mat_.rows();}
    
   /***********************************************************************************************
    * index_t size() const
    * 
    * Returns:
    *   index_t          Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename MAO>
    index_t _MATRIX_TRANSPOSE_<MAO>::size() const
    {return mat_.size();}

   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Get indexed result of expression
    *
    * Parameters:
    *     i                Row index
    *     j                Column index
    * Returns:
    *     dtype            Value of resultant matrix at index pair
    ***********************************************************************************************/
    template<typename MAO>
    typename _MATRIX_TRANSPOSE_<MAO>::dtype
    _MATRIX_TRANSPOSE_<MAO>::operator() (index_t i, index_t j) const
    {return this->mat_(j,i);} 

   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Get single-indexed result of expression
    *
    * Parameters:
    *     i                Index of the element to get
    ***********************************************************************************************/
    template<typename MAO>
    typename _MATRIX_TRANSPOSE_<MAO>::dtype
    _MATRIX_TRANSPOSE_<MAO>::operator() (index_t i) const
    {return this->mat_(i);} 
    
   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const&)
    * 
    * Check if evaluation target is safe
    *
    * Parameters:
    *     target       Candidate evalutation target matrix
    * Returns:
    *     'true' if direct evaluation corrupts memory 'false' otherwise
    ***********************************************************************************************/
    template<typename MAO>
    template<typename MAT>
    constexpr bool _MATRIX_TRANSPOSE_<MAO>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {
        // Pure transpose expressions are safe, however, arithmetic expressions involving the
        // transpose of the target matrix are not.
        return matl::unsafe_evaluation_target(target,this->mat_);
    }

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const*)
    *
    * Check if the expression references the address of the target matrix
    *
    * Parameters:
    *     mat         Matrix to check address of
    * Returns:
    *     'true' if the expression references the address of the target matrix
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO>
    template<typename MAT>
    constexpr bool _MATRIX_TRANSPOSE_<MAO>::references_address(_MATRIX_<MAT> const* mat) const
    {return matl::references_address(mat,this->mat_);}

   /***********************************************************************************************
    * void eval(_MATRIX_ const&) const
    * 
    * Evaluate the expression to the target matrix
    *
    * Parameters:
    *     target       Target matrix
    ***********************************************************************************************/
    template<typename MAO>
    template<typename MAT>
    void _MATRIX_TRANSPOSE_<MAO>::eval(_MATRIX_<MAT>& target) const
    {
        const index_t M = target.rows();
        const index_t N = target.cols();
        #ifdef MATL_BOUNDS_CHECK
        if ((M != this->mat_.cols()) || (N != this->mat_.rows()))
            throw std::out_of_range("Matrix transpose expression: target matrix dimensions do not match");
        #endif

        // Special case for square matrices
        // elements (i,j) and (j,i) evaluated simultaneously
        //
        // This covers cases in which the target matrix appears in the expression
        // (e.g. A = transpose(A+A)) as the target matrix and all matrices in the expression
        // are necessarily square. Any matrix multiplications (which may alter the size of the 
        // expression) will be evaluated during construction of the _MATRIX_TRANSPOSE_ object
        // or its constituents and so will be memory safe.
        if (M==N) 
        {
            for (index_t i = 1; i <= M; i++)
            {
                target(i,i) = (*this)(i,i);
                for (index_t j = i+1; j <= N; j++)
                {
                    dtype temp_ij = this->mat_(i,j);
                    dtype temp_ji = this->mat_(j,i);

                    target(i,j) = temp_ji;
                    target(j,i) = temp_ij;
                }   
            }        
        }
        else // (M!=N)
        {
            if (is_same_storage_order<MAT,MAO>::value)
            {
                const index_t K = this->size();
                for (index_t i = 1; i <= K; i++)
                    target(i) = (*this)(i);
            }
            else
            {
                // Copy elements to target matrix
                for (index_t i = 1; i <= M; i++)
                    for (index_t j = 1; j <= N; j++)
                        target(i,j) = (*this)(i,j);
            }
        }
    } 
}

#endif