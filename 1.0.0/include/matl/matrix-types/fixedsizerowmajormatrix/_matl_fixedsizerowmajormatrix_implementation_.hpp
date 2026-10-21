#ifndef _MATL_FIXEDSIZEROWMAJORMATRIX_IMPLEMENTATION__HPP
#define _MATL_FIXEDSIZEROWMAJORMATRIX_IMPLEMENTATION__HPP

namespace matl
{
   /***********************************************************************************************
    * FixedSizeRowMajorMatrix(_EXPRESSION_<EXPR> const&)
    * 
    * Expression constructor. Evaluates the expression into the new matrix. The expression must
    * have dimensions M x N.
    * 
    * Parameters:
    *     expr           Expression to evaluate
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    template<typename EXPR>
    FixedSizeRowMajorMatrix<T, M, N>::FixedSizeRowMajorMatrix(_EXPRESSION_<EXPR> const& expr)
    {*this = expr;}

   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::operator= (_MATRIX_&)
    * 
    * General matrix copy.
    * 
    * Parameters:
    *     _MATRIX_ mat        Matrix to copy
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    template<typename MAT>
    FixedSizeRowMajorMatrix<T, M, N>&
    FixedSizeRowMajorMatrix<T, M, N>::operator= (_MATRIX_<MAT>& mat)
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((this->rows() != mat.rows())|| (this->cols() != mat.cols()))
            throw std::out_of_range("FixedSizeRowMajorMatrix::operator= (_MATRIX_&): Matrix dimensions do not match");
        #endif
        
        for (index_t i = 1; i<=M; i++)
            for (index_t j = 1; j<=N; j++)
                (*this)(i,j) = mat(i,j);
        return *this;
    }
    

   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::operator= (FixedSizeRowMajorMatrix&)
    * 
    * Matrix copy (FixedSizeRowMajorMatrix specialization).
    * 
    * Parameters:
    *     FixedSizeRowMajorMatrix mat        Matrix to copy
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    FixedSizeRowMajorMatrix<T, M, N>& 
    FixedSizeRowMajorMatrix<T, M, N>::operator= (FixedSizeRowMajorMatrix<T, M, N>& mat)
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((this->rows() != mat.rows())|| (this->cols() != mat.cols()))
            throw std::out_of_range("FixedSizeRowMajorMatrix::operator= (FixedSizeRowMajorMatrix&): Matrix dimensions do not match");
        #endif
        for(index_t i=0; i<SIZE; i++)
            this->data_[i] = mat.data_[i];
        return *this;
    } 

   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::operator= (RowMajorMatrix&)
    * 
    * Matrix copy (RowMajorMatrix specialization).
    * 
    * Parameters:
    *     RowMajorMatrix mat        Matrix to copy
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    template<class T2>
    FixedSizeRowMajorMatrix<T, M, N>& 
    FixedSizeRowMajorMatrix<T, M, N>::operator= (RowMajorMatrix<T2>& mat)
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((this->rows() != mat.rows())|| (this->cols() != mat.cols()))
            throw std::out_of_range("FixedSizeRowMajorMatrix::operator= (RowMajorMatrix&): Matrix dimensions do not match");
        #endif
        for(index_t i=0; i<SIZE; i++)
            this->data_[i] = mat(i+1);
        return *this;
    } 

   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::operator= (_EXPRESSION_&)
    *
    * Trigger expression evaluation to this matrix.
    *
    * Parameters:
    *     _EXPRESSION_ expr      Matrix expression
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    template<typename EXPR>
    FixedSizeRowMajorMatrix<T, M, N>&
    FixedSizeRowMajorMatrix<T, M, N>::operator=(_EXPRESSION_<EXPR> const& expression)
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((this->rows() != expression.rows())|| (this->cols() != expression.cols()))
            throw std::out_of_range("FixedSizeRowMajorMatrix::operator= (_EXPRESSION_&): Matrix dimensions do not match");
        #endif

        // Check (at compile-time) if the evaluation is safe 
        if (expression.unsafe_evaluation_target(*this))
        {
            // Use temporary matrix to store the result
            FixedSizeRowMajorMatrix<T, M, N> temp;
            MemoryManager::allocate(temp, expression.rows(), expression.cols());
            expression.eval(temp);

            // Copy data to this matrix (I believe C-style fixed size arrays are not movable)
            this->copy(temp);
        }
        else // evaluation to this matrix is safe
        {
            // If Matrix is not allocated yet, allocate it
            if (this->size() == 0) 
                MemoryManager::allocate(*this,expression.rows(), expression.cols());
        
            // Evaluate expression    
            expression.eval(*this);
        }

        return *this;
    } 

    

   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::copy(_MATRIX_)
    *   
    * Copy contents of generic matrix to instance FixedSizeRowMajorMatrix. 
    * 
    * Parameters:
    *     _MATRIX_ mat        Matrix to copy
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    template<typename MAT>
    void FixedSizeRowMajorMatrix<T, M, N>::copy(const _MATRIX_<MAT>& mat)
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((this->rows() != mat.rows())|| (this->cols() != mat.cols()))
            throw std::out_of_range("FixedSizeRowMajorMatrix::copy(_MATRIX_): Matrix dimensions do not match");
        #endif
        
        for (index_t i = 1; i<=M; i++)
            for (index_t j = 1; j<=N; j++)
                (*this)(i,j) = mat(i,j);
    }

   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::operator() (index_t)
    * 
    * Single index getter.
    * 
    * Parameters:
    *     index_t i        Index of the element to get
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    T FixedSizeRowMajorMatrix<T, M, N>::operator() (index_t i) const
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((i < 1) || (i > SIZE))
            throw std::out_of_range("FixedSizeRowMajorMatrix::operator() (index_t): Index out of range");
        #endif
        return this->data_[i-1];
    }
    
   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::operator() (index_t, index_t)
    * 
    * Double index getter.
    * 
    * Parameters:
    *     index_t i        Row index of the element to get
    *     index_t j        Column index of the element to get
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    T FixedSizeRowMajorMatrix<T, M, N>::operator() (index_t i, index_t j) const
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((i < 1) || (i > M) || (j < 1) || (j > N))
            throw std::out_of_range("FixedSizeRowMajorMatrix::operator() (index_t, index_t): Index out of range");
        #endif
        return this->data_[(i-1)*N + (j-1)];
    } 

   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::operator() (index_t)
    *
    * Single index accessor.
    *
    * Parameters:
    *     index_t i        Index of the element to get
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    T& FixedSizeRowMajorMatrix<T, M, N>::operator() (index_t i)
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((i < 1) || (i > SIZE))
            throw std::out_of_range("FixedSizeRowMajorMatrix::operator() (index_t): Index out of range");
        #endif
        return this->data_[i-1];
    } 

   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::operator() (index_t, index_t)
    * 
    * Double index accessor.
    * 
    * Parameters:
    *     index_t i        Row index of the element to get
    *     index_t j        Column index of the element to get
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    T& FixedSizeRowMajorMatrix<T, M, N>::operator() (index_t i, index_t j)
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((i < 1) || (i > M) || (j < 1) || (j > N))
            throw std::out_of_range("FixedSizeRowMajorMatrix::operator() (index_t, index_t): Index out of range");
        #endif
        return this->data_[(i-1)*N + (j-1)];
    }
    
   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::rows()
    * 
    * Returns:
    *    index_t M        Number of rows
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    index_t FixedSizeRowMajorMatrix<T, M, N>::rows() const
    {return M;}
    
   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::cols()
    * 
    * Column getter.
    * 
    * Returns:
    *    index_t N        Number of columns
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    index_t FixedSizeRowMajorMatrix<T, M, N>::cols() const
    {return N;} 

   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::size()
    *
    * Returns:
    *   index_t SIZE        Number of elements
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    index_t FixedSizeRowMajorMatrix<T, M, N>::size() const
    {return SIZE;} 


   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::operator= (T)
    * 
    * Scalar data writing.
    * 
    * Parameters:
    *     T val        Value to write
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    template<class T2, enable_if_scalar<T2> >
    void FixedSizeRowMajorMatrix<T, M, N>::operator= (T2 val)
    {
        for (index_t i = 1; i <= SIZE; i++)
            (*this)(i) = val;
    }
     
}
#endif