#ifndef MATL_FIXEDSIZEROWMAJORMATRIX_IMPLEMENTATION_HPP
#define MATL_FIXEDSIZEROWMAJORMATRIX_IMPLEMENTATION_HPP

namespace matl
{
   /***********************************************************************************************
    * FixedSizeRowMajorMatrix::operator= (RowMajorMatrix<T2> const&)
    * 
    * Matrix copy (RowMajorMatrix specialization).
    * 
    * Parameters:
    *     RowMajorMatrix mat        Matrix to copy
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    template<class T2>
    FixedSizeRowMajorMatrix<T, M, N>& 
    FixedSizeRowMajorMatrix<T, M, N>::operator= (RowMajorMatrix<T2> const& mat)
    {
        #ifdef MATL_DIMENSION_CHECK
        if ((this->rows() != mat.rows())|| (this->cols() != mat.cols()))
            detail::assignment_dimensions(M, N, mat.rows(), mat.cols());
        #endif
        for(index_t i=0; i<SIZE; i++)
            this->data_[i] = mat(i+1);
        return *this;
    } 

   /***********************************************************************************************
    * void assign_expression(EXPR const&)
    *
    * Evaluates an expression into this matrix (used by the expression constructor and
    * assignment). The expression must have dimensions M x N.
    *
    * Parameters:
    *     expression      Matrix expression
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    template<ExpressionType EXPR>
    void FixedSizeRowMajorMatrix<T, M, N>::assign_expression(EXPR const& expression)
    {
        // The operands of every operation in the expression must fit
        #ifdef MATL_DIMENSION_CHECK
        detail::check_dimensions(expression);
        #endif

        #ifdef MATL_DIMENSION_CHECK
        if ((this->rows() != expression.rows())|| (this->cols() != expression.cols()))
            detail::assignment_dimensions(M, N, expression.rows(), expression.cols());
        #endif

        // Check (at compile-time) if the evaluation is safe 
        if (expression.unsafe_evaluation_target(*this))
        {
            // Use temporary matrix to store the result
            FixedSizeRowMajorMatrix<T, M, N> temp;
            allocate(temp, expression.rows(), expression.cols());
            expression.eval(temp);

            // Copy data to this matrix (I believe C-style fixed size arrays are not movable)
            this->copy_matrix(temp);
        }
        else // evaluation to this matrix is safe
        {
            // If Matrix is not allocated yet, allocate it
            if (this->size() == 0) 
                allocate(*this, expression.rows(), expression.cols());
        
            // Evaluate expression    
            expression.eval(*this);
        }
    } 

    

   /***********************************************************************************************
    * void copy_matrix(MAT const&)
    *   
    * Copy contents of generic matrix to instance FixedSizeRowMajorMatrix (used by copy and
    * matrix assignment). The matrix must have dimensions M x N.
    * 
    * Parameters:
    *     mat        Matrix to copy
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    template<MatrixType MAT>
    void FixedSizeRowMajorMatrix<T, M, N>::copy_matrix(MAT const& mat)
    {
        #ifdef MATL_DIMENSION_CHECK
        if ((this->rows() != mat.rows())|| (this->cols() != mat.cols()))
            detail::assignment_dimensions(M, N, mat.rows(), mat.cols());
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
            detail::index_out_of_range(i, SIZE);
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
            detail::index_out_of_range(i, j, M, N);
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
            detail::index_out_of_range(i, SIZE);
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
            detail::index_out_of_range(i, j, M, N);
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
    * FixedSizeRowMajorMatrix::operator= (S)
    * 
    * Scalar data writing.
    * 
    * Parameters:
    *     T val        Value to write
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    template<Scalar S>
    FixedSizeRowMajorMatrix<T, M, N>& FixedSizeRowMajorMatrix<T, M, N>::operator= (S val)
    {
        for (index_t i = 1; i <= SIZE; i++)
            (*this)(i) = val;
        return *this;
    }
     


   /***********************************************************************************************
    * void allocate(FixedSizeRowMajorMatrix<T, M, N>&, index_t, index_t)
    * 
    * The 'allocate' routine for 'FixedSizeRowMajorMatrix' objects is present to support more
    * seamless polymorphism when using mixed matrix types. Because the size of the matrix is 
    * constant and determined at compile time, it cannot be altered. This function checks that
    * the requested size is the same as the size of the matrix and returns if it is. If the
    * requested size is different, a matl::dimension_error is reported (see _matl_errors_.hpp).
    *
    * Parameters:
    *     matrix        FixedSizeRowMajorMatrix object
    *     rows          Number of rows in matrix
    *     cols          Number of columns in matrix
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    void allocate(FixedSizeRowMajorMatrix<T, M, N>& /*matrix*/, index_t rows, index_t cols)
    {
        if (rows != M || cols != N)
            detail::raise<dimension_error>("cannot allocate a " + detail::dimensions(M, N) + " fixed-size matrix as " +
                                           detail::dimensions(rows, cols));
    }
}
#endif