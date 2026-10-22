#ifndef MATL_FIXEDSIZEMATRIX_IMPLEMENTATION_HPP
#define MATL_FIXEDSIZEMATRIX_IMPLEMENTATION_HPP

namespace matl
{
   /***********************************************************************************************
    * FixedSizeMatrix::operator= (DynamicMatrix<T2, ORDER> const&)
    * 
    * Matrix copy (DynamicMatrix with the same storage order, copied in storage order).
    * 
    * Parameters:
    *     mat        Matrix to copy
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    template<class T2>
    FixedSizeMatrix<T, M, N, ORDER>& 
    FixedSizeMatrix<T, M, N, ORDER>::operator= (DynamicMatrix<T2, ORDER> const& mat)
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
    * FixedSizeMatrix(matrix_rows<T>), operator= (matrix_rows<T>)
    *
    * Initializer-list construction and assignment, one inner list per row, e.g.
    *
    *     Matrix<double, 2, 3> F = {{1, 2, 3},
    *                               {4, 5, 6}};
    *
    * The list must have M rows of N elements (always checked; reported as matl::dimension_error).
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    FixedSizeMatrix<T, M, N, ORDER>::FixedSizeMatrix(matrix_rows<T> rows)
    {this->copy_rows(rows);}

    template<class T, index_t M, index_t N, StorageOrder ORDER>
    FixedSizeMatrix<T, M, N, ORDER>& FixedSizeMatrix<T, M, N, ORDER>::operator= (matrix_rows<T> rows)
    {
        this->copy_rows(rows);
        return *this;
    }

    template<class T, index_t M, index_t N, StorageOrder ORDER>
    void FixedSizeMatrix<T, M, N, ORDER>::copy_rows(matrix_rows<T> rows)
    {
        index_t rows_in_list, cols_in_list;
        detail::initializer_list_dimensions(rows, rows_in_list, cols_in_list);
        if ((rows_in_list != M) || (cols_in_list != N))
            detail::raise<dimension_error>("cannot initialize a " + detail::dimensions(M, N) +
                                           " fixed-size matrix from a " +
                                           detail::dimensions(rows_in_list, cols_in_list) + " initializer list");
        detail::copy_initializer_list(*this, rows);
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
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    template<ExpressionType EXPR>
    void FixedSizeMatrix<T, M, N, ORDER>::assign_expression(EXPR const& expression)
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
            FixedSizeMatrix<T, M, N, ORDER> temp;
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
    * Copy contents of generic matrix to instance FixedSizeMatrix (used by copy and
    * matrix assignment). The matrix must have dimensions M x N.
    * 
    * Parameters:
    *     mat        Matrix to copy
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    template<MatrixType MAT>
    void FixedSizeMatrix<T, M, N, ORDER>::copy_matrix(MAT const& mat)
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
    * FixedSizeMatrix::operator() (index_t)
    * 
    * Single index getter.
    * 
    * Parameters:
    *     index_t i        Index of the element to get
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    T FixedSizeMatrix<T, M, N, ORDER>::operator() (index_t i) const
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((i < 1) || (i > SIZE))
            detail::index_out_of_range(i, SIZE);
        #endif
        return this->data_[i-1];
    }
    
   /***********************************************************************************************
    * FixedSizeMatrix::operator() (index_t, index_t)
    * 
    * Double index getter.
    * 
    * Parameters:
    *     index_t i        Row index of the element to get
    *     index_t j        Column index of the element to get
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    T FixedSizeMatrix<T, M, N, ORDER>::operator() (index_t i, index_t j) const
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((i < 1) || (i > M) || (j < 1) || (j > N))
            detail::index_out_of_range(i, j, M, N);
        #endif
        return this->data_[offset(i, j)];
    } 

   /***********************************************************************************************
    * FixedSizeMatrix::operator() (index_t)
    *
    * Single index accessor.
    *
    * Parameters:
    *     index_t i        Index of the element to get
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    T& FixedSizeMatrix<T, M, N, ORDER>::operator() (index_t i)
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((i < 1) || (i > SIZE))
            detail::index_out_of_range(i, SIZE);
        #endif
        return this->data_[i-1];
    } 

   /***********************************************************************************************
    * FixedSizeMatrix::operator() (index_t, index_t)
    * 
    * Double index accessor.
    * 
    * Parameters:
    *     index_t i        Row index of the element to get
    *     index_t j        Column index of the element to get
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    T& FixedSizeMatrix<T, M, N, ORDER>::operator() (index_t i, index_t j)
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((i < 1) || (i > M) || (j < 1) || (j > N))
            detail::index_out_of_range(i, j, M, N);
        #endif
        return this->data_[offset(i, j)];
    }
    
   /***********************************************************************************************
    * FixedSizeMatrix::rows()
    * 
    * Returns:
    *    index_t M        Number of rows
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    index_t FixedSizeMatrix<T, M, N, ORDER>::rows() const
    {return M;}
    
   /***********************************************************************************************
    * FixedSizeMatrix::cols()
    * 
    * Column getter.
    * 
    * Returns:
    *    index_t N        Number of columns
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    index_t FixedSizeMatrix<T, M, N, ORDER>::cols() const
    {return N;} 

   /***********************************************************************************************
    * FixedSizeMatrix::size()
    *
    * Returns:
    *   index_t SIZE        Number of elements
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    index_t FixedSizeMatrix<T, M, N, ORDER>::size() const
    {return SIZE;} 


   /***********************************************************************************************
    * FixedSizeMatrix::operator= (S)
    * 
    * Scalar data writing.
    * 
    * Parameters:
    *     T val        Value to write
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    template<Scalar S>
    FixedSizeMatrix<T, M, N, ORDER>& FixedSizeMatrix<T, M, N, ORDER>::operator= (S val)
    {
        for (index_t i = 1; i <= SIZE; i++)
            (*this)(i) = val;
        return *this;
    }
     


   /***********************************************************************************************
    * void allocate(FixedSizeMatrix<T, M, N, ORDER>&, index_t, index_t)
    * 
    * The 'allocate' routine for 'FixedSizeMatrix' objects is present to support more
    * seamless polymorphism when using mixed matrix types. Because the size of the matrix is 
    * constant and determined at compile time, it cannot be altered. This function checks that
    * the requested size is the same as the size of the matrix and returns if it is. If the
    * requested size is different, a matl::dimension_error is reported (see _matl_errors_.hpp).
    *
    * Parameters:
    *     matrix        FixedSizeMatrix object
    *     rows          Number of rows in matrix
    *     cols          Number of columns in matrix
    ***********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    void allocate(FixedSizeMatrix<T, M, N, ORDER>& /*matrix*/, index_t rows, index_t cols)
    {
        if (rows != M || cols != N)
            detail::raise<dimension_error>("cannot allocate a " + detail::dimensions(M, N) + " fixed-size matrix as " +
                                           detail::dimensions(rows, cols));
    }
}
#endif