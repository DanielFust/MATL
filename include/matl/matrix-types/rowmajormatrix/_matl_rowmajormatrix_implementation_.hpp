namespace matl
{

   /*=============================================================================================*
    *                                      CONSTRUCTORS
    *=============================================================================================*/



   /***********************************************************************************************
    * RowMajorMatrix()
    * 
    * Default constructor 
    * 
    * Template Parameters:
    *     T     Matrix data type
    * 
    * Parameters:
    * 
    * Returns:
    *     RowMajorMatrix      Row major matrix object with zero rows and zero columns
    ***********************************************************************************************/ 
    template<class T>
    RowMajorMatrix<T>::RowMajorMatrix() : rows_(0), 
                                          cols_(0), 
                                          data_(nullptr),
                                          external_mem_(false)
    {}



   /***********************************************************************************************
    * RowMajorMatrix(index_t rows, index_t cols)
    * 
    * Sized constructor 
    * 
    * Template Parameters:
    *     T     Matrix data type
    * 
    * Parameters:
    *     rows          number of rows
    *     cols          number of columns
    * 
    * Returns:
    *     RowMajorMatrix      Row major matrix object with specified rows and columns
    ***********************************************************************************************/ 
    template<class T>
    RowMajorMatrix<T>::RowMajorMatrix(index_t rows, index_t cols) : rows_(rows), 
                                                                    cols_(cols),
                                                                    data_(new T[rows*cols]),
                                                                    external_mem_(false)
    {}



   /***********************************************************************************************
    * RowMajorMatrix(index_t rows, index_t cols, T* data)
    * 
    * Preallocated constructor. Allows a RowMajorMatrix object to use preallocated memory (rows*cols
    * elements in row-major order).
    * 
    * This is intended to allow data for many persistent matrices to be allocated (or transferred
    * via MPI call) in a single contiguous block of memory.
    * 
    * Note: The memory remains owned by the caller and must outlive the matrix. It is not
    * deallocated by the RowMajorMatrix destructor.
    ***********************************************************************************************/
    template<class T>   
    RowMajorMatrix<T>::RowMajorMatrix(index_t rows, index_t cols, T* data) : rows_(rows), 
                                                                             cols_(cols), 
                                                                             data_(data),
                                                                             external_mem_(true)
    {}

   /***********************************************************************************************
    * RowMajorMatrix(MAT const&)
    * 
    * Copy constructor. Creates a new RowMajorMatrix object with the same dimensions and data as
    * the original matrix. The data is copied to a new array owned by the new matrix.
    * 
    * Parameters: 
    *    other          Matrix object to copy
    **************************************************************************************************/
    template<class T>
    template<MatrixType MAT>
    RowMajorMatrix<T>::RowMajorMatrix(MAT const& other) : rows_(other.rows()), 
                                                                    cols_(other.cols()),
                                                                    data_(new T[other.rows()*other.cols()]),
                                                                    external_mem_(false)
    {
        const index_t M = this->rows_;
        const index_t N = this->cols_;

        // Copy data from the original matrix
        for (index_t i = 1; i <= M; ++i)
            for (index_t j = 1; j <= N; ++j)
                (*this)(i,j) = other(i,j);
    } 

   /***********************************************************************************************
    * RowMajorMatrix(const RowMajorMatrix<T>&)
    * 
    * Copy constructor. Creates a new RowMajorMatrix object with the same dimensions and data as
    * the original matrix. The data is copied to a new array owned by the new matrix.
    * 
    * Parameters:
    *     other          RowMajorMatrix object to copy
    ***********************************************************************************************/ 
    template<class T>
    RowMajorMatrix<T>::RowMajorMatrix(const RowMajorMatrix<T>& other) : rows_(other.rows_), 
                                                                        cols_(other.cols_),
                                                                        data_(new T[other.rows_*other.cols_]),
                                                                        external_mem_(false)
    {
        // Copy data from the original matrix
        const index_t N = this->rows_ * this->cols_;
        for (index_t i = 0; i < N; ++i)
            this->data_[i] = other.data_[i];
    } 

   /***********************************************************************************************
    * RowMajorMatrix(EXPR const&)
    * 
    * Expression constructor. Creates a new RowMajorMatrix sized to the expression and evaluates
    * the expression into it.
    * 
    * Parameters:
    *     expr           Expression to evaluate
    ***********************************************************************************************/ 
    template<class T>
    template<ExpressionType EXPR>
    RowMajorMatrix<T>::RowMajorMatrix(EXPR const& expr) : rows_(0), 
                                                                        cols_(0),
                                                                        data_(nullptr),
                                                                        external_mem_(false)
    {*this = expr;}

   /***********************************************************************************************
    * RowMajorMatrix(matrix_rows<T>)
    * 
    * Initializer-list constructor: one inner list per row, e.g.
    * 
    *     Matrix<double> A = {{1, 2, 3},
    *                         {4, 5, 6}};      // 2x3
    * 
    * The rows must have the same length (always checked; reported as matl::dimension_error).
    * 
    * Parameters:
    *     rows           Rows of the matrix
    ***********************************************************************************************/ 
    template<class T>
    RowMajorMatrix<T>::RowMajorMatrix(matrix_rows<T> rows) : rows_(0), 
                                                            cols_(0),
                                                            data_(nullptr),
                                                            external_mem_(false)
    {
        index_t M, N;
        detail::initializer_list_dimensions(rows, M, N);
        allocate(*this, M, N);
        detail::copy_initializer_list(*this, rows);
    }

   /***********************************************************************************************
    * RowMajorMatrix(RowMajorMatrix<T>&&)
    * 
    * Move constructor. Takes the other matrix's data without copying (externally owned memory
    * stays externally owned). The other matrix is left empty (0x0).
    * 
    * Parameters:
    *     other          RowMajorMatrix object to move
    ***********************************************************************************************/ 
    template<class T>
    RowMajorMatrix<T>::RowMajorMatrix(RowMajorMatrix<T>&& other) noexcept : rows_(other.rows_), 
                                                                            cols_(other.cols_),
                                                                            data_(other.data_),
                                                                            external_mem_(other.external_mem_)
    {
        // Nullify other matrix so the data is not deleted
        other.data_         = nullptr;
        other.rows_         = 0;
        other.cols_         = 0;
        other.external_mem_ = false;
    } 



   /***********************************************************************************************
    * Destructor
    * 
    * Deletes the matrix data, unless it is externally owned.
    ***********************************************************************************************/ 
    template<class T>
    RowMajorMatrix<T>::~RowMajorMatrix()
    {if (!this->external_mem_) delete[] this->data_;}





   /*=============================================================================================*
    *                                   Assignment and Copying
    *=============================================================================================*/
  

  

   /***********************************************************************************************
    * RowMajorMatrix<T>& operator= (EXPR const&)
    * 
    * Triggers the evaluation of an expression and assigns the result to the RowMajorMatrix. If the
    * RowMajorMatrix is not allocated, it will be allocated to the size of the expression, 
    * otherwise it must have the same dimensions as the expression.
    * 
    * Parameters:
    *   expression          Expression to evaluate
    ***********************************************************************************************/
    template<class T>
    template<ExpressionType EXPR>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (EXPR const& expression)
    {
        #ifdef MATL_DEBUG
        std::cout << "RowMajorMatrix<T>& operator= (EXPR const&)" << std::endl;
        #endif

        // Dimension check debug flag:
        // The operands of every operation in the expression must fit. Only allow dynamic
        // allocation if the RowMajorMatrix has a size of 0. This is intended to discourage
        // repeated hidden reallocations 
        #ifdef MATL_DIMENSION_CHECK
        detail::check_dimensions(expression);
        if (  ((this->rows_ != expression.rows()) || (this->cols_ != expression.cols()))
            && (this->size() != 0) ) 
        {
            detail::assignment_dimensions(this->rows_, this->cols_, expression.rows(), expression.cols());
        }    
        #endif

        // Check (at compile-time) if the evaluation is safe 
        if (expression.unsafe_evaluation_target(*this))
        {
            // Evaluate into a temporary, then copy the result. The temporary's memory is
            // recycled (see EvaluationTemporary), so it is copied rather than moved, and this
            // matrix keeps its own memory.
            detail::EvaluationTemporary< RowMajorMatrix<T> > temp(expression.rows(), expression.cols());
            expression.eval(temp.matrix());
            this->copy(temp.matrix());
        }
        else // evaluation to this matrix is safe
        {
            // If Matrix is not allocated yet, allocate it
            if (this->size() == 0) 
                allocate(*this, expression.rows(), expression.cols());
        
            // Evaluate expression    
            expression.eval(*this);
        }
                
        return *this;    
    }






   /***********************************************************************************************
    * RowMajorMatrix<T>& void operator=(RowMajorMatrix<T>&&)
    * 
    * Same type RowMajorMatrix rvalue assignment. Takes ownership rvalue's data without copying.
    * 
    * Template Parameters:
    * Parameters:
    *     mat          RowMajorMatrix
    * 
    ***********************************************************************************************/
    template<class T>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T>&& mat) noexcept
    {
        #ifdef MATL_DEBUG
        std::cout << "RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T>&&)" 
                    << std::endl;
        #endif
        // A matrix over external memory keeps its buffer, so the data is copied into it and the
        // dimensions must match. Otherwise the matrix takes the other matrix's data, whatever its
        // size (as when matrices are moved around in a std::vector). This function is noexcept,
        // so a dimension error here terminates the program.
        #ifdef MATL_DIMENSION_CHECK
        if (this->external_mem_ && ((this->rows_ != mat.rows()) || (this->cols_ != mat.cols())))
            detail::assignment_dimensions(this->rows_, this->cols_, mat.rows(), mat.cols());
        #endif

        // Check for self assignment
        if (this == &mat) return *this; 

        
        // If memory is externally owned, the pointer should remain the same and
        // data should be copied rather than invoking move semantics
        if (this->external_mem()) 
        {
            const index_t N = this->size();
            for (index_t i = 0; i < N; ++i)
                this->data_[i] = mat.data_[i];   
                
            // Nullify other matrix and return early
            deallocate(mat);
    
            return *this;
        }

        // Deallocate this matrix  
        deallocate(*this);

        // Assign data pointer to pointer of other matrix data
        this->data_         = mat.data_;
        this->rows_         = mat.rows_;
        this->cols_         = mat.cols_;
        this->external_mem_ = mat.external_mem_;
        #ifdef MATL_DEBUG
        std::cout << "RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T>&&) - data moved" << std::endl;
        #endif
        
        // Nullify other matrix so the data is not deleted
        mat.data_         = nullptr;
        mat.rows_         = 0;
        mat.cols_         = 0;
        mat.external_mem_ = false;

        return *this;
    }



   /***********************************************************************************************
    *  RowMajorMatrix<T>& operator=(RowMajorMatrix<T> const&)
    * 
    * Copy assignment for RowMajorMatrix of the same type. Copies the data into this matrix,
    * allocating it first if it is empty.
    * 
    * Parameters:
    *     mat         Other RowMajorMatrix
    ***********************************************************************************************/
    template<class T>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T> const& mat)
    {
        #ifdef MATL_DEBUG
        std::cout << "RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T> const&) called" << std::endl;
        #endif

        // Check for self assignment (copy checks the dimensions)
        if (this != &mat) this->copy(mat);

        return *this;
    }



   /***********************************************************************************************
    *  RowMajorMatrix<T>& operator=(matrix_rows<T>)
    * 
    * Initializer-list assignment, one inner list per row: A = {{1, 2}, {3, 4}}. An empty matrix is
    * allocated to the size of the list; otherwise the list must have the matrix's dimensions.
    * A = {} empties the matrix. The dimensions are always checked (reported as
    * matl::dimension_error).
    * 
    * Parameters:
    *     rows        Rows of the matrix
    ***********************************************************************************************/
    template<class T>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (matrix_rows<T> rows)
    {
        if (rows.size() == 0)
        {
            deallocate(*this);
            return *this;
        }

        index_t M, N;
        detail::initializer_list_dimensions(rows, M, N);
        if (this->size() == 0)
            allocate(*this, M, N);
        else if ((M != this->rows_) || (N != this->cols_))
            detail::assignment_dimensions(this->rows_, this->cols_, M, N);

        detail::copy_initializer_list(*this, rows);
        return *this;
    }



   /***********************************************************************************************
    *  RowMajorMatrix<T>& operator=(RowMajorMatrix<T2> const&)
    * 
    * Copy assignment for RowMajorMatrix of another element type
    * 
    * Template Parameters:
    *     T2          Data type of other matrix
    * Parameters:
    *     mat         Other RowMajorMatrix
    * 
    ***********************************************************************************************/
    template<class T>
    template<class T2>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T2> const& mat)
    {
        #ifdef MATL_DIMENSION_CHECK
        if (this->data_ != nullptr)
            if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
            {
                detail::assignment_dimensions(this->rows_, this->cols_, mat.rows(), mat.cols());
            }    
        #endif

        #ifdef MATL_DEBUG
        std::cout << "RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T2> const&) called" << std::endl;
        #endif

        // Check for self assignment
        if (static_cast<void const*>(this) != static_cast<void const*>(&mat)) this->copy(mat);

        return *this;
    }



   /***********************************************************************************************
    *  RowMajorMatrix<T>& operator=(MAT const&)
    * 
    * Copy assignment for general Matrix
    * 
    * Template Parameters:
    *     MAT         Other matrix type
    * Parameters:
    *     mat         Other RowMajorMatrix
    * 
    ***********************************************************************************************/
    template<class T>
    template<MatrixType MAT>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (MAT const& mat)
    {
        #ifdef MATL_DIMENSION_CHECK
        if (this->data_ != nullptr)
            if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
            {
                detail::assignment_dimensions(this->rows_, this->cols_, mat.rows(), mat.cols());
            }    
        #endif

        // If no self-assignment, copy data
        if (static_cast<void const*>(this) != static_cast<void const*>(&mat)) this->copy(mat);

        return *this;
    }



   /***********************************************************************************************
    * RowMajorMatrix<T>& operator= (S)
    * 
    * Scalar assignement shortcut   
    * 
    * Template Parameters:
    *    S           Scalar data type
    * Parameters:
    *    val         Scalar value
    ***********************************************************************************************/
    template<class T>
    template<Scalar S>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (S val)
    {
        const index_t N = this->size();
        for (index_t i = 0; i < N; ++i) this->data_[i] = val;
        return *this;
    } 


   /***********************************************************************************************
    * RowMajorMatrix<T>& operator+= (S)
    * 
    * Addition assignment (by scalar)
    * 
    * Parameters:
    *   val         Scalar value
    ***********************************************************************************************/
    template<class T>
    template<Scalar S>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator+= (S val)
    {
        const index_t N = this->size();
        for (index_t i = 0; i < N; ++i) this->data_[i] += val;
        return *this;
    }

   /***********************************************************************************************
    * RowMajorMatrix<T>& operator-= (S)
    * 
    * Subtraction assignment (by scalar)
    * 
    * Parameters:
    *   val         Scalar value
    ***********************************************************************************************/
    template<class T>
    template<Scalar S>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator-= (S val)
    {
        const index_t N = this->size();
        for (index_t i = 0; i < N; ++i) this->data_[i] -= val;
        return *this;
    } 

   /***********************************************************************************************
    * RowMajorMatrix<T>& operator*= (S)
    * 
    * Multiplication assignment (by scalar)
    * 
    * Parameters:
    *   val         Scalar value
    ***********************************************************************************************/
    template<class T>
    template<Scalar S>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator*= (S val)
    {
        const index_t N = this->size();
        for (index_t i = 0; i < N; ++i) this->data_[i] *= val;
        return *this;
    }
    
   /***********************************************************************************************
    * RowMajorMatrix<T>& operator/= (S)
    * 
    * Division assignment (by scalar)
    * 
    * Parameters:
    *   val         Scalar value
    ***********************************************************************************************/
    template<class T>
    template<Scalar S>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator/= (S val)
    {
        const index_t N = this->size();
        if constexpr (std::is_integral<dtype>::value)
        {
            // Integer division cannot be replaced with multiplication by a reciprocal
            for (index_t i = 0; i < N; ++i) this->data_[i] /= val;
        }
        else
        {
            dtype reciprocal = 1.0 / val;
            for (index_t i = 0; i < N; ++i) this->data_[i] *= reciprocal;
        }
            
        return *this;
    } 


   /***********************************************************************************************
    * RowMajorMatrix<T>& operator+= (MAO const&)
    * 
    * Addition assignment (by matrix or expression)
    * 
    * Parameters:
    *     mao               Matrix or expression
    * Returns:
    *     RowMajorMatrix    reference to self
    ***********************************************************************************************/
    template<class T>
    template<MatrixAlgebraObject MAO>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator+= (MAO const& mao)
    {
        *this = *this + mao;
        return *this;
    } 

   /***********************************************************************************************
    * RowMajorMatrix<T>& operator-= (MAO const&)
    *
    * Subtraction assignment (by matrix or expression)
    *
    * Parameters:
    *     mao               Matrix or expression
    * Returns:
    *     RowMajorMatrix    reference to self
    ***********************************************************************************************/
    template<class T>
    template<MatrixAlgebraObject MAO>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator-= (MAO const& mao)
    {
        *this = *this - mao;
        return *this;
    }
   
   /***********************************************************************************************
    * RowMajorMatrix<T>& operator*= (MAO const&)
    * 
    * Multiplication assignment (matrix product, by matrix or expression)
    * 
    * Parameters:
    *     mao               Matrix or expression
    * Returns:
    *     RowMajorMatrix    reference to self
    ***********************************************************************************************/
    template<class T>
    template<MatrixAlgebraObject MAO>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator*= (MAO const& mao)
    {
        *this = *this * mao;
        return *this;
    }
    

    

    
 



   /***********************************************************************************************
    *  void copy(MAT const&)
    * 
    * Copy contents of generic matrix to instance RowMajorMatrix. Allocates new RowMajorMatrix if 
    * it is not already allocated.
    * 
    * Template Parameters:
    *     MAT          Matrix type
    * 
    * Parameters:
    *     mat          Matrix
    * 
    ***********************************************************************************************/ 
    template<class T>
    template<MatrixType MAT>
    void RowMajorMatrix<T>::copy(MAT const& mat)
    {
        const index_t block = 32;//block_size(); // change later
        const index_t M = mat.rows();
        const index_t N = mat.cols();

        // Allocate if not aleady allocated
        if (this->data_ == nullptr) 
            allocate(*this,mat.rows(),mat.cols());
            
        // Dimension check debug flag:
        // Only allow dynamic allocation if the RowMajorMatrix has a size of 0. This is intended
        // to discourage repeated hidden reallocations
        else
        {
            #ifdef MATL_DIMENSION_CHECK
            if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
            {
                detail::assignment_dimensions(this->rows_, this->cols_, mat.rows(), mat.cols());
            }    
            #endif
        }


        // Cache blocking to mitigate miss-rate for intermingled row and column major ordering
        for (index_t ii=1; ii<(M+1); ii+=block)
            for (index_t jj=1; jj<(N+1); jj+=block)
            {
                const index_t istop = std::min(ii+block, M+1);
                const index_t jstop = std::min(jj+block, N+1);
                for (index_t i=ii; i!=istop; ++i)
                    for (index_t j=jj; j!=jstop; ++j)
                        (*this)(i,j) = mat(i,j);
            }

        return;    
    }




   /***********************************************************************************************
    *  void copy(const RowMajorMatrix<T2>&)
    * 
    * Copy contents of other RowMajorMatrix to instance RowMajorMatrix. Allocates new 
    * RowMajorMatrix if it is not already allocated
    * 
    * Template Parameters:
    *     T2          Argument data type
    * 
    * Parameters:
    *     mat         RowMajorMatrix
    * 
    ***********************************************************************************************/ 
    template<class T>
    template<class T2>
    void RowMajorMatrix<T>::copy(RowMajorMatrix<T2> const& mat)
    {
        #ifdef MATL_DEBUG
        std::cout << "void RowMajorMatrix<T>::copy(RowMajorMatrix<T2> const& mat)" << std::endl;
        #endif
        // Allocate if not aleady allocated
        if (this->data_ == nullptr) 
            allocate(*this, mat.rows(), mat.cols());
 
        // Dimension check debug flag:
        // Only allow dynamic allocation if the RowMajorMatrix has a size of 0. This is intended
        // to discourage repeated hidden reallocations
        
        #ifdef MATL_DIMENSION_CHECK
        if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
        {
            detail::assignment_dimensions(this->rows_, this->cols_, mat.rows(), mat.cols());
        }    
        #endif
        
        // Copy data
        const index_t N = this->size();
        for (index_t i=0; i<N; ++i) this->data_[i] = mat.data_[i];

        return;    
    }




   /*=============================================================================================*
    *                                         DATA ACCESS
    *=============================================================================================*/
    
   
   /***********************************************************************************************
    * T operator() (index_t) const
    * 
    * Single index data getter
    * 
    * Parameters:
    *     i       index
    * 
    * Returns;
    *     T       1-indexed data corresponding to i in Row Major order 
    ***********************************************************************************************/ 
    template<class T>
    inline T RowMajorMatrix<T>::operator() (index_t i) const
    {
        // Bounds checking debug flag
        #ifdef MATL_BOUNDS_CHECK
        if ((i<1)||(i>this->size()))
        {
            detail::index_out_of_range(i, this->size());
        }
        #endif

        return this->data_[i-1];
    }


   /***********************************************************************************************
    * T operator() (index_t, index_t) const
    * 
    * Single index data getter
    * 
    * Parameters:
    *     i       row index
    *     j       column index
    * 
    * Returns;
    *     T       1-indexed data corresponding to (i,j) index pair 
    ***********************************************************************************************/ 
    template<class T>
    inline T RowMajorMatrix<T>::operator() (index_t i, index_t j) const
    {
        // Bounds checking debug flag
        #ifdef MATL_BOUNDS_CHECK
        if ( ((i<1)||(i>this->rows())) || ((j<1)||(j>this->cols())) )
        {
            detail::index_out_of_range(i, j, this->rows(), this->cols());
        }
        #endif

        return this->data_[ this->cols_*(i-1) + (j-1) ];
    }



   /***********************************************************************************************
    * T& operator() (index_t)
    * 
    * Single index data accessor.
    * 
    * Parameters:
    *     i       index
    * 
    * Returns;
    *     T       1-indexed data corresponding to i in Row Major order 
    ***********************************************************************************************/ 
    template<class T>
    inline T& RowMajorMatrix<T>::operator() (index_t i)
    {
        // Bounds checking debug flag
        #ifdef MATL_BOUNDS_CHECK
        if ((i<1)||(i>this->size()))
        {
            detail::index_out_of_range(i, this->size());
        }
        #endif
        return this->data_[i-1];
    }



   /***********************************************************************************************
    * T& operator() (index_t, index_t)
    * 
    * Single index data getter
    * 
    * Parameters:
    *     i       row index
    *     j       column index
    * 
    * Returns;
    *     T       1-indexed data corresponding to (i,j) index pair 
    ***********************************************************************************************/ 
    template<class T>
    inline T& RowMajorMatrix<T>::operator() (index_t i, index_t j)
    {
        // Bounds checking debug flag
        #ifdef MATL_BOUNDS_CHECK
        if ( ( (i<1)||(i>this->rows()) ) || ( (j<1)||(j>this->cols()) ) )
        {
            detail::index_out_of_range(i, j, this->rows(), this->cols());
        }
        #endif

        return this->data_[ this->cols_*(i-1) + (j-1) ];
    }







   /*=============================================================================================*
    *                                          UTILITIES
    *=============================================================================================*/





   /***********************************************************************************************
    * index_t rows() const
    * 
    * Number of rows of the matrix.
    * 
    * Parameters:
    * Returns:
    *     index_t      Number of rows
    ***********************************************************************************************/ 
    template<class T>
    inline index_t RowMajorMatrix<T>::rows() const
    {return this->rows_;}

   /***********************************************************************************************
    * index_t cols() const
    * 
    * Number of columns of the matrix.
    * 
    * Parameters:
    * Returns:
    *     index_t      Number of columns
    ***********************************************************************************************/
    template<class T>
    inline index_t RowMajorMatrix<T>::cols() const
    {return this->cols_;}

   /***********************************************************************************************
    * index_t size() const
    * 
    * Number of elements of the matrix.
    * 
    * Parameters:
    * Returns:
    *     index_t      Number of elements
    ***********************************************************************************************/
    template<class T>
    inline index_t RowMajorMatrix<T>::size() const
    {return this->rows_ * this->cols_;}


   /***********************************************************************************************
    * void resize(index_t, index_t)
    * 
    * Resize matrix to specified shape.
    * 
    * Parameters:
    *     rows         New number of rows
    *     cols         New number of columns
    * Returns:
    ***********************************************************************************************/ 
    /*
    template<class T>
    void RowMajorMatrix<T>::resize(index_t rows, index_t cols)
    {
        index_t new_size = rows*cols;
        index_t old_size = this->size();

        if (new_size <= (old_size + this->buffer()) )
        {
            this->rows_ = rows;
            this->cols_ = cols;
            this->buffer_ = old_size + this->buffer() - new_size;

            for (index_t i = old_size; i < new_size; ++i) {this->data_[i] = 0;}
        }
        else
        {
            T* new_data = new T[new_size];
            for (index_t i = 0; i < old_size; ++i) {new_data[i] = this->data_[i];}

            this->rows_ = rows;
            this->cols_ = cols;

            if (!this->external_mem_) delete[] this->data_;

            this->data_ = new_data;  
        }
        return;
    }
    */


   /***********************************************************************************************
    * bool external_mem() const
    * 
    * Returns:
    *    bool      True if memory is externally owned, false otherwise
    ***********************************************************************************************/
    template<class T>
    bool RowMajorMatrix<T>::external_mem() const
    {return this->external_mem_;}
    


   /***********************************************************************************************
    * Prints the RowMajorMatrix in human readable form
    ***********************************************************************************************/
    template<class T>
    void RowMajorMatrix<T>::print() const
    {
        for (index_t i = 1; i <= this->rows(); ++i)
        {
            for (index_t j = 1; j <= this->cols(); ++j)
            {
                std::cout << (*this)(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }


   /***********************************************************************************************
    * Prints the RowMajorMatrix in human readable form with a trailing newline
    ***********************************************************************************************/
    template<class T>
    void RowMajorMatrix<T>::println() const
    {
        for (index_t i = 1; i <= this->rows(); ++i)
        {
            for (index_t j = 1; j <= this->cols(); ++j)
            {
                std::cout << (*this)(i,j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }


   /*=============================================================================================*
    *                                  ALLOCATION AND DEALLOCATION
    *=============================================================================================*/


   /***********************************************************************************************
    * void allocate(RowMajorMatrix<T>&, index_t, index_t)
    * 
    * Allocate memory for a RowMajorMatrix on the heap, owned by the matrix. Any data the matrix
    * held before is released first (see deallocate).
    * 
    * Parameters:
    *     matrix        RowMajorMatrix object to allocate memory for
    *     rows          Number of rows in matrix
    *     cols          Number of columns in matrix
    ***********************************************************************************************/  
    template<class T>
    void allocate(RowMajorMatrix<T>& matrix, index_t rows, index_t cols)
    {
        #ifdef MATL_DEBUG
        std::cout << "In allocate(RowMajorMatrix<T>&, index_t, index_t)..." << std::endl;
        #endif
        deallocate(matrix);
        matrix.data_ = new T[rows*cols];
        matrix.rows_ = rows;
        matrix.cols_ = cols;
    }


   /***********************************************************************************************
    * void deallocate(RowMajorMatrix<T>&)
    * 
    * Deallocate the memory of a RowMajorMatrix and leave it empty (0x0). Externally owned memory
    * is not deleted; the matrix just stops using it.
    * 
    * Parameters:
    *     matrix        RowMajorMatrix object to deallocate memory for
    ***********************************************************************************************/
    template<class T>
    void deallocate(RowMajorMatrix<T>& matrix)
    {
        #ifdef MATL_DEBUG
        std::cout << "In deallocate(RowMajorMatrix<T>&)..." << std::endl;
        #endif
        if (!matrix.external_mem_) delete[] matrix.data_;

        matrix.data_         = nullptr;
        matrix.rows_         = 0;
        matrix.cols_         = 0;
        matrix.external_mem_ = false;
    }
}
