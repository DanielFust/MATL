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
    RowMajorMatrix<T>::RowMajorMatrix() : ScratchMemoryBorrower(),
                                          rows_(0), 
                                          cols_(0), 
                                          data_(nullptr)

    {
        this->set_scratch_container(nullptr);
        this->set_external_mem(false);
    }



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
    RowMajorMatrix<T>::RowMajorMatrix(index_t rows, index_t cols) : ScratchMemoryBorrower(),
                                                                    rows_(rows), 
                                                                    cols_(cols),
                                                                    data_(nullptr)
    {
        this->set_external_mem(false);

        // Check if memory should be allocated in scratch memory or on the heap
        MemoryManager& memory_manager = MemoryManager::get_memory_manager();
        if (memory_manager.enabled() && memory_manager.ALLOCATE_MATRICES_IN_SCRATCH)
            memory_manager.assign_memory(*this, rows, cols);
        else
        {
            this->data_ = new T[rows*cols];
            this->set_scratch_container(nullptr);
        }
    }



   /***********************************************************************************************
    * RowMajorMatrix(index_t rows, index_t cols, T* data)
    * 
    * Preallocated constructor. Allows a RowMajorMatrix object to take ownership of preallocated
    * memory. 
    * 
    * This is intended to allows data for many persistant matrices to be allocated (or transferred
    * via MPI call) in a single contiguous block of memory without additional overhead of many
    * ScratchArray containers.
    * 
    * Note: This memory is not owned by the MemoryManager and will not be deallocated by
    * the RowMajorMatrix destructor.
    ***********************************************************************************************/
    template<class T>   
    RowMajorMatrix<T>::RowMajorMatrix(index_t rows, index_t cols, T* data) : ScratchMemoryBorrower(true),
                                                                             rows_(rows), 
                                                                             cols_(cols), 
                                                                             data_(data)
    {
        this->set_external_mem(true);
        this->set_scratch_container(nullptr);
    }

   /***********************************************************************************************
    * RowMajorMatrix(_MATRIX_<MAT> const&)
    * 
    * Copy constructor. Creates a new RowMajorMatrix object with the same dimensions and data as
    * the original matrix. The data is copied to a new array owned by the new matrix.
    * 
    * Parameters: 
    *    other          Matrix object to copy
    **************************************************************************************************/
    template<class T>
    template<typename MAT>
    RowMajorMatrix<T>::RowMajorMatrix(_MATRIX_<MAT> const& other) : ScratchMemoryBorrower(),
                                                                    rows_(other.rows()), 
                                                                    cols_(other.cols()),
                                                                    data_(nullptr)
    {
        const index_t M = this->rows_;
        const index_t N = this->cols_;

        this->set_external_mem(false);

        // Allocate memory for the new matrix
        allocate(*this, M, N);

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
    RowMajorMatrix<T>::RowMajorMatrix(const RowMajorMatrix<T>& other) : ScratchMemoryBorrower(),
                                                                        rows_(other.rows_), 
                                                                        cols_(other.cols_),
                                                                        data_(nullptr)
    {
        this->set_external_mem(false);

        // Allocate memory for the new matrix
        allocate(*this, rows_, cols_);

        // Copy data from the original matrix
        const index_t N = this->rows_ * this->cols_;
        for (index_t i = 0; i < N; ++i)
            this->data_[i] = other.data_[i];
    } 

   /***********************************************************************************************
    * RowMajorMatrix(_EXPRESSION_<EXPR> const&)
    * 
    * Expression constructor. Creates a new RowMajorMatrix sized to the expression and evaluates
    * the expression into it. Memory is allocated according to the MemoryManager settings.
    * 
    * Parameters:
    *     expr           Expression to evaluate
    ***********************************************************************************************/ 
    template<class T>
    template<typename EXPR>
    RowMajorMatrix<T>::RowMajorMatrix(_EXPRESSION_<EXPR> const& expr) : ScratchMemoryBorrower(),
                                                                        rows_(0), 
                                                                        cols_(0),
                                                                        data_(nullptr)
    {*this = expr;}

   /***********************************************************************************************
    * RowMajorMatrix(RowMajorMatrix<T>&&)
    * 
    * Move constructor. Takes ownership of the other matrix's data without copying. If the data is
    * borrowed from the MemoryManager, the ScratchContainer pairing is transferred to this matrix.
    * The other matrix is left empty (0x0).
    * 
    * Parameters:
    *     other          RowMajorMatrix object to move
    ***********************************************************************************************/ 
    template<class T>
    RowMajorMatrix<T>::RowMajorMatrix(RowMajorMatrix<T>&& other) : ScratchMemoryBorrower(other.external_mem()),
                                                                   rows_(other.rows_), 
                                                                   cols_(other.cols_),
                                                                   data_(other.data_)
    {
        // Transfer the ScratchContainer pairing if the memory is managed by the MemoryManager
        if (other.get_scratch_container() != nullptr)
            ScratchMemoryBorrower::operator=(static_cast<ScratchMemoryBorrower&&>(other));

        // Nullify other matrix so the data is not deleted
        other.data_   = nullptr;
        other.rows_   = 0;
        other.cols_   = 0;
        other.scratch_ = nullptr;
        other.set_external_mem(false);
    } 



   /***********************************************************************************************
    * Destructor
    * 
    * Template Parameters:
    *     T        (implicit) Matrix data type
    * 
    * Parameters:
    ***********************************************************************************************/ 
    template<class T>
    RowMajorMatrix<T>::~RowMajorMatrix()
    {deallocate(*this);}
   /*
    // If data is owned by memory manager, free it (do not delete)
        if (this->get_scratch_container() != nullptr) 
            MemoryManager::get_memory_manager().free_memory(*this->get_scratch_container());
        // Else memory is owned by matrix, delete it
        else if (!this->external_mem()) delete[] this->data_;
    }
    */





   /*=============================================================================================*
    *                                   Assignment and Copying
    *=============================================================================================*/
  

  

   /***********************************************************************************************
    * RowMajorMatrix<T>& operator= (_EXPRESSION_ const&)
    * 
    * Triggers the evaluation of an expression and assigns the result to the RowMajorMatrix. If the
    * RowMajorMatrix is not allocated, it will be allocated to the size of the expression, 
    * otherwise it must have the same dimensions as the expression.
    * 
    * Parameters:
    *   expression          Expression to evaluate
    ***********************************************************************************************/
    template<class T>
    template<typename EXPR>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (_EXPRESSION_<EXPR> const& expression)
    {
        #ifdef MATL_DEBUG
        std::cout << "RowMajorMatrix<T>& operator= (_EXPRESSION_ const&)" << std::endl;
        #endif

        // Dimension check debug flag:
        // Only allow dynamic allocation if the RowMajorMatrix has a size of 0. This is intended
        // to discourage repeated hidden reallocations 
        #ifdef MATL_DIMENSION_CHECK
        if (  ((this->rows_ != expression.rows()) || (this->cols_ != expression.cols()))
            && (this->size() != 0) ) 
        {
            std::cerr << "Error in 'void RowMajorMatrix<T>::operator= (const _EXPRESSION_<EXPR> "
                          "&expression)'"
                      << "\n    RowMajorMatrix with dimensions " << this->rows_ << "x" << this->cols_
                      << " incompatible with _EXPRESSION_ with dimensions " << expression.rows()
                      << "x" << expression.cols() 
                      << std::endl;
            abort();          
        }    
        #endif

        // Check (at compile-time) if the evaluation is safe 
        if (expression.unsafe_evaluation_target(*this))
        {
            // Use temporary matrix to store the result
            RowMajorMatrix<T> temp;
            MemoryManager::allocate(temp, expression.rows(), expression.cols());
            expression.eval(temp);

            // Memory should have the same type of ownership after evaluation.
            // If the memory is privately or externally owned, the result must be copied.
            // If it is managed by the MemoryManager, the result can be moved.
            if (this->get_scratch_container()!=nullptr) *this = std::move(temp);    
            else this->copy(temp);    
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


    #if 0
   /***********************************************************************************************
    * RowMajorMatrix<T>& operator= (_MATRIX_TIMES_MATRIX_ const&)
    * 
    * Specialized assignment for matrix-matrix assigment to multiplication expressions. This
    * method provides additional logic to handle the case where the is RowMajorMatrix is also
    * referenced in the expression.
    * 
    * When the RowMajorMatrix is also referenced in the expression, the expression is evaluated
    * to a temporary matrix, then moved or copied to this matrix.
    *
    * Parameters:
    *    expr          Matrix-matrix multiplication expression
    ***********************************************************************************************/
    template<class T>
    template<typename MAO1, typename MAO2>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (_MATRIX_TIMES_MATRIX_<MAO1,MAO2> const& expr)
    {
        std::cout << "RowMajorMatrix<T>& operator= (_MATRIX_TIMES_MATRIX_ const&)" << std::endl;

        // If the RowMajorMatrix is also referenced in the expression, evaluate to a temporary matrix
        if ((expr.mat1_address() == this) || (expr.mat2_address() == this))
        {
            RowMajorMatrix<T> temp;
            if (MemoryManager::ALLOCATE_NEW_MATRICES_IN_SCRATCH)
                MemoryManager::allocate(temp, expr.rows(), expr.cols());
            else
                temp = RowMajorMatrix<T>(expr.rows(), expr.cols());

            // Evaluate to temporary matrix
            expr.eval(temp);
            
            // If this matrix is using external memory, the result must be copied
            if (this->external_mem())
                this->copy(temp);
            // Otherwise, case the temporary matrix as rvalue reference to move data
            else
                *this = std::move(temp);
        }
        // Otherwise, evaluate directly to this matrix
        else
        {
            // If matrix is not allocated, allocate it
            if (this->size() == 0) 
                MemoryManager::allocate(*this, expr.rows(), expr.cols());

            // Evaluate matrix multiplication to this matrix    
            expr.eval(*this);
        }
            
        return *this;
    } 
    #endif




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
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T>&& mat)
    {
        #ifdef MATL_DEBUG
        std::cout << "RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T>&&)" 
                    << std::endl;
        #endif
        #ifdef MATL_DIMENSION_CHECK
        if (this->data_ != nullptr)
            if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
            {
                std::cerr << "Error in 'void RowMajorMatrix<T>::operator= (RowMajorMatrix<T>&&)'"
                        << "\n    RowMajorMatrix with dimensions " << this->rows_ << "x" 
                        << this->cols_
                        << " incompatible with Matrix with dimensions " << mat.rows()
                        << "x" << mat.cols() 
                        << std::endl;
                abort();          
            }    
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
        this->data_    = mat.data_;
        this->rows_    = mat.rows_;
        this->cols_    = mat.cols_;
        this->set_external_mem(mat.external_mem());
        #ifdef MATL_DEBUG
        std::cout << "RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T>&&) - data moved" << std::endl;
        #endif
        
        // Invoke move assignment on the ScratchMemoryBorrower base class if the memory is managed by the MemoryManager
        if (mat.get_scratch_container() != nullptr)
        {
            //this->ScratchMemoryBorrower::operator=(static_cast<ScratchMemoryBorrower&&>(mat));
            ScratchMemoryBorrower::operator=(static_cast<ScratchMemoryBorrower&&>(mat));
        }
        else
            (*this).scratch_ = nullptr;
            //static_cast<ScratchMemoryBorrower&>(*this).scratch_ = nullptr;
        
        // Nullify other matrix so the data is not deleted
        mat.data_      = nullptr;
        mat.rows_      = 0;
        mat.cols_      = 0;
        mat.set_external_mem(false);

        return *this;
    }



   /***********************************************************************************************
    *  RowMajorMatrix<T>& operator=(RowMajorMatrix&)
    * 
    * Copy assignment for RowMajorMatrix
    * 
    * Template Parameters:
    *     T2          Data type of other matrix
    * Parameters:
    *     mat         Other RowMajorMatrix
    * 
    ***********************************************************************************************/
    template<class T>
    template<class T2>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (RowMajorMatrix<T2>& mat)
    {
        #ifdef MATL_DIMENSION_CHECK
        if (this->data_ != nullptr)
            if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
            {
                std::cerr << "Error in 'void RowMajorMatrix<T>::operator= (const RowMajorMatrix<T2>&)'"
                            << "\n    RowMajorMatrix with dimensions " << this->rows_ << "x" 
                            << this->cols_
                            << " incompatible with RowMajorMatrix with dimensions " << mat.rows()
                            << "x" << mat.cols() 
                            << std::endl;
                abort();          
            }    
        #endif

        #ifdef MATL_DEBUG
        std::cout << "RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (const RowMajorMatrix<T2>&) called" << std::endl;
        #endif

        // Check for self assignment
        if (this!=&mat) this->copy(mat);

        return *this;
    }



   /***********************************************************************************************
    *  RowMajorMatrix<T>& operator=(const Matrix&)
    * 
    * Copy assignment for general Matrix
    * 
    * Template Parameters:
    *     MAT_T       Other matrix type
    * Parameters:
    *     mat         Other RowMajorMatrix
    * 
    ***********************************************************************************************/
    template<class T>
    template<typename MAT_T>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (_MATRIX_<MAT_T>& mat)
    {
        #ifdef MATL_DIMENSION_CHECK
        if (this->data_ != nullptr)
            if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
            {
                std::cerr << "Error in 'void RowMajorMatrix<T>::operator= "
                            << "(const Matrix<MAT_T>& mat)'"
                            << "\n    RowMajorMatrix with dimensions " << this->rows_ << "x" 
                            << this->cols_
                            << " incompatible with Matrix with dimensions " << mat.rows()
                            << "x" << mat.cols() 
                            << std::endl;
                abort();          
            }    
        #endif

        //std::cout << "void RowMajorMatrix<T>::operator= (const Matrix<MAT_T>&) called" << std::endl;

        // If no self-assignment, copy data
        if (this != &mat) this->copy(mat);

        return *this;
    }



   /***********************************************************************************************
    * RowMajorMatrix<T>& operator= (const T2)
    * 
    * Scalar assignement shortcut   
    * 
    * Template Parameters:
    *    T2          Scalar data type
    * Parameters:
    *    val         Scalar value
    ***********************************************************************************************/
    template<class T>
    template<class T2, enable_if_scalar<T2>>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator= (const T2 val)
    {
        const index_t N = this->size();
        for (index_t i = 0; i < N; ++i) this->data_[i] = val;
        return *this;
    } 


   /***********************************************************************************************
    * RowMajorMatrix<T>& operator+= (const T2)
    * 
    * Addition assignment (by scalar)
    * 
    * Parameters:
    *   val         Scalar value
    ***********************************************************************************************/
    template<class T>
    template<class T2, enable_if_scalar<T2>>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator+= (const T2 val)
    {
        const index_t N = this->size();
        for (index_t i = 0; i < N; ++i) this->data_[i] += val;
        return *this;
    }

   /***********************************************************************************************
    * RowMajorMatrix<T>& operator-= (const T2)
    * 
    * Subtraction assignment (by scalar)
    * 
    * Parameters:
    *   val         Scalar value
    ***********************************************************************************************/
    template<class T>
    template<class T2, enable_if_scalar<T2>>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator-= (const T2 val)
    {
        const index_t N = this->size();
        for (index_t i = 0; i < N; ++i) this->data_[i] -= val;
        return *this;
    } 

   /***********************************************************************************************
    * RowMajorMatrix<T>& operator*= (const T2)
    * 
    * Multiplication assignment (by scalar)
    * 
    * Parameters:
    *   val         Scalar value
    ***********************************************************************************************/
    template<class T>
    template<class T2, enable_if_scalar<T2>>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator*= (const T2 val)
    {
        const index_t N = this->size();
        for (index_t i = 0; i < N; ++i) this->data_[i] *= val;
        return *this;
    }
    
   /***********************************************************************************************
    * RowMajorMatrix<T>& operator/= (const T2)
    * 
    * Division assignment (by scalar)
    * 
    * Parameters:
    *   val         Scalar value
    ***********************************************************************************************/
    template<class T>
    template<class T2, enable_if_scalar<T2>>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator/= (const T2 val)
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
    * RowMajorMatrix<T>& operator+= (_MATRIX_<MAT_T>&)
    * 
    * Addition assignment (by matrix)
    * 
    * Parameters:
    *     mat               Matrix
    * Returns:
    *     RowMajorMatrix    reference to self
    ***********************************************************************************************/
    template<class T>
    template<typename MAT>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator+= (_MATRIX_<MAT>& mat)
    {
        *this = *this + mat;
        return *this;
    } 

   /***********************************************************************************************
    * RowMajorMatrix<T>& operator-= (_MATRIX_<MAT_T>&)
    *
    * Subtraction assignment (by matrix)
    *
    * Parameters:
    *     mat               Matrix
    * Returns:
    *     RowMajorMatrix    reference to self
    ***********************************************************************************************/
    template<class T>
    template<typename MAT>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator-= (_MATRIX_<MAT>& mat)
    {
        *this = *this - mat;
        return *this;
    }
   
   /***********************************************************************************************
    * RowMajorMatrix<T>& operator*= (_MATRIX_<MAT_T>&)
    * 
    * Multiplication assignment (by matrix)
    * 
    * Parameters:
    *     mat               Matrix
    * Returns:
    *     RowMajorMatrix    reference to self
    ***********************************************************************************************/
    template<class T>
    template<typename MAT>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator*= (_MATRIX_<MAT>& mat)
    {
        *this = *this * mat;
        return *this;
    }
    
   /***********************************************************************************************
    * RowMajorMatrix<T>& operator+= (_EXPRESSION_<EXPR> const&)
    *
    * Addition assignment (by matrix expression)
    *
    * Parameters:
    *     expr              Matrix expression
    * Returns:
    *     RowMajorMatrix    reference to self
    ***********************************************************************************************/
    template<class T>
    template<typename EXPR>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator+= (_EXPRESSION_<EXPR> const& expr)
    {
        *this = *this + expr;
        return *this;
    }
    
   /***********************************************************************************************
    * RowMajorMatrix<T>& operator-= (_EXPRESSION_<EXPR> const&)
    *
    * Subtraction assignment (by matrix expression)
    *
    * Parameters:
    *     expr              Matrix expression
    * Returns:
    *     RowMajorMatrix    reference to self
    ***********************************************************************************************/
    template<class T>
    template<typename EXPR>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator-= (_EXPRESSION_<EXPR> const& expr)
    {
        *this = *this - expr;
        return *this;
    }
    
   /***********************************************************************************************
    * RowMajorMatrix<T>& operator*= (_EXPRESSION_<EXPR> const&)
    *
    * Multiplication assignment (by matrix expression)
    *
    * Parameters:
    *     expr              Matrix expression
    * Returns:
    *     RowMajorMatrix    reference to self
    ***********************************************************************************************/
    template<class T>
    template<typename EXPR>
    RowMajorMatrix<T>& RowMajorMatrix<T>::operator*= (_EXPRESSION_<EXPR> const& expr)
    {
        *this = *this * expr;
        return *this;
    } 



   /***********************************************************************************************
    *  void copy(const _MATRIX_<MAT_T>&)
    * 
    * Copy contents of generic matrix to instance RowMajorMatrix. Allocates new RowMajorMatrix if 
    * it is not already allocated.
    * 
    * Template Parameters:
    *     MAT_T        Matrix type
    * 
    * Parameters:
    *     mat          Matrix
    * 
    ***********************************************************************************************/ 
    template<class T>
    template<typename MAT>
    void RowMajorMatrix<T>::copy(const _MATRIX_<MAT>& mat)
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
                std::cerr << "Error in 'void RowMajorMatrix<T>::copy() (const Matrix<MAT_T>&)'"
                        << "\n    RowMajorMatrix with dimensions " << this->rows_ << "x" << this->cols_
                        << " incompatible with Matrix with dimensions " << mat.rows()
                        << "x" << mat.cols() 
                        << std::endl;
                abort();          
            }    
            #endif
        }


        // Cache blocking to mitigate miss-rate for intermingled row and column major ordering
        for (index_t ii=1; ii<(M+1); ii+=block)
            for (index_t jj=1; jj<(N+1); jj+=block)
            {
                const index_t istop = branchless_min(ii+block, M+1);
                const index_t jstop = branchless_min(jj+block, N+1);
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
            std::cerr << "Error in 'void RowMajorMatrix<T>::copy() (const RowMajorMatrix<T2>&)'"
                    << "\n    RowMajorMatrix with dimensions " << this->rows_ << "x" << this->cols_
                    << " incompatible with Matrix with dimensions " << mat.rows()
                    << "x" << mat.cols() 
                    << std::endl;
            abort();          
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
            std::cerr << "\nError in 'inline T RowMajorMatrix<T>::operator() (unsigned long i) const' : "
                      << "\n    Index i=" << i << " out of bounds of matrix"
                         " with size " << this->size() << std::endl;
            abort();          
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
            std::cerr << "\nError in 'inline T RowMajorMatrix<T>::operator() (unsigned long, unsigned " 
                      << "long) const' : "
                      << "\n    Index pair i=" << i << ", j=" << j << " out of bounds of matrix"
                         " with dimensions " << this->rows() << "x" << this->cols() << std::endl;
            abort();          
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
            std::cerr << "\nError in 'inline T& RowMajorMatrix<T>::operator() (unsigned long i)' : "
                      << "\n    Index i=" << i << " out of bounds of matrix"
                      << " with size " << this->size() << std::endl;
            abort();          
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
            std::cerr << "\nError in 'inline T& RowMajorMatrix<T>::operator() (unsigned long, unsigned "
                         "long)' : "
                      << "\n    Index pair i=" << i << ", j=" << j << " out of bounds of matrix"
                         " with dimensions " << this->rows() << "x" << this->cols() << std::endl;
            abort();          
        }
        #endif

        return this->data_[ this->cols_*(i-1) + (j-1) ];
    }







   /*=============================================================================================*
    *                                          UTILITIES
    *=============================================================================================*/


   /***********************************************************************************************
    * void operator= (const T2)
    * 
    * Scalar data writing
    * 
    * Template Parameters:
    *     T2      (implicit) scalar data type
    * 
    * Parameters:
    *     val     Scalar value of all matrix elements
    * 
    * Returns:
    ***********************************************************************************************/ 
    template<class T1>
    template<class T2, enable_if_scalar<T2> >
    void RowMajorMatrix<T1>::operator= (const T2 val)
    {
        const index_t N = this->size();
        for (index_t i = 0; i < N; ++i)
            this->data_[i] = val;
        return;    
    }


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

            if (this->scratch_ != nullptr) 
                MemoryManager::get_memory_manager().free_memory(*this->scratch_);
            else delete[] this->data_;

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
    {return ScratchMemoryBorrower::external_mem();}

   /***********************************************************************************************
    * void set_external_mem(bool)
    * 
    * Set external memory flag
    * 
    * Parameters:
    *    ext_mem      True if memory is externally owned, false otherwise
    ***********************************************************************************************/
    template<class T>
    void RowMajorMatrix<T>::set_external_mem(bool ext_mem)
    {ScratchMemoryBorrower::set_external_mem(ext_mem);} 
    


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


   /***********************************************************************************************
    * void set_scratch_container(ScratchArray*)
    * 
    * Saves the address of the scratch memory container this object is borrowing from.
    ***********************************************************************************************/
    template<class T>
    void RowMajorMatrix<T>::set_scratch_container(ScratchContainer* scratch)
    {
        // If the container is the same as the current container, return
        if (scratch == this->get_scratch_container()) return;

        // If the matrix is allocated data may need to be copied to new ScratchContainer
        if (this->data_ != nullptr)
        {
            // If the current data pointer does not point to the data in the scratch container, 
            // then copy the data. This should be more robust than checking if the scratch
            // is different, as the order of operations in MemoryManager methods may cause
            // temporary misalignment between data pointers and 
            // ScratchContainer-ScratchMemoryBorrower pairing pointers.
            std::byte* data = reinterpret_cast<std::byte*>(this->data_);
            if (data != scratch->data())
            {
                index_t nbytes = this->size()*sizeof(T);
                if (nbytes>scratch->size())
                {
                    std::cerr << "Error in 'void RowMajorMatrix<T>::set_scratch_container(ScratchContainer*)' : "
                              << "\n    ScratchContainer is too small to hold matrix data."
                              << "\n    Deallocating matrix prior to assigning scratch memory is strongly recommended."
                              << std::endl;
                    abort();
                }

                // Copy data to new ScratchContainer
                for (index_t i = 0; i < nbytes; ++i) scratch->data()[i] = data[i];

                // If the memory was previously owned by the MemoryManager, free it
                if (this->get_scratch_container() != nullptr) 
                    MemoryManager::get_memory_manager().free_memory(*this->get_scratch_container());
                // If the memory was previously owned by the matrix, deallocate it
                else if (this->external_mem() == false) delete[] this->data_;

                // External memory flag will always be false after this operation
                this->set_external_mem(false);

                // Set new data pointer
                this->data_ = reinterpret_cast<T*>(scratch->data());
            }
        };
        // Finally, set the scratch container pointer
        ScratchMemoryBorrower::set_scratch_container(scratch);
    }

   /***********************************************************************************************
    * ScratchContainer* get_scratch_container()
    * 
    * Returns the address of the scratch memory container this object is borrowing from.
    ***********************************************************************************************/
    template<class T>
    ScratchContainer* RowMajorMatrix<T>::get_scratch_container()
    {return ScratchMemoryBorrower::get_scratch_container();}  

   /***********************************************************************************************
    * void copy_memory_to_private()
    * 
    * Copies the data from shared scratch memory to newly allocated array owned by the 
    * RowMajorMatrix.
    ***********************************************************************************************/ 
    template<class T>
    void RowMajorMatrix<T>::copy_memory_to_private()
    {
        #ifdef MATL_DEBUG
        std::cout << "In RowMajorMatrix<T>::copy_memory_to_private()..." << std::endl;
        #endif
        // Return if memory is not owned by memory manager
        ScratchContainer* scratch = this->get_scratch_container();
        if (scratch == nullptr) return;

        // Save dimensions (free_memory nullifies the borrower, resetting them to 0)
        const index_t rows = this->rows_;
        const index_t cols = this->cols_;
        const index_t N    = rows*cols;

        // Allocate new data array and copy data from shared memory to new memory
        T* new_data = new T[N];
        for (index_t i = 0; i < N; ++i) {new_data[i] = this->data_[i];}

        // Free shared memory (unpairs this matrix from the container)
        MemoryManager::get_memory_manager().free_memory(*scratch);

        // Assign privately owned data array to object
        this->data_ = new_data;
        this->rows_ = rows;
        this->cols_ = cols;
        this->scratch_ = nullptr;
        this->set_external_mem(false);
    }

   /***********************************************************************************************
    * void nullify_borrowed_memory()
    * 
    * Sets the matrix size to zero and nullifies the data pointers if the memory is borrowed from
    * scratch owned by the MemoryManager 
    ***********************************************************************************************/
    template<class T>
    void RowMajorMatrix<T>::nullify_borrowed_memory()
    {
        #ifdef MATL_DEBUG
        std::cout << "In RowMajorMatrix<T>::nullify_borrowed_memory()..." << std::endl;
        #endif
        // If the memory is not regulated by the MemoryManager, return
        if (ScratchMemoryBorrower::external_mem()) return;
        if (ScratchMemoryBorrower::get_scratch_container() == nullptr) return;

        // Set size to zero and nullify pointers
        this->rows_ = 0;
        this->cols_ = 0;
        this->data_ = nullptr;
        ScratchMemoryBorrower::nullify_borrowed_memory();
    } 

   /***********************************************************************************************
    * std::byte* data() const
    * 
    * Returns: 
    *     Address of the first element of the matrix data
    ***********************************************************************************************/
    template<class T>
    std::byte* RowMajorMatrix<T>::data() const
    {return reinterpret_cast<std::byte*>(this->data_);}  

   /***********************************************************************************************
    * index_t bytes_used() const
    * 
    * Returns:
    *     The number of bytes used by the matrix
    ***********************************************************************************************/ 
    template<class T>
    index_t RowMajorMatrix<T>::bytes_used() const
    {return this->size()*sizeof(T);}
}
