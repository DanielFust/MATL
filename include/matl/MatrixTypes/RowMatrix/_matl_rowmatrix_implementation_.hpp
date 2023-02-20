#ifndef _MATL_ROWMATRIX_IMPLEMENTATION__HPP
#define _MATL_ROWMATRIX_IMPLEMENTATION__HPP

//#include <iostream>
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"
//#include "_matl_matrix_.hpp"



namespace matl
{

   /*=============================================================================================*
    *                                      CONSTRUCTORS
    *=============================================================================================*/



   /***********************************************************************************************
    * RowMatrix()
    * 
    * Default constructor 
    * 
    * Template Parameters:
    *     T     Matrix data type
    * 
    * Parameters:
    * 
    * Returns:
    *     RowMatrix      Row major matrix object with zero rows and zero columns
    ***********************************************************************************************/ 
    template<class T>
    RowMatrix<T>::RowMatrix() : rows_(0), cols_(0), buffer_(0), data_(nullptr) {}



   /***********************************************************************************************
    * RowMatrix(index_t rows, index_t cols)
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
    *     RowMatrix      Row major matrix object with specified rows and columns
    ***********************************************************************************************/ 
    template<class T>
    RowMatrix<T>::RowMatrix(index_t rows, index_t cols) : rows_(rows), 
                                                          cols_(cols), 
                                                          buffer_(0)
    {this->data_ = new T[rows*cols];}


   /***********************************************************************************************
    * RowMatrix(index_t rows, index_t cols, index_t buffer)
    * 
    * Sized buffered constructor 
    *     Allocates extra contiguous data storage space
    * 
    * Template Parameters:
    *     T     Matrix data type
    * 
    * Parameters:
    *     rows          number of rows
    *     cols          number of columns
    *     buffer        buffer size
    * 
    * Returns:
    *     RowMatrix      Row major matrix object with specified rows, columns, and data buffer
    ***********************************************************************************************/ 
    template<class T>
    RowMatrix<T>::RowMatrix(index_t rows, index_t cols, index_t buffer) : rows_(rows), 
                                                                          cols_(cols), 
                                                                          buffer_(buffer)
    {this->data_ = new T[rows*cols + buffer];}



   



   /***********************************************************************************************
    * Destructor
    * 
    * Template Parameters:
    *     T        (implicit) Matrix data type
    * 
    * Parameters:
    ***********************************************************************************************/ 
    template<class T>
    RowMatrix<T>::~RowMatrix()
    {
        std::cout << "RowMatrix<T> Destructor called" << std::endl;
        //this->println();
        delete[] this->data_;
    }





    


   /*=============================================================================================*
    *                                   Assignment and Copying
    *=============================================================================================*/
  

  

   /***********************************************************************************************
    * Constructs or repopulates a RowMatrix from an expression template
    ***********************************************************************************************/
    template<class T>
    template<typename EXP_T>
    void RowMatrix<T>::operator= (_MATEXPRESSION_<EXP_T> const& expression)
    {
        std::cout << "Commencing RowMatrix::operator=" << std::endl;
        // Dimension check debug flag:
        // Only allow dynamic allocation if the RowMatrix has a size of 0. This is intended
        // to discourage repeated hidden reallocations 
        #ifdef MATL_DIM_CHECK
        if (  ((this->rows_ != expression.rows()) || (this->cols_ != expression.cols()))
            && (this->size() != 0) ) 
        {
            std::cerr << "Error in 'void RowMatrix<T1>::operator= (const _MATEXPRESSION_<EXP_T> "
                          "&expression)'"
                      << "\n    RowMatrix with dimensions " << this->rows_ << "x" << this->cols_
                      << " incompatible with _MATEXPRESSION_ with dimensions " << expression.rows()
                      << "x" << expression.cols() 
                      << std::endl;
            abort();          
        }    
        #endif
 
        // Matrix is not allocated yet
        if (this->size() == 0) {
            this->resize(expression.rows(), expression.cols());
        }

        // Check for memory conflict
        if (expression.memory_conflict(*this))
        {
            //std::cout << "expression.memory_conflict(*this) = " << expression.memory_conflict(*this) << std::endl;
            RowMatrix<T> mat(this->rows(), this->cols());
            //std::cout<< "mat.print():"<<std::endl;
            //mat.print();
            expression.eval(mat);
            //std::cout<< "mat.print():"<<std::endl;
            //mat.print();
            *this = static_cast<RowMatrix<T>&&>(mat);
        } 
        else {expression.eval(*this);}
                
        std::cout << "Returning from RowMatrix::operator=" << std::endl;
        return;    
    }







   /***********************************************************************************************
    *  void operator=(RowMatrix<T>&&)
    * 
    * Same type RowMatrix rvalue assignment. Takes ownership rvalue's data without copying.
    * 
    * Template Parameters:
    * Parameters:
    *     mat          RowMatrix
    * 
    ***********************************************************************************************/
    template<class T>
    void RowMatrix<T>::operator= (RowMatrix<T>&& mat)
    {
        #ifdef MATL_DIM_CHECK
            if (this->data_ != nullptr)
                if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
                {
                    std::cerr << "Error in 'void RowMatrix<T>::operator= (RowMatrix<T>&&)'"
                            << "\n    RowMatrix with dimensions " << this->rows_ << "x" 
                            << this->cols_
                            << " incompatible with Matrix with dimensions " << mat.rows()
                            << "x" << mat.cols() 
                            << std::endl;
                    abort();          
                }    
        #endif

        // Delete any existing data
        if (this->data_ != nullptr) delete[] this->data_;

        // Assign data pointer to pointer of other matrix data
        this->data_ = mat.data_;
        this->rows_ = mat.rows_;
        this->cols_ = mat.cols_;
        this->buffer_ = mat.buffer_;

        mat.data_ = nullptr;
        mat.rows_ = 0;
        mat.cols_ = 0;
        mat.buffer_ = 0;


        return;
    }






   /***********************************************************************************************
    *  void operator=(RowMatrix&)
    * 
    * Copy assignment for RowMatrix
    * 
    * Template Parameters:
    *     T2          Data type of other matrix
    * Parameters:
    *     mat         Other RowMatrix
    * 
    ***********************************************************************************************/
    template<class T>
    template<class T2>
    void RowMatrix<T>::operator= (RowMatrix<T2>& mat)
    {
        #ifdef MATL_DIM_CHECK
            if (this->data_ != nullptr)
                if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
                {
                    std::cerr << "Error in 'void RowMatrix<T>::operator= (const RowMatrix<T2>&)'"
                            << "\n    RowMatrix with dimensions " << this->rows_ << "x" 
                            << this->cols_
                            << " incompatible with RowMatrix with dimensions " << mat.rows()
                            << "x" << mat.cols() 
                            << std::endl;
                    abort();          
                }    
        #endif

        std::cout << "void RowMatrix<T>::operator= (const RowMatrix<T2>&) called" << std::endl;

        this->copy(mat);

        return;
    }



   /***********************************************************************************************
    *  void operator=(const Matrix&)
    * 
    * Copy assignment for general Matrix
    * 
    * Template Parameters:
    *     MAT_T       Other matrix type
    * Parameters:
    *     mat         Other RowMatrix
    * 
    ***********************************************************************************************/
    template<class T>
    template<typename MAT_T>
    void RowMatrix<T>::operator= (Matrix<MAT_T>& mat)
    {
        #ifdef MATL_DIM_CHECK
            if (this->data_ != nullptr)
                if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
                {
                    std::cerr << "Error in 'void RowMatrix<T>::operator= "
                              << "(const Matrix<MAT_T>& mat)'"
                              << "\n    RowMatrix with dimensions " << this->rows_ << "x" 
                              << this->cols_
                              << " incompatible with Matrix with dimensions " << mat.rows()
                              << "x" << mat.cols() 
                              << std::endl;
                    abort();          
                }    
        #endif

        //std::cout << "void RowMatrix<T>::operator= (const Matrix<MAT_T>&) called" << std::endl;

        this->copy(mat);

        return;
    }





   /***********************************************************************************************
    *  void copy(const Matrix<MAT_T>&)
    * 
    * Copy contents of generic matrix to instance RowMatrix. Allocates new RowMatrix if it is not
    * already allocated
    * 
    * Template Parameters:
    *     MAT_T        Matrix type
    * 
    * Parameters:
    *     mat          Matrix
    * 
    ***********************************************************************************************/ 
    template<class T>
    template<typename MAT_T>
    void RowMatrix<T>::copy(const Matrix<MAT_T>& mat)
    {
        const index_t block = block_size();
        const index_t M = mat.rows();
        const index_t N = mat.cols();

        // Allocate if not aleady allocated
        if (this->data_ == nullptr) 
        {
            this->rows_ = mat.rows();
            this->cols_ = mat.cols();
            this->buffer_ = mat.buffer();

            this->data_ = new T[this->rows() + this->cols() + this->buffer()];
        }
        // Dimension check debug flag:
        // Only allow dynamic allocation if the RowMatrix has a size of 0. This is intended
        // to discourage repeated hidden reallocations
        else
        {
            #ifdef MATL_DIM_CHECK
            if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
            {
                std::cerr << "Error in 'void RowMatrix<T>::copy() (const Matrix<MAT_T>&)'"
                        << "\n    RowMatrix with dimensions " << this->rows_ << "x" << this->cols_
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
                const index_t istop = branchless_min(ii+block, M);
                const index_t jstop = branchless_min(jj+block, N);
                for (index_t i=ii; i!=istop; ++i)
                    for (index_t j=jj; j!=jstop; ++j)
                        (*this)(i,j) = mat(i,j);
            }

        return;    
    }






   /***********************************************************************************************
    *  void copy(const RowMatrix<T2>&)
    * 
    * Copy contents of other RowMatrix to instance RowMatrix. Allocates new RowMatrix if it is not
    * already allocated
    * 
    * Template Parameters:
    *     T2          Argument data type
    * 
    * Parameters:
    *     mat         RowMatrix
    * 
    ***********************************************************************************************/ 
    template<class T>
    template<class T2>
    void RowMatrix<T>::copy(RowMatrix<T2> const& mat)
    {
        // Allocate if not aleady allocated
        if (this->data_ == nullptr) 
        {
            this->rows_ = mat.rows();
            this->cols_ = mat.cols();
            this->buffer_ = mat.buffer();

            this->data_ = new T[this->rows() + this->cols() + this->buffer()];
        }
        // Dimension check debug flag:
        // Only allow dynamic allocation if the RowMatrix has a size of 0. This is intended
        // to discourage repeated hidden reallocations
        
        #ifdef MATL_DIM_CHECK
        if ( (this->rows_ != mat.rows()) || (this->cols_ != mat.cols()) )
        {
            std::cerr << "Error in 'void RowMatrix<T>::copy() (const RowMatrix<T2>&)'"
                    << "\n    RowMatrix with dimensions " << this->rows_ << "x" << this->cols_
                    << " incompatible with Matrix with dimensions " << mat.rows()
                    << "x" << mat.cols() 
                    << std::endl;
            abort();          
        }    
        #endif
        

        const index_t N = this->size() + 1;
        for (index_t i=1; i!=N; ++i) (*this)(i) = mat(i);

        return;    
    }




   /***********************************************************************************************
    *  RowMatrix blank_copy() const
    * 
    * Creates a blank copy of the matrix with the same dimensions and size
    * 
    * Template Parameters:
    *     T            Matrix data type  
    * 
    * Parameters:
    * Returns:
    *     RowMatrix    blank copy of this RowMatrix
    ***********************************************************************************************/ 
    template<class T>
    RowMatrix<T> RowMatrix<T>::blank_copy() const
    {
        RowMatrix<T> mat(this->rows(), this->cols(), this->buffer());
        return mat;
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
    inline T RowMatrix<T>::operator() (index_t i) const
    {
        // Bounds checking debug flag
        #ifdef MATL_BOUNDS_CHECK
        if ((i<1)||(i>this->size()))
        {
            std::cerr << "\nError in 'inline T RowMatrix<T>::operator() (unsigned long i) const' : "
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
    inline T RowMatrix<T>::operator() (index_t i, index_t j) const
    {
        // Bounds checking debug flag
        #ifdef MATL_BOUNDS_CHECK
        if ( ((i<1)||(i>this->rows())) || ((j<1)||(j>this->cols())) )
        {
            std::cerr << "\nError in 'inline T RowMatrix<T>::operator() (unsigned long, unsigned " 
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
    inline T& RowMatrix<T>::operator() (index_t i)
    {
        // Bounds checking debug flag
        #ifdef MATL_BOUNDS_CHECK
        if ((i<1)||(i>this->size()))
        {
            std::cerr << "\nError in 'inline T& RowMatrix<T>::operator() (unsigned long i)' : "
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
    inline T& RowMatrix<T>::operator() (index_t i, index_t j)
    {
        // Bounds checking debug flag
        #ifdef MATL_BOUNDS_CHECK
        if ( ( (i<1)||(i>this->rows()) ) || ( (j<1)||(j>this->cols()) ) )
        {
            std::cerr << "\nError in 'inline T& RowMatrix<T>::operator() (unsigned long, unsigned "
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
    void RowMatrix<T1>::operator= (const T2 val)
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
    inline index_t RowMatrix<T>::rows() const
    {
        return this->rows_;
    }



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
    inline index_t RowMatrix<T>::cols() const
    {
        return this->cols_;
    }



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
    inline index_t RowMatrix<T>::size() const
    {
        return this->rows_ * this->cols_;
    }


   
   /***********************************************************************************************
    * index_t buffer() const
    * 
    * Size of data buffer.
    * 
    * Parameters:
    * Returns:
    *     index_t      Data buffer size
    ***********************************************************************************************/
    template<class T>
    inline index_t RowMatrix<T>::buffer() const
    {
        return this->buffer_;
    }


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
    template<class T>
    void RowMatrix<T>::resize(index_t rows, index_t cols)
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

            delete[] this->data_;
            this->data_ = new_data;  
        }
        return;
    }



   /***********************************************************************************************
    * Prints the RowMatrix in human readable form
    ***********************************************************************************************/
    template<class T>
    void RowMatrix<T>::print() const
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
    * Prints the RowMatrix in human readable form with a trailing newline
    ***********************************************************************************************/
    template<class T>
    void RowMatrix<T>::println() const
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
    * index_t block_size()
    *
    * Computes a cache block size based on current thread's chunk size
    ***********************************************************************************************/ 
    template<class T>
    inline index_t RowMatrix<T>::block_size()
    {return chunk_size()/sizeof(T);}


}



#endif