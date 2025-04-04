#ifndef _MATL_ROWMATRIX__HPP
#define _MATL_ROWMATRIX__HPP

#include <iostream>
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"
#include "_matl_utils_.hpp"
#include "_matl_matrix_.hpp"



namespace matl
{
   /***********************************************************************************************
    * class RowMatrix
    * 
    * Row major matrix class  
    * 
    * Template Arguments:
    *     T     Matrix data type     
    ***********************************************************************************************/ 
    template<class T>
    class RowMatrix : public Matrix< RowMatrix<T> >//_ROWMAJOR_<RowMatrix<T>, T>   //: public _MATEXPRESSION_<RowMatrix<T>, T>
    {
      private:
        index_t rows_;
        index_t cols_;
        index_t buffer_;
        T*      data_;

        index_t block_size();

      public:
        
        typedef T dtype;

        /********************* Constructors *********************/
        
        // Default constructor
        RowMatrix();

        // Sized constructor
        RowMatrix(index_t, index_t);

        // Sized buffered constructor
        RowMatrix(index_t, index_t, index_t);

        // Destructor
        ~RowMatrix();

        /***************** Assignment and Copying **************/

        // Expression execution trigger
        template<typename EXPR_T>
        void operator= (_MATEXPRESSION_<EXPR_T> const&);

        // Same type RowMatrix rvalue assignment
        void operator= (RowMatrix<T>&&);

        // Copy assignment/constructor for RowMatrix
        template<class T2>
        void operator= (RowMatrix<T2>&);

        // Copy assignment/constructor for General Matrix
        template<typename MAT_T>
        void operator= (Matrix<MAT_T>&);

        // General matrix copy
        template<typename MAT_T>
        void copy(const Matrix<MAT_T>&);

        // RowMatrix copy
        template<class T2>
        void copy(const RowMatrix<T2>&);


        /********************* Data Access *********************/

        // Single index getter
        dtype operator() (index_t) const;

        // Double index getter
        dtype operator() (index_t, index_t) const;

        // Single index accessor
        dtype& operator() (index_t);

        // Double index accessor
        dtype& operator() (index_t, index_t);

        /********************** Utilities **********************/

        // Data writing shortcut
        void operator= (const T);

        // Row Getter
        index_t rows() const;

        // Coloumn Getter
        index_t cols() const;

        // Size Getter
        index_t size() const;

        // Buffer size
        index_t buffer() const;

        // Print
        void print() const;

        // Print with newline
        void println() const;
    };











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
    {
        this->data_ = new T[rows*cols];
    }


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
    {delete[] this->data_;}





    


   /*=============================================================================================*
    *                                   Assignment and Copying
    *=============================================================================================*/
  

    #if 0
   /***********************************************************************************************
    * Constructs or repopulates a RowMatrix from an expression template
    ***********************************************************************************************/
    template<class T>
    template<typename EXP_T>
    void RowMatrix<T>::operator= (_MATEXPRESSION_<EXP_T> &expression)
    {
        std::cout << "void RowMatrix<T>::operator= (const _MATEXPRESSION_<EXP_T> &) called" << std::endl;

        // Dimension check debug flag:
        // Only allow dynamic allocation if the RowMatrix has a size of 0. This is intended
        // to discourage repeated hidden reallocations 
        #ifdef DIM_CHECK
        if (  ((this->rows_ != expression.rows()) || (this->rows_ != 0))
            ||((this->cols_ != expression.cols()) || (this->cols_ != 0)) )
        {
            std::cerr << "Error in 'void RowMatrix<T1>::operator= (const _MATEXPRESSION_<E,T2> "
                          "&expression)'"
                      << "\n    RowMatrix with dimensions " << this->rows_ << "x" << this->cols_
                      << " incompatible with _MATEXPRESSION_ with dimensions " << expression.rows()
                      << "x" << expression.cols() 
                      << std::endl;
            abort();          
        }    
        #endif
        // Case 1: RowMatrix has same shape as expression. No allocation required
        if (this->size() == expression.size())
        {
            // Fallthrough to assignment
        }
        // Case 2: RowMatrix has dimensions 0x0 (other cases should be filtered with debug flags)
        else if (this->data_ == nullptr)
        {
            this->rows_ = expression.rows();
            this->cols_ = expression.cols();
            this->data_ = new T[this->size() + this->buffer_];
        }
        // Case 3: Buffered size is enough to accomodate expression size
        // (THIS CASE SHOULD BE FILTERED OUT BY DEBUG FLAGS)
        else if ((this->size()+this->buffer_) > expression.size())
        {
            this->buffer_ = this->size()+this->buffer_ - expression.size();
            this->rows_ = expression.rows();
            this->cols_ = expression.cols();
            this->data_ = new T[expression.size()];
        }
        // Case 4: Buffered size is not enough to accomodate expression size
        // (THIS CASE SHOULD BE FILTERED OUT BY DEBUG FLAGS)
        else 
        {
            delete[] this->data_;
            this->rows_ = expression.rows();
            this->cols_ = expression.cols();
            this->data_ = new T[expression.size()];
        }
            
        static_cast<EXP_T&>(expression).eval(*this);
        return;    
    }
    #endif





   /***********************************************************************************************
    * Constructs or repopulates a RowMatrix from an expression template
    ***********************************************************************************************/
    template<class T>
    template<typename EXP_T>
    void RowMatrix<T>::operator= (_MATEXPRESSION_<EXP_T> const& expression)
    {
        // Dimension check debug flag:
        // Only allow dynamic allocation if the RowMatrix has a size of 0. This is intended
        // to discourage repeated hidden reallocations 
        #ifdef DIM_CHECK
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
        // Case 1: RowMatrix has same shape as expression. No allocation required
        if (this->size() == expression.size())
        {
            // Fallthrough to assignment
        }
        // Case 2: RowMatrix has dimensions 0x0 (other cases should be filtered with debug flags)
        else if (this->data_ == nullptr)
        {
            this->rows_ = expression.rows();
            this->cols_ = expression.cols();
            this->data_ = new T[this->size() + this->buffer()];
        }
        // Case 3: Buffered size is enough to accomodate expression size
        // (THIS CASE SHOULD BE FILTERED OUT BY DEBUG FLAGS)
        else if ((this->size()+this->buffer()) > expression.size())
        {
            this->buffer_ = this->size()+this->buffer_ - expression.size();
            this->rows_ = expression.rows();
            this->cols_ = expression.cols();
            this->data_ = new T[expression.size()];
        }
        // Case 4: Buffered size is not enough to accomodate expression size
        // (THIS CASE SHOULD BE FILTERED OUT BY DEBUG FLAGS)
        else 
        {
            delete[] this->data_;
            this->rows_ = expression.rows();
            this->cols_ = expression.cols();
            this->data_ = new T[expression.size()];
        }
                
        expression.eval(*this);
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
        #ifdef DIM_CHECK
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

        std::cout << "void RowMatrix<T>::operator= (RowMatrix<T>&&) called" << std::endl;

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
        #ifdef DIM_CHECK
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
        #ifdef DIM_CHECK
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

        std::cout << "void RowMatrix<T>::operator= (const Matrix<MAT_T>&) called" << std::endl;

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
            #ifdef DIM_CHECK
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
    void RowMatrix<T>::copy(const RowMatrix<T2>& mat)
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
        else
        {
            #ifdef DIM_CHECK
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
        }


        // Cache blocking to mitigate miss-rate for intermingled row and column major ordering
        const index_t N = this->size() + 1;
        for (index_t i=1; i!=N; ++i) (*this)(i) = mat(i);

        return;    
    }











   /*=============================================================================================*
    *                                         DATA ACCESS
    *=============================================================================================*/
    
   
   /***********************************************************************************************
    * Single index data getter
    ***********************************************************************************************/ 
    template<class T>
    inline T RowMatrix<T>::operator() (index_t i) const
    {
        // Bounds checking debug flag
        #ifdef BOUNDS_CHECK
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
    * Double index data getter
    ***********************************************************************************************/
    template<class T>
    inline T RowMatrix<T>::operator() (index_t i, index_t j) const
    {
        // Bounds checking debug flag
        #ifdef BOUNDS_CHECK
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
    * Single index data accessor
    ***********************************************************************************************/
    template<class T>
    inline T& RowMatrix<T>::operator() (index_t i)
    {
        // Bounds checking debug flag
        #ifdef BOUNDS_CHECK
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
    * Double index data accessor
    ***********************************************************************************************/
    template<class T>
    inline T& RowMatrix<T>::operator() (index_t i, index_t j)
    {
        // Bounds checking debug flag
        #ifdef BOUNDS_CHECK
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
    * Overwites all RowMatrix elements with the assigned value
    ***********************************************************************************************/    
    template<class T>
    void RowMatrix<T>::operator= (const T val)
    {
        const index_t N = this->size();
        for (index_t i = 0; i < N; ++i)
            this->data_[i] = val;
        return;    
    }



   /***********************************************************************************************
    * Returns number of rows in the RowMatrix
    ***********************************************************************************************/
    template<class T>
    inline index_t RowMatrix<T>::rows() const
    {
        return this->rows_;
    }



   /***********************************************************************************************
    * Returns number of columns in the RowMatrix
    ***********************************************************************************************/
    template<class T>
    inline index_t RowMatrix<T>::cols() const
    {
        return this->cols_;
    }



   /***********************************************************************************************
    * Returns number of elements in the RowMatrix
    ***********************************************************************************************/
    template<class T>
    inline index_t RowMatrix<T>::size() const
    {
        return this->rows_ * this->cols_;
    }


   
   /***********************************************************************************************
    * Returns the size of the buffer
    ***********************************************************************************************/
    template<class T>
    inline index_t RowMatrix<T>::buffer() const
    {
        return this->buffer_;
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