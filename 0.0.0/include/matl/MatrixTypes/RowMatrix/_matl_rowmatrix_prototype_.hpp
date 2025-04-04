//#ifndef _MATL_ROWMATRIX_PROTOTYPE__HPP
//#define _MATL_ROWMATRIX_PROTOTYPE__HPP

//#include <iostream>
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"
//#include "_matl_matrix_.hpp"



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
    class RowMatrix : public Matrix< RowMatrix<T> >
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

        // Blank copy
        RowMatrix<T> blank_copy() const;


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
        template<class T2, enable_if_scalar<T2> = true>
        void operator= (const T2);

        // Row Getter
        index_t rows() const;

        // Coloumn Getter
        index_t cols() const;

        // Size Getter
        index_t size() const;

        // Buffer size
        index_t buffer() const;

        // Resize
        void resize(index_t, index_t);

        // Print
        void print() const;

        // Print with newline
        void println() const;
    };
}

//#endif