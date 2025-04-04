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
    * class RowMajorMatrix
    * 
    * Row major matrix class  
    * 
    * Template Arguments:
    *     T     Matrix data type     
    ***********************************************************************************************/ 
    template<class T>
    class RowMajorMatrix : public MATRIX< RowMajorMatrix<T> >
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
        RowMajorMatrix();

        // Sized constructor
        RowMajorMatrix(index_t, index_t);

        // Sized buffered constructor
        RowMajorMatrix(index_t, index_t, index_t);

        // Destructor
        ~RowMajorMatrix();

        /***************** Assignment and Copying **************/

        // Expression execution trigger
        template<typename EXPR_T>
        void operator= (_MATEXPRESSION_<EXPR_T> const&);

        // Same type RowMajorMatrix rvalue assignment
        void operator= (RowMajorMatrix<T>&&);

        // Copy assignment/constructor for RowMajorMatrix
        template<class T2>
        void operator= (RowMajorMatrix<T2>&);

        // Copy assignment/constructor for General Matrix
        template<typename MAT_T>
        void operator= (MATRIX<MAT_T>&);

        // General matrix copy
        template<typename MAT_T>
        void copy(const MATRIX<MAT_T>&);

        // RowMajorMatrix copy
        template<class T2>
        void copy(const RowMajorMatrix<T2>&);

        // Blank copy
        RowMajorMatrix<T> blank_copy() const;


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