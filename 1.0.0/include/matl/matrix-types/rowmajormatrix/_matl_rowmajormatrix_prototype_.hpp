namespace matl
{
   /***********************************************************************************************
    * class RowMajorMatrix
    * 
    * Row major matrix class. The matrix owns its data (allocated on the heap), unless it was
    * constructed over externally owned memory, which it uses but never deletes.
    * 
    * Template Arguments:
    *     T     Matrix data type     
    ***********************************************************************************************/ 
    template<class T>
    class RowMajorMatrix : public _MATRIX_< RowMajorMatrix<T> >
    {
      private:
        index_t rows_;
        index_t cols_;
        T*      data_;

        // Is the memory externally owned (used, but not deleted, by the matrix)
        bool    external_mem_;

      public:
        
        typedef T dtype;

        /********************* Constructors *********************/
        
        // Default constructor
        RowMajorMatrix();

        // Sized constructor
        RowMajorMatrix(index_t, index_t);

        // Preallocated constructor
        RowMajorMatrix(index_t, index_t, T*);

        // Copy constructor (General matrix)
        template<typename MAT>
        RowMajorMatrix(_MATRIX_<MAT> const&);

        // Copy constructor (Same type)
        RowMajorMatrix(RowMajorMatrix<T> const&);

        // Expression constructor (evaluates the expression; implicit, so that
        // 'Matrix<double> m = A + B;' evaluates to a new matrix)
        template<typename EXPR>
        RowMajorMatrix(_EXPRESSION_<EXPR> const&);

        // Move constructor
        RowMajorMatrix(RowMajorMatrix<T>&&) noexcept;

        // Destructor
        ~RowMajorMatrix();

        /************* Assignment and Copying ******************/

        // Expression execution trigger
        template<typename EXPR>
        RowMajorMatrix<T>& operator= (_EXPRESSION_<EXPR> const&);


        // RowMajorMatrix move assignment
        RowMajorMatrix<T>& operator= (RowMajorMatrix<T>&&) noexcept;

        // RowMajorMatrix copy assignment (same type)
        RowMajorMatrix<T>& operator= (RowMajorMatrix<T> const&);

        // Copy assignment/constructor for RowMajorMatrix
        template<class T2>
        RowMajorMatrix<T>& operator= (RowMajorMatrix<T2>&);

        // Copy assignment/constructor for General Matrix
        template<typename MAT>
        RowMajorMatrix<T>& operator= (_MATRIX_<MAT>&);

        // Scalar assignement shortcut
        template<class T2, enable_if_scalar<T2> >
        RowMajorMatrix<T>& operator= (const T2);

        // Addition assignment (by scalar)
        template<class T2, enable_if_scalar<T2> = true >
        RowMajorMatrix<T>& operator+= (const T2);

        // Subtraction assignment (by scalar)
        template<class T2, enable_if_scalar<T2> = true >
        RowMajorMatrix<T>& operator-= (const T2);

        // Multiplication assignment (by scalar)
        template<class T2, enable_if_scalar<T2> = true >
        RowMajorMatrix<T>& operator*= (const T2);

        // Division assignment (by scalar)
        template<class T2, enable_if_scalar<T2> = true >
        RowMajorMatrix<T>& operator/= (const T2);

        // Addition assignment (by matrix)
        template<typename MAT>
        RowMajorMatrix<T>& operator+= (_MATRIX_<MAT>&);

        // Subtraction assignment (by matrix)
        template<typename MAT>
        RowMajorMatrix<T>& operator-= (_MATRIX_<MAT>&);

        // Multiplication assignment (by matrix)
        template<typename MAT>
        RowMajorMatrix<T>& operator*= (_MATRIX_<MAT>&);

        // Addition assignment (by matrix expression)
        template<typename EXPR>
        RowMajorMatrix<T>& operator+= (_EXPRESSION_<EXPR> const&);

        // Subtraction assignment (by matrix expression)
        template<typename EXPR>
        RowMajorMatrix<T>& operator-= (_EXPRESSION_<EXPR> const&);

        // Multiplication assignment (by matrix expression)
        template<typename EXPR>
        RowMajorMatrix<T>& operator*= (_EXPRESSION_<EXPR> const&);

        // General matrix copy
        template<typename MAT>
        void copy(const _MATRIX_<MAT>&);

        // RowMajorMatrix copy
        template<class T2>
        void copy(const RowMajorMatrix<T2>&);

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

        // Resize (This might be antithecal to the design)
        //void resize(index_t, index_t);

        // Print
        void print() const;

        // Print with newline
        void println() const;

        // Is the memory externally owned (not deleted by the matrix)
        bool external_mem() const;

        // Allocation and deallocation of the matrix data
        template<class U> friend void allocate(RowMajorMatrix<U>&, index_t, index_t);
        template<class U> friend void deallocate(RowMajorMatrix<U>&);
    };


   /***********************************************************************************************
    * Allocating and Deallocating Matrices
    * 
    * 'allocate' and 'deallocate' are the general purpose procedures for allocating and
    * deallocating matrices (loosely resembling the Fortran90 allocate/deallocate vernacular).
    * Matrix memory is allocated on the heap and owned by the matrix.
    ***********************************************************************************************/

    // Matrix allocation
    template<class T>
    void allocate(RowMajorMatrix<T>&, index_t, index_t);

    // Matrix deallocation
    template<class T>
    void deallocate(RowMajorMatrix<T>&);
}
