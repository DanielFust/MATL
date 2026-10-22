namespace matl
{
    namespace detail
    {
        // Traits of RowMajorMatrix (see _matl_type_traits_.hpp)
        template<class T>
        struct traits< RowMajorMatrix<T> >
        {
            typedef T dtype;
            static constexpr bool         fixed_size    = false;
            static constexpr index_t      rows          = 0;
            static constexpr index_t      cols          = 0;
            static constexpr StorageOrder storage_order = StorageOrder::RowMajor;
        };
    }

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
    class RowMajorMatrix : public MatrixBase< RowMajorMatrix<T> >
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
        template<MatrixType MAT>
        RowMajorMatrix(MAT const&);

        // Copy constructor (Same type)
        RowMajorMatrix(RowMajorMatrix<T> const&);

        // Expression constructor (evaluates the expression; implicit, so that
        // 'Matrix<double> m = A + B;' evaluates to a new matrix)
        template<ExpressionType EXPR>
        RowMajorMatrix(EXPR const&);

        // Move constructor
        RowMajorMatrix(RowMajorMatrix<T>&&) noexcept;

        // Destructor
        ~RowMajorMatrix();

        /************* Assignment and Copying ******************/

        // Expression execution trigger
        template<ExpressionType EXPR>
        RowMajorMatrix<T>& operator= (EXPR const&);


        // RowMajorMatrix move assignment
        RowMajorMatrix<T>& operator= (RowMajorMatrix<T>&&) noexcept;

        // RowMajorMatrix copy assignment (same type)
        RowMajorMatrix<T>& operator= (RowMajorMatrix<T> const&);

        // Copy assignment for RowMajorMatrix of another element type
        template<class T2>
        RowMajorMatrix<T>& operator= (RowMajorMatrix<T2> const&);

        // Copy assignment for General Matrix
        template<MatrixType MAT>
        RowMajorMatrix<T>& operator= (MAT const&);

        // Scalar assignment (sets every element)
        template<Scalar S>
        RowMajorMatrix<T>& operator= (S);

        // Addition assignment (by scalar)
        template<Scalar S>
        RowMajorMatrix<T>& operator+= (S);

        // Subtraction assignment (by scalar)
        template<Scalar S>
        RowMajorMatrix<T>& operator-= (S);

        // Multiplication assignment (by scalar)
        template<Scalar S>
        RowMajorMatrix<T>& operator*= (S);

        // Division assignment (by scalar)
        template<Scalar S>
        RowMajorMatrix<T>& operator/= (S);

        // Addition assignment (by matrix or expression)
        template<MatrixAlgebraObject MAO>
        RowMajorMatrix<T>& operator+= (MAO const&);

        // Subtraction assignment (by matrix or expression)
        template<MatrixAlgebraObject MAO>
        RowMajorMatrix<T>& operator-= (MAO const&);

        // Multiplication assignment (matrix product, by matrix or expression)
        template<MatrixAlgebraObject MAO>
        RowMajorMatrix<T>& operator*= (MAO const&);

        // General matrix copy
        template<MatrixType MAT>
        void copy(MAT const&);

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

        // RowMajorMatrix types of other element types (to copy their data)
        template<class U> friend class RowMajorMatrix;

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
