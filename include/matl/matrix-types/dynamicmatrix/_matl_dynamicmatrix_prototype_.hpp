#ifndef MATL_DYNAMICMATRIX_PROTOTYPE_HPP
#define MATL_DYNAMICMATRIX_PROTOTYPE_HPP
namespace matl
{
    namespace detail
    {
        // Traits of DynamicMatrix (see _matl_type_traits_.hpp)
        template<class T, StorageOrder ORDER>
        struct traits< DynamicMatrix<T, ORDER> >
        {
            typedef T dtype;
            static constexpr bool         fixed_size    = false;
            static constexpr index_t      rows          = 0;
            static constexpr index_t      cols          = 0;
            static constexpr StorageOrder storage_order = ORDER;
        };
    }

   /***********************************************************************************************
    * class DynamicMatrix<T, ORDER>
    * 
    * Matrix whose dimensions are set at run time. The matrix owns its data (allocated on the heap),
    * unless it was constructed over externally owned memory, which it uses but never deletes.
    * Usually named through the aliases RowMajorMatrix<T> (= Matrix<T>) and ColMajorMatrix<T>.
    * 
    * The elements are stored row by row (StorageOrder::RowMajor) or column by column
    * (StorageOrder::ColMajor). Element access is the same for both: A(i,j) is row i, column j.
    * The single index A(i) runs through the elements in storage order.
    * 
    * Template Arguments:
    *     T         Matrix data type     
    *     ORDER     Storage order (RowMajor or ColMajor)
    ***********************************************************************************************/ 
    template<class T, StorageOrder ORDER>
    class DynamicMatrix : public MatrixBase< DynamicMatrix<T, ORDER> >
    {
        static_assert(ORDER == StorageOrder::RowMajor || ORDER == StorageOrder::ColMajor,
                      "DynamicMatrix: the storage order must be RowMajor or ColMajor");

      private:
        index_t rows_;
        index_t cols_;
        T*      data_;

        // Is the memory externally owned (used, but not deleted, by the matrix)
        bool    external_mem_;

        // Position of element (i,j) (1-based) in the data array
        index_t offset(index_t i, index_t j) const
        {
            if constexpr (ORDER == StorageOrder::RowMajor) return this->cols_ * (i - 1) + (j - 1);
            else                                           return this->rows_ * (j - 1) + (i - 1);
        }

        // Elements from one row to the next, and from one column to the next (for views)
        std::ptrdiff_t row_stride() const
        {return (ORDER == StorageOrder::RowMajor) ? static_cast<std::ptrdiff_t>(this->cols_) : 1;}

        std::ptrdiff_t col_stride() const
        {return (ORDER == StorageOrder::RowMajor) ? 1 : static_cast<std::ptrdiff_t>(this->rows_);}

      public:
        
        typedef T dtype;

        /********************* Constructors *********************/
        
        // Default constructor
        DynamicMatrix();

        // Sized constructor
        DynamicMatrix(index_t, index_t);

        // Preallocated constructor (rows*cols elements in the matrix's storage order)
        DynamicMatrix(index_t, index_t, T*);

        // Initializer-list constructor, one inner list per row: Matrix<double> A = {{1, 2}, {3, 4}};
        DynamicMatrix(matrix_rows<T>);

        // Copy constructor (General matrix)
        template<MatrixType MAT>
        DynamicMatrix(MAT const&);

        // Copy constructor (Same type)
        DynamicMatrix(DynamicMatrix const&);

        // Expression constructor (evaluates the expression; implicit, so that
        // 'Matrix<double> m = A + B;' evaluates to a new matrix)
        template<ExpressionType EXPR>
        DynamicMatrix(EXPR const&);

        // Move constructor
        DynamicMatrix(DynamicMatrix&&) noexcept;

        // Destructor
        ~DynamicMatrix();

        /************* Assignment and Copying ******************/

        // Expression execution trigger
        template<ExpressionType EXPR>
        DynamicMatrix& operator= (EXPR const&);


        // Move assignment
        DynamicMatrix& operator= (DynamicMatrix&&) noexcept;

        // Copy assignment (same type)
        DynamicMatrix& operator= (DynamicMatrix const&);

        // Initializer-list assignment (A = {} empties the matrix)
        DynamicMatrix& operator= (matrix_rows<T>);

        // Copy assignment for a DynamicMatrix of another element type (same storage order)
        template<class T2>
        DynamicMatrix& operator= (DynamicMatrix<T2, ORDER> const&);

        // Copy assignment for General Matrix
        template<MatrixType MAT>
        DynamicMatrix& operator= (MAT const&);

        // Scalar assignment (sets every element)
        template<Scalar S>
        DynamicMatrix& operator= (S);

        // Addition assignment (by scalar)
        template<Scalar S>
        DynamicMatrix& operator+= (S);

        // Subtraction assignment (by scalar)
        template<Scalar S>
        DynamicMatrix& operator-= (S);

        // Multiplication assignment (by scalar)
        template<Scalar S>
        DynamicMatrix& operator*= (S);

        // Division assignment (by scalar)
        template<Scalar S>
        DynamicMatrix& operator/= (S);

        // Addition assignment (by matrix or expression)
        template<MatrixAlgebraObject MAO>
        DynamicMatrix& operator+= (MAO const&);

        // Subtraction assignment (by matrix or expression)
        template<MatrixAlgebraObject MAO>
        DynamicMatrix& operator-= (MAO const&);

        // Multiplication assignment (matrix product, by matrix or expression)
        template<MatrixAlgebraObject MAO>
        DynamicMatrix& operator*= (MAO const&);

        // General matrix copy
        template<MatrixType MAT>
        void copy(MAT const&);

        // Copy of a DynamicMatrix with the same storage order (element by element in storage order)
        template<class T2>
        void copy(DynamicMatrix<T2, ORDER> const&);

        /********************* Data Access *********************/

        // Single index getter
        dtype operator() (index_t) const;

        // Double index getter
        dtype operator() (index_t, index_t) const;

        // Single index accessor
        dtype& operator() (index_t);

        // Double index accessor
        dtype& operator() (index_t, index_t);

        // Element access with last: A(last, 1), A(2, last - 1) (ref-qualified like the views
        // below, so that the overloads are distinguished on every compiler)
        template<typename R, typename C> requires detail::ElementIndices<R, C>
        dtype& operator() (R const& r, C const& c) &
        {return (*this)(detail::element_index(r, this->rows_), detail::element_index(c, this->cols_));}

        template<typename R, typename C> requires detail::ElementIndices<R, C>
        dtype operator() (R const& r, C const& c) const&
        {return (*this)(detail::element_index(r, this->rows_), detail::element_index(c, this->cols_));}

        // Views of rows and columns (see MatrixView): A(range(1, 3), all), A(all, 2), ...
        template<typename R, typename C> requires detail::ViewIndices<R, C>
        MatrixView<T, ORDER> operator() (R const& r, C const& c) &
        {return detail::make_view<ORDER>(this->data_, this->rows_, this->cols_, this->row_stride(), this->col_stride(), r, c);}

        // Read-only view of a const matrix
        template<typename R, typename C> requires detail::ViewIndices<R, C>
        MatrixView<const T, ORDER> operator() (R const& r, C const& c) const&
        {
            return detail::make_view<ORDER>(static_cast<const T*>(this->data_), this->rows_, this->cols_,
                                            this->row_stride(), this->col_stride(), r, c);
        }

        // Read-only view of a temporary matrix, which the view keeps alive
        template<typename R, typename C> requires detail::ViewIndices<R, C>
        MatrixView<const T, ORDER> operator() (R const& r, C const& c) &&
        {
            const auto owner = std::make_shared<const DynamicMatrix>(std::move(*this));
            return detail::make_view<ORDER>(owner->data(), owner->rows(), owner->cols(),
                                            owner->row_stride(), owner->col_stride(), r, c, owner);
        }

        // Elements in storage order (e.g. to pass the matrix to a C or Fortran library)
        T*       data()       {return this->data_;}
        T const* data() const {return this->data_;}

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

        // DynamicMatrix types of other element types and storage orders (to copy their data)
        template<class U, StorageOrder O> friend class DynamicMatrix;

        // Allocation and deallocation of the matrix data
        template<class U, StorageOrder O> friend void allocate(DynamicMatrix<U, O>&, index_t, index_t);
        template<class U, StorageOrder O> friend void deallocate(DynamicMatrix<U, O>&);
    };


   /***********************************************************************************************
    * Allocating and Deallocating Matrices
    * 
    * 'allocate' and 'deallocate' are the general purpose procedures for allocating and
    * deallocating matrices (loosely resembling the Fortran90 allocate/deallocate vernacular).
    * Matrix memory is allocated on the heap and owned by the matrix.
    ***********************************************************************************************/

    // Matrix allocation
    template<class T, StorageOrder ORDER>
    void allocate(DynamicMatrix<T, ORDER>&, index_t, index_t);

    // Matrix deallocation
    template<class T, StorageOrder ORDER>
    void deallocate(DynamicMatrix<T, ORDER>&);
}
#endif
