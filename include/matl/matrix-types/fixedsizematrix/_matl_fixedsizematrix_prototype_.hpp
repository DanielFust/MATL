#ifndef MATL_FIXEDSIZEMATRIX_PROTOTYPE_HPP
#define MATL_FIXEDSIZEMATRIX_PROTOTYPE_HPP

namespace matl 
{
    namespace detail
    {
        // Traits of FixedSizeMatrix (see _matl_type_traits_.hpp)
        template<class T, index_t M, index_t N, StorageOrder ORDER>
        struct traits< FixedSizeMatrix<T, M, N, ORDER> >
        {
            typedef T dtype;
            static constexpr bool         fixed_size    = true;
            static constexpr index_t      rows          = M;
            static constexpr index_t      cols          = N;
            static constexpr StorageOrder storage_order = ORDER;
        };
    }

   /***********************************************************************************************
    * class FixedSizeMatrix<T, M, N, ORDER>
    * 
    * Matrix whose dimensions M x N are known at compile time. The elements are stored in the
    * object itself (no heap memory), row by row (StorageOrder::RowMajor) or column by column
    * (StorageOrder::ColMajor). Usually named through the aliases FixedSizeMatrix<T,M,N>
    * (= Matrix<T,M,N>) and FixedSizeColMajorMatrix<T,M,N>.
    *
    * Template Arguments:
    *     T         Matrix data type
    *     M, N      Number of rows and columns
    *     ORDER     Storage order (RowMajor or ColMajor)
    **********************************************************************************************/
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    class FixedSizeMatrix : public MatrixBase< FixedSizeMatrix<T, M, N, ORDER> >
    {
        static_assert(ORDER == StorageOrder::RowMajor || ORDER == StorageOrder::ColMajor,
                      "FixedSizeMatrix: the storage order must be RowMajor or ColMajor");

        // Position of element (i,j) (1-based) in the data array
        static constexpr index_t offset(index_t i, index_t j)
        {
            if constexpr (ORDER == StorageOrder::RowMajor) return N * (i - 1) + (j - 1);
            else                                           return M * (j - 1) + (i - 1);
        }

      private:

        // Matrix data
        T data_[M*N];  

      public:  

        // Static row getter
        static constexpr index_t ROWS = M;

        // Static col getter
        static constexpr index_t COLS = N;

        // Static size getter
        static constexpr index_t SIZE = M*N;

        // Matrix datatype
        typedef T dtype;

        /********************* Constructors *********************/

        // Default constructor (elements are uninitialized)
        FixedSizeMatrix() = default;

        // Copy constructor
        FixedSizeMatrix(FixedSizeMatrix const&) = default;

        // Initializer-list constructor, one inner list per row (must be M rows of N elements)
        FixedSizeMatrix(matrix_rows<T>);

        // The members constrained by SizesMatch below accept only matrices and expressions that
        // can have dimensions M x N (others do not compile; dynamic dimensions are checked at run
        // time). They are defined here and forward to the out-of-line assign_expression and
        // copy_matrix, because Clang before version 16 cannot match an out-of-line definition
        // whose constraints use the class's template parameters.

        // Expression constructor (evaluates the expression; implicit, so that
        // 'Matrix<double,M,N> m = A + B;' evaluates to a new matrix)
        template<ExpressionType EXPR> requires SizesMatch<FixedSizeMatrix, EXPR>
        FixedSizeMatrix(EXPR const& expr) {this->assign_expression(expr);}

        /**************** Assignment & Copying *****************/
        
        // Matrix copy (general)
        template<MatrixType MAT> requires SizesMatch<FixedSizeMatrix, MAT>
        FixedSizeMatrix& operator= (MAT const& mat) {this->copy_matrix(mat); return *this;}

        // Matrix copy (FixedSizeMatrix specialization)
        FixedSizeMatrix& operator= (FixedSizeMatrix const&) = default;

        // Initializer-list assignment, one inner list per row (must be M rows of N elements)
        FixedSizeMatrix& operator= (matrix_rows<T>);

        // Matrix copy (DynamicMatrix with the same storage order: copied in storage order)
        template<class T2>
        FixedSizeMatrix& operator= (DynamicMatrix<T2, ORDER> const&);

        // Expression evaluation trigger
        template<ExpressionType EXPR> requires SizesMatch<FixedSizeMatrix, EXPR>
        FixedSizeMatrix& operator=(EXPR const& expr) {this->assign_expression(expr); return *this;}

        // Copy 
        template<MatrixType MAT> requires SizesMatch<FixedSizeMatrix, MAT>
        void copy(MAT const& mat) {this->copy_matrix(mat);}


        /********************* Data Access *********************/

        // Single index getter
        dtype operator() (index_t) const;

        // Double index getter
        dtype operator() (index_t, index_t) const;

        // Single index accessor
        dtype& operator() (index_t);

        // Double index accessor
        dtype& operator() (index_t, index_t);

        /********************* Utilities *********************/

        // Row Getter
        index_t rows() const;

        // Coloumn Getter
        index_t cols() const;

        // Size Getter
        index_t size() const;

        // Scalar assignment (sets every element)
        template<Scalar S>
        FixedSizeMatrix& operator= (S val);

      private:

        // Evaluate an expression into this matrix
        template<ExpressionType EXPR>
        void assign_expression(EXPR const&);

        // Copy the elements of a matrix with the same dimensions
        template<MatrixType MAT>
        void copy_matrix(MAT const&);

        // Copy the elements of an initializer list with the same dimensions
        void copy_rows(matrix_rows<T>);
    }; 

    // Fixed-size matrix "allocation" (checks the requested dimensions match the fixed size)
    template<class T, index_t M, index_t N, StorageOrder ORDER>
    void allocate(FixedSizeMatrix<T, M, N, ORDER>&, index_t, index_t);
}
#endif