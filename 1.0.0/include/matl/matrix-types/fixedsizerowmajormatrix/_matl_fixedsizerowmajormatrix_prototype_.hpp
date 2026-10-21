#ifndef MATL_FIXEDSIZEROWMAJORMATRIX_PROTOTYPE_HPP
#define MATL_FIXEDSIZEROWMAJORMATRIX_PROTOTYPE_HPP

namespace matl 
{
    namespace detail
    {
        // Traits of FixedSizeRowMajorMatrix (see _matl_type_traits_.hpp)
        template<class T, index_t M, index_t N>
        struct traits< FixedSizeRowMajorMatrix<T, M, N> >
        {
            typedef T dtype;
            static constexpr bool         fixed_size    = true;
            static constexpr index_t      rows          = M;
            static constexpr index_t      cols          = N;
            static constexpr StorageOrder storage_order = StorageOrder::RowMajor;
        };
    }

   /***********************************************************************************************
    * class FixedSizeRowMajorMatrix
    * 
    * Fixed size row major matrix class. This class is an interface for raw arrays, with size
    * MxN known at compile time, which are treated as matrices with row-major ordering.
    **********************************************************************************************/
    template<class T, index_t M, index_t N>
    class FixedSizeRowMajorMatrix : public MatrixBase< FixedSizeRowMajorMatrix<T, M, N> >
    {
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
        FixedSizeRowMajorMatrix() = default;

        // Copy constructor
        FixedSizeRowMajorMatrix(FixedSizeRowMajorMatrix const&) = default;

        // Expression constructor (evaluates the expression; implicit, so that
        // 'Matrix<double,M,N> m = A + B;' evaluates to a new matrix). Expressions with other
        // compile-time dimensions are rejected (SizesMatch); dynamic ones are checked at run time.
        template<ExpressionType EXPR> requires SizesMatch<FixedSizeRowMajorMatrix, EXPR>
        FixedSizeRowMajorMatrix(EXPR const&);

        /**************** Assignment & Copying *****************/
        
        // Matrix copy (general)
        template<MatrixType MAT> requires SizesMatch<FixedSizeRowMajorMatrix, MAT>
        FixedSizeRowMajorMatrix& operator= (MAT const&);

        // Matrix copy (FixedSizeRowMajorMatrix specialization)
        FixedSizeRowMajorMatrix& operator= (FixedSizeRowMajorMatrix const&) = default;

        // Matrix copy (RowMajorMatrix specialization)
        template<class T2>
        FixedSizeRowMajorMatrix& operator= (RowMajorMatrix<T2> const&);

        // Expression evaluation trigger
        template<ExpressionType EXPR> requires SizesMatch<FixedSizeRowMajorMatrix, EXPR>
        FixedSizeRowMajorMatrix& operator=(EXPR const& expr);

        // Copy 
        template<MatrixType MAT> requires SizesMatch<FixedSizeRowMajorMatrix, MAT>
        void copy(MAT const& mat);


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
        FixedSizeRowMajorMatrix& operator= (S val);
    }; 

    // Fixed-size matrix "allocation" (checks the requested dimensions match the fixed size)
    template<class T, index_t M, index_t N>
    void allocate(FixedSizeRowMajorMatrix<T, M, N>&, index_t, index_t);
}
#endif