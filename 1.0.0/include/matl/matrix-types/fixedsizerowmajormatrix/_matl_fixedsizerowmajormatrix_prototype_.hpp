#ifndef _MATL_FIXEDSIZEROWMAJORMATRIX_PROTOTYPE__HPP
#define _MATL_FIXEDSIZEROWMAJORMATRIX_PROTOTYPE__HPP

namespace matl 
{
   /***********************************************************************************************
    * class FixedSizeRowMajorMatrix
    * 
    * Fixed size row major matrix class. This class is an interface for raw arrays, with size
    * MxN known at compile time, which are treated as matrices with row-major ordering.
    **********************************************************************************************/
    template<class T, index_t M, index_t N>
    class FixedSizeRowMajorMatrix : public _FIXED_SIZE_MATRIX_<FixedSizeRowMajorMatrix<T, M, N> >
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
        // 'Matrix<double,M,N> m = A + B;' evaluates to a new matrix)
        template<typename EXPR>
        FixedSizeRowMajorMatrix(_EXPRESSION_<EXPR> const&);

        /**************** Assignment & Copying *****************/
        
        // Matrix copy (general)
        template<typename MAT>
        FixedSizeRowMajorMatrix& operator= (_MATRIX_<MAT>&);

        // Matrix copy (FixedSizeRowMajorMatrix specialization)
        FixedSizeRowMajorMatrix& operator= (FixedSizeRowMajorMatrix<T, M, N>&);

        // Matrix copy (RowMajorMatrix specialization)
        template<class T2>
        FixedSizeRowMajorMatrix& operator= (RowMajorMatrix<T2>&);

        // Expression evaluation trigger
        template<typename EXPR>
        FixedSizeRowMajorMatrix& operator=(_EXPRESSION_<EXPR> const& expr);


        // Copy 
        template<typename MAT>
        void copy(const _MATRIX_<MAT>& mat);


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

        // Scalar data writing
        template<class T2, enable_if_scalar<T2> = true>
        void operator= (T2 val);
    }; 

    // Fixed-size matrix "allocation" (checks the requested dimensions match the fixed size)
    template<class T, index_t M, index_t N>
    void allocate(FixedSizeRowMajorMatrix<T, M, N>&, index_t, index_t);
}
#endif