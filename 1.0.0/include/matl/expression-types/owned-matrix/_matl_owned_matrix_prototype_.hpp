#ifndef _MATL_OWNED_MATRIX_PROTOTYPE__HPP
#define _MATL_OWNED_MATRIX_PROTOTYPE__HPP

namespace matl
{
   /***********************************************************************************************
    * class _OWNED_MATRIX_
    *
    * Elementwise expression that owns a matrix. Temporary (rvalue) matrices used as operands, e.g.
    * 'A + eye<double>(3)', are moved into an _OWNED_MATRIX_ so that the expression remains valid
    * after the temporary would otherwise have been destroyed.
    *
    * The matrix is held through a shared pointer to const, so copies of the expression (and of
    * expressions containing it) share the matrix instead of copying its data. The matrix is
    * never modified after it is moved in.
    *
    * Template Parameters:
    *     MAT        Matrix type (e.g. RowMajorMatrix<T>, FixedSizeRowMajorMatrix<T,M,N>)
    ***********************************************************************************************/
    template<typename MAT>
    class _OWNED_MATRIX_ : public _EXPRESSION_<_OWNED_MATRIX_<MAT> >
    {
      public:
        // Typedefs
        typedef typename type_traits<MAT>::dtype dtype;

      private:
        // Owned matrix
        std::shared_ptr<const MAT> mat_;

      public:

        // Constructor (takes ownership of the matrix)
        explicit _OWNED_MATRIX_(MAT mat);

        // Check if evaluation target is safe
        template<typename MAT2>
        constexpr bool unsafe_evaluation_target(_MATRIX_<MAT2> const&) const;

        // Check if expression references address of target matrix
        template<typename MAT2>
        constexpr bool references_address(_MATRIX_<MAT2> const*) const;

        // Evaluate the expression to the target matrix
        template<typename MAT2>
        void eval(_MATRIX_<MAT2>&) const;

        // Get the number of rows
        index_t rows() const;

        // Get the number of columns
        index_t cols() const;

        // Get the number of elements
        index_t size() const;

        // Get the element at the i-th row and j-th column
        dtype operator() (index_t, index_t) const;

        // Get the element at the i-th index
        dtype operator() (index_t) const;
    };
}

#endif
