/**************************************************************************************************
 * File: _matl_matrixview_prototype_.hpp
 *
 * Description: Views: rectangular selections of a matrix's rows and columns that read and write
 *              the matrix's own elements (nothing is copied).
 *
 *     A(range(1, 3), all)            rows 1 to 3                     MATLAB A(1:3, :)
 *     A(all, 2)                      column 2                        MATLAB A(:, 2)
 *     A(4, range(2, last))           row 4, columns 2 to the end     MATLAB A(4, 2:end)
 *     A(range(1, 2, last), all)      every other row                 MATLAB A(1:2:end, :)
 *
 *     A(all, 1) = 0.0;                           writes into A
 *     A(range(1, 2), range(1, 2)) = B * C;       an expression evaluated into part of A
 *     double s = sum(A(all, j));                 a view is used like a matrix
 *
 * A view stores the address of its element (1,1), its dimensions, and two strides: the distance
 * in memory from one row of the view to the next, and from one column to the next. Element (i,j)
 * is at origin[(i-1)*row_stride + (j-1)*col_stride]. A matrix is the special case with strides
 * (cols, 1) (row-major) or (1, rows) (column-major); a view of rows 2 to 3 and columns 2 to 4 of
 * a row-major 4x5 matrix starts at A(2,2) with strides (5, 1), and a step multiplies the stride.
 *
 * A view refers to the matrix's memory, like a stored expression: it is valid only while the
 * matrix keeps its memory (it must not be destroyed, moved from, or reallocated). A view of a
 * temporary matrix, e.g. ones<double>(3, 3)(all, 1), keeps the matrix alive itself. Views of
 * const matrices (and of temporaries) are read-only.
 *
 * Copying a view gives a second view of the same elements; assigning to a view writes elements:
 *
 *     auto v = A(all, 1);            // v is a view of column 1 of A
 *     v = 0.0;                       // sets column 1 of A to 0
 *     Matrix<double> c = A(all, 1);  // c is a new matrix (a copy of column 1)
 **************************************************************************************************/
#ifndef MATL_MATRIXVIEW_PROTOTYPE_HPP
#define MATL_MATRIXVIEW_PROTOTYPE_HPP
namespace matl
{
    namespace detail
    {
        // Traits of MatrixView (see _matl_type_traits_.hpp): dynamic dimensions, and the storage
        // order of the viewed matrix (the order in which its elements are cheapest to visit)
        template<class T, StorageOrder ORDER>
        struct traits< MatrixView<T, ORDER> >
        {
            typedef std::remove_const_t<T> dtype;
            static constexpr bool         fixed_size    = false;
            static constexpr index_t      rows          = 0;
            static constexpr index_t      cols          = 0;
            static constexpr StorageOrder storage_order = ORDER;
        };

        template<class T, StorageOrder ORDER>
        inline constexpr bool is_view_v< MatrixView<T, ORDER> > = true;

        template<class T, StorageOrder ORDER>
        inline constexpr bool is_writable_view_v< MatrixView<T, ORDER> > = !std::is_const_v<T>;

        template<class T, StorageOrder ORDER>
        inline constexpr bool linear_index_v< MatrixView<T, ORDER> > = false;


       /*******************************************************************************************
        * struct Footprint
        *
        * Where the elements of a matrix or view are in memory. An assignment compares the
        * footprint of its target with those of the matrices and views it reads, to find out
        * whether it would overwrite elements that it still has to read (see unsafe_operand).
        *******************************************************************************************/
        struct Footprint
        {
            std::uintptr_t origin       = 0;    // address of element (1,1)
            index_t        rows         = 0;
            index_t        cols         = 0;
            std::ptrdiff_t row_stride   = 0;    // bytes from element (i,j) to (i+1,j)
            std::ptrdiff_t col_stride   = 0;    // bytes from element (i,j) to (i,j+1)
            std::size_t    element_size = 0;

            bool empty() const {return this->rows == 0 || this->cols == 0;}

            // Lowest address of an element
            std::uintptr_t low() const
            {
                std::ptrdiff_t offset = 0;
                if (this->row_stride < 0) offset += this->row_stride * static_cast<std::ptrdiff_t>(this->rows - 1);
                if (this->col_stride < 0) offset += this->col_stride * static_cast<std::ptrdiff_t>(this->cols - 1);
                return this->origin + static_cast<std::uintptr_t>(offset);
            }

            // One past the last byte of the element with the highest address
            std::uintptr_t high() const
            {
                std::ptrdiff_t offset = 0;
                if (this->row_stride > 0) offset += this->row_stride * static_cast<std::ptrdiff_t>(this->rows - 1);
                if (this->col_stride > 0) offset += this->col_stride * static_cast<std::ptrdiff_t>(this->cols - 1);
                return this->origin + static_cast<std::uintptr_t>(offset) + this->element_size;
            }
        };

        // 'true' if a and b may share elements. The test compares the address ranges, so it is
        // conservative: e.g. two columns of a row-major matrix interleave without sharing
        // elements, but count as overlapping (the assignment then uses a temporary).
        inline bool overlap(Footprint const& a, Footprint const& b)
        {return !a.empty() && !b.empty() && a.low() < b.high() && b.low() < a.high();}

        // 'true' if element (i,j) of a is element (i,j) of b, for every (i,j)
        inline bool same_elements(Footprint const& a, Footprint const& b)
        {
            return a.origin == b.origin && a.rows == b.rows && a.cols == b.cols &&
                   a.element_size == b.element_size &&
                   (a.rows <= 1 || a.row_stride == b.row_stride) &&
                   (a.cols <= 1 || a.col_stride == b.col_stride);
        }

        // Footprint of a matrix or view. Only the address and dimensions are read, not the
        // elements. The parameter is a forwarding reference so that an assignment target is
        // passed as non-const: GCC's -Wmaybe-uninitialized assumes a function reads an object
        // passed by const reference, and a fixed-size target is uninitialized before its first
        // assignment (Matrix<double, 3, 3> F; F = A + B;).
        template<typename X>
        Footprint footprint(X&& x)
        {
            typedef std::remove_cvref_t<X> type;
            if constexpr (is_view_v<type>)
                return x.footprint();
            else
            {
                const std::ptrdiff_t size = sizeof(dtype_t<type>);
                const bool row_major = (traits<type>::storage_order == StorageOrder::RowMajor);
                return Footprint{reinterpret_cast<std::uintptr_t>(x.data()), x.rows(), x.cols(),
                                 row_major ? static_cast<std::ptrdiff_t>(x.cols()) * size : size,
                                 row_major ? size : static_cast<std::ptrdiff_t>(x.rows()) * size,
                                 sizeof(dtype_t<type>)};
            }
        }

       /*******************************************************************************************
        * MatrixView<T, ORDER> make_view(origin, rows, cols, row_stride, col_stride, r, c, owner)
        *
        * View of the rows selected by r and the columns selected by c (see _matl_indexing_.hpp)
        * of the rows x cols elements at 'origin' with the given strides (in elements). Used by
        * the matrices and by views (for views of views).
        *
        * Parameters:
        *     owner        Keeps a temporary matrix alive for the view (empty otherwise)
        *******************************************************************************************/
        template<StorageOrder ORDER, class T, typename R, typename C>
        MatrixView<T, ORDER> make_view(T* origin, index_t rows, index_t cols,
                                       std::ptrdiff_t row_stride, std::ptrdiff_t col_stride,
                                       R const& r, C const& c, std::shared_ptr<const void> owner = nullptr)
        {
            const Selection sr = select_indices(r, rows, "row");
            const Selection sc = select_indices(c, cols, "column");

            // An empty view has no element (1,1) (its first index may be outside the matrix)
            T* first = origin;
            if (sr.count > 0 && sc.count > 0)
                first = origin + (sr.first - 1) * row_stride + (sc.first - 1) * col_stride;

            return MatrixView<T, ORDER>(first, sr.count, sc.count, sr.step * row_stride, sc.step * col_stride,
                                        std::move(owner));
        }
    }


   /***********************************************************************************************
    * class MatrixView<T, ORDER>
    *
    * View of rows and columns of a matrix (see the description above). Made by indexing a matrix
    * (or a view) with all, a range, or last; not constructed directly.
    *
    * A view reads like an expression (it can be an operand of any operation, and is stored by
    * value in expressions, which is cheap) and writes like a matrix (it can be assigned to; an
    * assignment must have the view's dimensions, checked when MATL_DIMENSION_CHECK is defined).
    *
    * Template Parameters:
    *     T         Element type (const T for a read-only view)
    *     ORDER     Storage order of the viewed matrix (RowMajor or ColMajor)
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    class MatrixView : public ExpressionBase< MatrixView<T, ORDER> >
    {
        static_assert(ORDER == StorageOrder::RowMajor || ORDER == StorageOrder::ColMajor,
                      "MatrixView: the storage order must be RowMajor or ColMajor");

      public:
        // Element type
        typedef std::remove_const_t<T> dtype;

      private:
        T*             origin_;          // element (1,1)
        index_t        rows_;
        index_t        cols_;
        std::ptrdiff_t row_stride_;      // elements from (i,j) to (i+1,j)
        std::ptrdiff_t col_stride_;      // elements from (i,j) to (i,j+1)

        // The viewed matrix, if the view keeps it alive (views of temporary matrices)
        std::shared_ptr<const void> owner_;

        // Calls f(element) for every element, in storage order
        template<typename F>
        void for_each(F&& f) const;

      public:

        /********************* Constructors *********************/

        // View of the given elements (see detail::make_view)
        MatrixView(T* origin, index_t rows, index_t cols, std::ptrdiff_t row_stride, std::ptrdiff_t col_stride,
                   std::shared_ptr<const void> owner = nullptr);

        // Copy: another view of the same elements
        MatrixView(MatrixView const&) = default;

        /********************* Assignment ***********************/
        // Assignments write the viewed elements; they do not make the view refer to other
        // elements. Views of const matrices cannot be assigned to.

        // Assignment from a view of the same type (writes its elements)
        MatrixView& operator= (MatrixView const&);

        // Assignment from a matrix, view or expression
        template<MatrixAlgebraObject X>
        MatrixView& operator= (X const&);

        // Scalar assignment (sets every element)
        template<Scalar S>
        MatrixView& operator= (S);

        // Arithmetic assignment (by scalar)
        template<Scalar S> MatrixView& operator+= (S);
        template<Scalar S> MatrixView& operator-= (S);
        template<Scalar S> MatrixView& operator*= (S);
        template<Scalar S> MatrixView& operator/= (S);

        // Arithmetic assignment (by matrix, view or expression; *= is the matrix product)
        template<MatrixAlgebraObject X> MatrixView& operator+= (X const&);
        template<MatrixAlgebraObject X> MatrixView& operator-= (X const&);
        template<MatrixAlgebraObject X> MatrixView& operator*= (X const&);

        /********************* Data Access *********************/
        // Element access does not depend on the view's constness (as for std::span): a const
        // view of a non-const matrix still writes the matrix.

        // Element (i,j)
        T& operator() (index_t, index_t) const;

        // Element k in storage order (for vectors: element k)
        T& operator() (index_t) const;

        // Element with last: v(last, 1), v(2, last - 1)
        template<typename R, typename C> requires detail::ElementIndices<R, C>
        T& operator() (R const& r, C const& c) const
        {return (*this)(detail::element_index(r, this->rows_), detail::element_index(c, this->cols_));}

        // View of this view: v(range(1, 2), all)
        template<typename R, typename C> requires detail::ViewIndices<R, C>
        MatrixView operator() (R const& r, C const& c) const
        {
            return detail::make_view<ORDER>(this->origin_, this->rows_, this->cols_,
                                            this->row_stride_, this->col_stride_, r, c, this->owner_);
        }

        /********************** Utilities **********************/

        // Number of rows
        index_t rows() const;

        // Number of columns
        index_t cols() const;

        // Number of elements
        index_t size() const;

        /********************* Expression interface **************/
        // (see ExpressionBase)

        // Where the elements are in memory
        detail::Footprint footprint() const;

        // Check if evaluating the view into the target would overwrite elements still to be read
        bool unsafe_evaluation_target(detail::Footprint const&) const;

        // Check if the view shares memory with the given footprint
        bool reads(detail::Footprint const&) const;

        // Copy the elements into the target
        template<WritableMatrix MAT>
        void eval(MAT&) const;
    };


   /***********************************************************************************************
    * swap(a, b)
    *
    * Exchanges the elements of two views with the same dimensions, e.g. two rows:
    *
    *     swap(A(1, all), A(3, all));
    *
    * std::swap(a, b) would instead exchange copies of the views (copying a view does not copy
    * elements), leaving both with b's elements, so call swap unqualified or as matl::swap. The
    * overload for two views of the same type is found instead of std::swap by
    * 'using std::swap; swap(a, b);'. If the views overlap, the result is that of t = a; a = b;
    * b = t.
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    void swap(MatrixView<T, ORDER>& a, MatrixView<T, ORDER>& b);

    template<class T, StorageOrder ORDER_A, StorageOrder ORDER_B>
    void swap(MatrixView<T, ORDER_A> a, MatrixView<T, ORDER_B> b);
}
#endif
