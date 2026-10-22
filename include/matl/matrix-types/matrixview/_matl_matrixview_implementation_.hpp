/**************************************************************************************************
 * File: _matl_matrixview_implementation_.hpp
 *
 * Description: Implementation of MatrixView (see _matl_matrixview_prototype_.hpp).
 **************************************************************************************************/
#ifndef MATL_MATRIXVIEW_IMPLEMENTATION_HPP
#define MATL_MATRIXVIEW_IMPLEMENTATION_HPP
namespace matl
{
    namespace detail
    {
       /*******************************************************************************************
        * void assign_to(TARGET& target, X const& x)
        *
        * Writes a matrix, view or expression into a target that cannot be reallocated (a view):
        * the dimensions must be the target's. If the assignment would overwrite elements of the
        * target that it still has to read (e.g. A(range(2, last), all) = A(range(1, last - 1), all)
        * shifts rows that overlap), the value is first evaluated into a temporary.
        *******************************************************************************************/
        template<WritableMatrix TARGET, MatrixAlgebraObject X>
        void assign_to(TARGET& target, X const& x)
        {
            #ifdef MATL_DIMENSION_CHECK
            check_dimensions(x);
            if ((target.rows() != x.rows()) || (target.cols() != x.cols()))
                assignment_dimensions(target.rows(), target.cols(), x.rows(), x.cols());
            #endif

            // Copies a matrix, or evaluates an expression, into a target
            const auto write = [&x](auto& to)
            {
                if constexpr (MatrixType<X>) assign_elementwise(to, x);
                else                         x.eval(to);
            };

            const Footprint target_footprint = footprint(target);
            bool unsafe;
            if constexpr (MatrixType<X>) unsafe = unsafe_operand(target_footprint, x);
            else                         unsafe = x.unsafe_evaluation_target(target_footprint);

            if (unsafe)
            {
                constexpr StorageOrder order = (traits<X>::storage_order == StorageOrder::ColMajor)
                                               ? StorageOrder::ColMajor : StorageOrder::RowMajor;
                EvaluationTemporary< DynamicMatrix<dtype_t<X>, order> > temp(x.rows(), x.cols());
                write(temp.matrix());
                assign_elementwise(target, temp.matrix());
            }
            else
                write(target);
        }
    }


    /*=============================================================================================*
     *                                         CONSTRUCTORS
     *=============================================================================================*/

   /***********************************************************************************************
    * MatrixView(origin, rows, cols, row_stride, col_stride, owner)
    *
    * View of rows x cols elements: element (i,j) is origin[(i-1)*row_stride + (j-1)*col_stride].
    *
    * Parameters:
    *     origin       Address of element (1,1)
    *     rows, cols   Dimensions
    *     row_stride   Elements from one row of the view to the next
    *     col_stride   Elements from one column of the view to the next
    *     owner        The viewed matrix, if the view keeps it alive (empty otherwise)
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    MatrixView<T, ORDER>::MatrixView(T* origin, index_t rows, index_t cols,
                                     std::ptrdiff_t row_stride, std::ptrdiff_t col_stride,
                                     std::shared_ptr<const void> owner) :
        origin_(origin),
        rows_(rows),
        cols_(cols),
        row_stride_(row_stride),
        col_stride_(col_stride),
        owner_(std::move(owner))
    {}


    /*=============================================================================================*
     *                                         ASSIGNMENT
     *=============================================================================================*/

   /***********************************************************************************************
    * MatrixView& operator= (MatrixView const&)
    *
    * Writes the elements of another view of the same type into the viewed elements.
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    MatrixView<T, ORDER>& MatrixView<T, ORDER>::operator= (MatrixView const& other)
    {
        static_assert(!std::is_const_v<T>, "MatrixView: a view of a const matrix is read-only");
        detail::assign_to(*this, other);
        return *this;
    }

   /***********************************************************************************************
    * MatrixView& operator= (X const&)
    *
    * Writes a matrix, view or expression (with the view's dimensions) into the viewed elements.
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    template<MatrixAlgebraObject X>
    MatrixView<T, ORDER>& MatrixView<T, ORDER>::operator= (X const& x)
    {
        static_assert(!std::is_const_v<T>, "MatrixView: a view of a const matrix is read-only");
        detail::assign_to(*this, x);
        return *this;
    }

   /***********************************************************************************************
    * MatrixView& operator= (S)
    *
    * Sets every viewed element to a scalar.
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    template<Scalar S>
    MatrixView<T, ORDER>& MatrixView<T, ORDER>::operator= (S value)
    {
        static_assert(!std::is_const_v<T>, "MatrixView: a view of a const matrix is read-only");
        this->for_each([value](T& x) {x = value;});
        return *this;
    }

   /***********************************************************************************************
    * Arithmetic assignment by a scalar (every viewed element)
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    template<Scalar S>
    MatrixView<T, ORDER>& MatrixView<T, ORDER>::operator+= (S value)
    {
        static_assert(!std::is_const_v<T>, "MatrixView: a view of a const matrix is read-only");
        this->for_each([value](T& x) {x += value;});
        return *this;
    }

    template<class T, StorageOrder ORDER>
    template<Scalar S>
    MatrixView<T, ORDER>& MatrixView<T, ORDER>::operator-= (S value)
    {
        static_assert(!std::is_const_v<T>, "MatrixView: a view of a const matrix is read-only");
        this->for_each([value](T& x) {x -= value;});
        return *this;
    }

    template<class T, StorageOrder ORDER>
    template<Scalar S>
    MatrixView<T, ORDER>& MatrixView<T, ORDER>::operator*= (S value)
    {
        static_assert(!std::is_const_v<T>, "MatrixView: a view of a const matrix is read-only");
        this->for_each([value](T& x) {x *= value;});
        return *this;
    }

    // As for matrices, division by a floating-point scalar is multiplication by its reciprocal
    template<class T, StorageOrder ORDER>
    template<Scalar S>
    MatrixView<T, ORDER>& MatrixView<T, ORDER>::operator/= (S value)
    {
        static_assert(!std::is_const_v<T>, "MatrixView: a view of a const matrix is read-only");
        if constexpr (std::is_integral_v<dtype>)
            this->for_each([value](T& x) {x /= value;});
        else
        {
            const dtype reciprocal = dtype(1) / value;
            this->for_each([reciprocal](T& x) {x *= reciprocal;});
        }
        return *this;
    }

   /***********************************************************************************************
    * Arithmetic assignment by a matrix, view or expression (*= is the matrix product, so X must
    * be square)
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    template<MatrixAlgebraObject X>
    MatrixView<T, ORDER>& MatrixView<T, ORDER>::operator+= (X const& x)
    {return *this = *this + x;}

    template<class T, StorageOrder ORDER>
    template<MatrixAlgebraObject X>
    MatrixView<T, ORDER>& MatrixView<T, ORDER>::operator-= (X const& x)
    {return *this = *this - x;}

    template<class T, StorageOrder ORDER>
    template<MatrixAlgebraObject X>
    MatrixView<T, ORDER>& MatrixView<T, ORDER>::operator*= (X const& x)
    {return *this = *this * x;}

   /***********************************************************************************************
    * void for_each(F&& f) const
    *
    * Calls f(element) for every viewed element, in the storage order of the viewed matrix (so
    * that consecutive elements are as close in memory as the view allows).
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    template<typename F>
    void MatrixView<T, ORDER>::for_each(F&& f) const
    {
        // Outer and inner dimension, in storage order
        const bool     row_major    = (ORDER == StorageOrder::RowMajor);
        const index_t  outer        = row_major ? this->rows_ : this->cols_;
        const index_t  inner        = row_major ? this->cols_ : this->rows_;
        const std::ptrdiff_t outer_stride = row_major ? this->row_stride_ : this->col_stride_;
        const std::ptrdiff_t inner_stride = row_major ? this->col_stride_ : this->row_stride_;

        for (index_t a = 0; a < outer; ++a)
        {
            T* line = this->origin_ + static_cast<std::ptrdiff_t>(a) * outer_stride;
            for (index_t b = 0; b < inner; ++b)
                f(line[static_cast<std::ptrdiff_t>(b) * inner_stride]);
        }
    }


    /*=============================================================================================*
     *                                         DATA ACCESS
     *=============================================================================================*/

   /***********************************************************************************************
    * T& operator() (index_t i, index_t j) const
    *
    * Returns:
    *     Element (i,j) of the view (1-based)
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    inline T& MatrixView<T, ORDER>::operator() (index_t i, index_t j) const
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((i < 1) || (i > this->rows_) || (j < 1) || (j > this->cols_))
            detail::index_out_of_range(i, j, this->rows_, this->cols_);
        #endif
        return this->origin_[static_cast<std::ptrdiff_t>(i - 1) * this->row_stride_ +
                             static_cast<std::ptrdiff_t>(j - 1) * this->col_stride_];
    }

   /***********************************************************************************************
    * T& operator() (index_t k) const
    *
    * Returns:
    *     Element k of the view in storage order (1-based); for a row or column vector, element k
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    inline T& MatrixView<T, ORDER>::operator() (index_t k) const
    {
        #ifdef MATL_BOUNDS_CHECK
        if ((k < 1) || (k > this->size()))
            detail::index_out_of_range(k, this->size());
        #endif
        if (this->rows_ == 1) return (*this)(1, k);
        if (this->cols_ == 1) return (*this)(k, 1);
        if constexpr (ORDER == StorageOrder::RowMajor)
            return (*this)((k - 1) / this->cols_ + 1, (k - 1) % this->cols_ + 1);
        else
            return (*this)((k - 1) % this->rows_ + 1, (k - 1) / this->rows_ + 1);
    }


    /*=============================================================================================*
     *                                           UTILITIES
     *=============================================================================================*/

    template<class T, StorageOrder ORDER>
    inline index_t MatrixView<T, ORDER>::rows() const
    {return this->rows_;}

    template<class T, StorageOrder ORDER>
    inline index_t MatrixView<T, ORDER>::cols() const
    {return this->cols_;}

    template<class T, StorageOrder ORDER>
    inline index_t MatrixView<T, ORDER>::size() const
    {return this->rows_ * this->cols_;}


    /*=============================================================================================*
     *                                     EXPRESSION INTERFACE
     *=============================================================================================*/

   /***********************************************************************************************
    * detail::Footprint footprint() const
    *
    * Returns:
    *     Where the viewed elements are in memory (strides in bytes)
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    detail::Footprint MatrixView<T, ORDER>::footprint() const
    {
        const std::ptrdiff_t size = sizeof(T);
        return detail::Footprint{reinterpret_cast<std::uintptr_t>(this->origin_), this->rows_, this->cols_,
                                 this->row_stride_ * size, this->col_stride_ * size, sizeof(T)};
    }

   /***********************************************************************************************
    * bool unsafe_evaluation_target(detail::Footprint const& target) const
    *
    * Copying the view into the target (given by its footprint) is safe unless the view shares elements with the target
    * at other positions (see detail::unsafe_operand).
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    bool MatrixView<T, ORDER>::unsafe_evaluation_target(detail::Footprint const& target) const
    {return detail::unsafe_operand(target, *this);}

   /***********************************************************************************************
    * bool reads(detail::Footprint const& target) const
    *
    * Returns:
    *     'true' if the view may share elements with the target
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    bool MatrixView<T, ORDER>::reads(detail::Footprint const& target) const
    {return detail::overlap(target, this->footprint());}

   /***********************************************************************************************
    * void eval(MAT& target) const
    *
    * Copies the viewed elements into the target (a safe target with the view's dimensions).
    ***********************************************************************************************/
    template<class T, StorageOrder ORDER>
    template<WritableMatrix MAT>
    void MatrixView<T, ORDER>::eval(MAT& target) const
    {detail::assign_elementwise(target, *this);}


    /*=============================================================================================*
     *                                             SWAP
     *=============================================================================================*/

    namespace detail
    {
       /*******************************************************************************************
        * void swap_views(MatrixView a, MatrixView b)
        *
        * Exchanges the elements of two views with the same dimensions (see swap). Views that do
        * not share memory exchange their elements directly; overlapping views go through a
        * temporary (t = a; a = b; b = t).
        *******************************************************************************************/
        template<class T, StorageOrder ORDER_A, StorageOrder ORDER_B>
        void swap_views(MatrixView<T, ORDER_A> a, MatrixView<T, ORDER_B> b)
        {
            static_assert(!std::is_const_v<T>, "swap: a view of a const matrix is read-only");

            const index_t M = a.rows();
            const index_t N = a.cols();
            #ifdef MATL_DIMENSION_CHECK
            if ((M != b.rows()) || (N != b.cols()))
                incompatible_dimensions("swap", M, N, b.rows(), b.cols());
            #endif

            const Footprint footprint_a = a.footprint();
            const Footprint footprint_b = b.footprint();
            if (same_elements(footprint_a, footprint_b)) return;

            if (overlap(footprint_a, footprint_b))
            {
                EvaluationTemporary< DynamicMatrix<T, ORDER_A> > temp(M, N);
                assign_elementwise(temp.matrix(), a);
                a = b;
                b = temp.matrix();
                return;
            }

            // In the storage order of a
            if constexpr (ORDER_A == StorageOrder::RowMajor)
            {
                for (index_t i = 1; i <= M; ++i)
                    for (index_t j = 1; j <= N; ++j)
                        std::swap(a(i,j), b(i,j));
            }
            else
            {
                for (index_t j = 1; j <= N; ++j)
                    for (index_t i = 1; i <= M; ++i)
                        std::swap(a(i,j), b(i,j));
            }
        }
    }

    // Two views of the same type (preferred to std::swap for 'using std::swap; swap(a, b);')
    template<class T, StorageOrder ORDER>
    void swap(MatrixView<T, ORDER>& a, MatrixView<T, ORDER>& b)
    {detail::swap_views(a, b);}

    // Any two views of the same element type, including temporaries: swap(A(1, all), A(3, all))
    template<class T, StorageOrder ORDER_A, StorageOrder ORDER_B>
    void swap(MatrixView<T, ORDER_A> a, MatrixView<T, ORDER_B> b)
    {detail::swap_views(a, b);}


    namespace detail
    {
       /*******************************************************************************************
        * struct evaluator<MatrixView>
        *
        * Reads elements of a view by (row, column). Views have no single-index access in storage
        * order (see linear_index_v), so neither does this evaluator.
        *******************************************************************************************/
        template<class T, StorageOrder ORDER>
        struct evaluator< MatrixView<T, ORDER> >
        {
            typedef std::remove_const_t<T> dtype;

            MatrixView<T, ORDER> const& view_;

            explicit evaluator(MatrixView<T, ORDER> const& view) : view_(view) {}

            dtype operator() (index_t i, index_t j) const {return view_(i,j);}
        };
    }
}
#endif
