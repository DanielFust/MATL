/**************************************************************************************************
 * File:   _matl_type_traits_.hpp
 *
 * Description: Concepts and type traits for the MATL library.
 *
 * Every matrix and expression type describes itself with one specialization of detail::traits,
 * defined next to the type (e.g. traits<RowMajorMatrix<T>> in the RowMajorMatrix prototype):
 *
 *     dtype           Element type
 *     fixed_size      'true' if the dimensions are known at compile time
 *     rows, cols      Dimensions when fixed_size, otherwise 0
 *     storage_order   Order of the elements in memory (matrices), or the order in which an
 *                     expression's elements are cheapest to visit
 *
 * Scalars have traits too, so that they can be used as expression operands: they are fixed-size
 * and have no storage order of their own (StorageOrder::Any).
 **************************************************************************************************/
#ifndef MATL_TYPE_TRAITS_HPP
#define MATL_TYPE_TRAITS_HPP

namespace matl
{
    /****************************************** Concepts ******************************************/

    // Scalar: an arithmetic value (e.g. int, float, double)
    template<typename T>
    concept Scalar = std::is_arithmetic_v<T>;

    // MatrixType: a matrix (not an expression), ignoring references and cv-qualifiers
    template<typename T>
    concept MatrixType = std::derived_from< std::remove_cvref_t<T>, MatrixBase<std::remove_cvref_t<T>> >;

    // ExpressionType: an expression (formula) built from matrices, expressions and scalars
    template<typename T>
    concept ExpressionType = std::derived_from< std::remove_cvref_t<T>, ExpressionBase<std::remove_cvref_t<T>> >;

    // MatrixAlgebraObject: a matrix or an expression
    template<typename T>
    concept MatrixAlgebraObject = MatrixType<T> || ExpressionType<T>;

    // TemporaryMatrix: a deduced forwarding-reference type that binds a temporary (rvalue) matrix.
    // For 'template<TemporaryMatrix M> f(M&& m)', M is deduced as a non-reference type only for
    // rvalue arguments.
    template<typename T>
    concept TemporaryMatrix = MatrixType<T> && !std::is_reference_v<T>;

    // Operand: anything that can appear in a matrix algebra expression (matrix, expression or scalar)
    template<typename T>
    concept Operand = MatrixAlgebraObject<T> || Scalar<std::remove_cvref_t<T>>;


    namespace detail
    {
        /************************************* Type traits ****************************************/

        // Storage order (declared in matl, see _matl_matrix_types_declarations_.hpp)
        using StorageOrder = matl::StorageOrder;

        // Storage order of an expression combining operands with storage orders a and b
        constexpr StorageOrder combine(StorageOrder a, StorageOrder b)
        {
            if (a == StorageOrder::Any) return b;
            if (b == StorageOrder::Any) return a;
            return (a == b) ? a : StorageOrder::Mixed;
        }

        // Storage order of a transpose (element (i,j) of the transpose is element (j,i))
        constexpr StorageOrder transposed(StorageOrder order)
        {
            if (order == StorageOrder::RowMajor) return StorageOrder::ColMajor;
            if (order == StorageOrder::ColMajor) return StorageOrder::RowMajor;
            return order;
        }

        // Traits of a matrix, expression or scalar (specialized next to each type)
        template<typename T> struct traits;

        // Scalars
        template<Scalar S>
        struct traits<S>
        {
            typedef S dtype;
            static constexpr bool         fixed_size    = true;
            static constexpr index_t      rows          = 0;
            static constexpr index_t      cols          = 0;
            static constexpr StorageOrder storage_order = StorageOrder::Any;
        };

        // Element type
        template<typename T>
        using dtype_t = typename traits< std::remove_cvref_t<T> >::dtype;

        // Type promotion (the type of a sum of the given types, e.g. promote_t<int,double> = double)
        template<typename FIRST, typename... OTHERS>
        using promote_t = decltype( (std::declval<FIRST>() + ... + std::declval<OTHERS>()) );

        // 'true' if A and B have the same storage order, so their elements can be visited together
        // in storage order
        template<typename A, typename B>
        inline constexpr bool same_storage_order_v =
            (traits<A>::storage_order == traits<B>::storage_order) &&
            (traits<A>::storage_order == StorageOrder::RowMajor || traits<A>::storage_order == StorageOrder::ColMajor);

        // Expression kinds (specialized in the expression type prototypes)
        template<typename T> inline constexpr bool is_transpose_v = false;
        template<typename T> inline constexpr bool is_product_v   = false;
        template<typename T> inline constexpr bool is_owned_v     = false;

        // How an expression stores an operand of type T: matrices by reference (so that a stored
        // expression is a live formula), expressions and scalars by value
        template<typename T>
        using stored_t = std::conditional_t< MatrixType<T>, T const&, T >;


        /*********************************** Default matrix types *********************************/

        // Default matrix type: dynamic (no dimensions) or fixed-size (M x N)
        template<class T, index_t ... DIMS> struct default_matrix;

        template<class T>
        struct default_matrix<T>
        {typedef RowMajorMatrix<T> type;};

        template<class T, index_t M, index_t N>
        struct default_matrix<T, M, N>
        {typedef FixedSizeRowMajorMatrix<T, M, N> type;};

        // Matrix type that holds the value of T: T itself for matrices, otherwise a dynamic or
        // fixed-size matrix with T's element type and dimensions, column-major if T is
        // column-major and row-major otherwise
        template<typename T>
        struct result_matrix
        {
            static constexpr StorageOrder order = (traits<T>::storage_order == StorageOrder::ColMajor)
                                                  ? StorageOrder::ColMajor : StorageOrder::RowMajor;
            typedef std::conditional_t< traits<T>::fixed_size,
                                        FixedSizeMatrix<dtype_t<T>, traits<T>::rows, traits<T>::cols, order>,
                                        DynamicMatrix<dtype_t<T>, order> > type;
        };

        template<MatrixType MAT>
        struct result_matrix<MAT>
        {typedef MAT type;};

        template<typename T>
        using result_matrix_t = typename result_matrix<T>::type;
    }

    // Default matrix types: Matrix<T> (dynamic) and Matrix<T,M,N> (fixed-size)
    template<class T, index_t ... DIMS>
    using Matrix = typename detail::default_matrix<T, DIMS...>::type;

    // FixedSize: a matrix or expression whose dimensions are known at compile time
    template<typename T>
    concept FixedSize = MatrixAlgebraObject<T> && detail::traits< std::remove_cvref_t<T> >::fixed_size;


    /************************************* Dimension concepts *************************************/
    // These reject operands whose dimensions are known at compile time and do not fit, e.g. adding
    // a 2x3 and a 3x2 fixed-size matrix does not compile. Dynamic dimensions are only known at
    // run time, so operands with dynamic dimensions always satisfy them (they are checked when the
    // expression is evaluated, see MATL_DIMENSION_CHECK). A scalar fits any dimensions.

    namespace detail
    {
        // Dimensions of A and B are known at compile time (neither is a scalar)
        template<typename A, typename B>
        inline constexpr bool both_fixed_size_v = !Scalar<A> && !Scalar<B> &&
                                                  traits<A>::fixed_size && traits<B>::fixed_size;
    }

    // SizesMatch: A and B can have the same dimensions (elementwise operations, assignment)
    template<typename A, typename B>
    concept SizesMatch = !detail::both_fixed_size_v< std::remove_cvref_t<A>, std::remove_cvref_t<B> > ||
                         (detail::traits< std::remove_cvref_t<A> >::rows == detail::traits< std::remove_cvref_t<B> >::rows &&
                          detail::traits< std::remove_cvref_t<A> >::cols == detail::traits< std::remove_cvref_t<B> >::cols);

    // InnerSizesMatch: A can be multiplied by B (columns of A = rows of B)
    template<typename A, typename B>
    concept InnerSizesMatch = !detail::both_fixed_size_v< std::remove_cvref_t<A>, std::remove_cvref_t<B> > ||
                              detail::traits< std::remove_cvref_t<A> >::cols == detail::traits< std::remove_cvref_t<B> >::rows;
}

#endif
