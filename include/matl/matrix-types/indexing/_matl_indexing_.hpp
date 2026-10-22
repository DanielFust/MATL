/**************************************************************************************************
 * File: _matl_indexing_.hpp
 *
 * Description: Indices that select the rows and columns of a view (see MatrixView):
 *
 *     range(first, last)          first, first+1, ..., last              MATLAB first:last
 *     range(first, step, last)    first, first+step, ..., up to last     MATLAB first:step:last
 *     all                         every row or column                    MATLAB :
 *     last, last - k              the last row or column, k before it    MATLAB end, end-k
 *
 * Indices are 1-based and a range includes its last index, as in MATLAB. The step may be
 * negative: range(last, -1, 1) runs backwards. A range whose step points away from its last index
 * is empty (e.g. range(3, 2)), as in MATLAB. A step of 0 is always reported (matl::index_error);
 * indices outside the matrix are reported when MATL_BOUNDS_CHECK is defined.
 **************************************************************************************************/
#ifndef MATL_INDEXING_HPP
#define MATL_INDEXING_HPP
namespace matl
{
    // Type of 'all'
    struct AllIndices {};

    // Every row or column (MATLAB ':')
    inline constexpr AllIndices all{};

    // Type of 'last' and 'last - k': the index 'offset' places before the last row or column
    struct LastIndex
    {
        std::ptrdiff_t offset = 0;
    };

    // The last row or column (MATLAB 'end')
    inline constexpr LastIndex last{};

    // last - k
    template<std::integral I>
    constexpr LastIndex operator- (LastIndex l, I k)
    {return LastIndex{l.offset + static_cast<std::ptrdiff_t>(k)};}

    // An end point of a range: an index, or a position counted back from the last index
    struct IndexBound
    {
        std::ptrdiff_t value     = 0;       // the index, or the offset from the last index
        bool           from_last = false;

        template<std::integral I>
        constexpr IndexBound(I i) : value(static_cast<std::ptrdiff_t>(i)), from_last(false) {}

        constexpr IndexBound(LastIndex l) : value(l.offset), from_last(true) {}

        // The index in a dimension of size n
        constexpr std::ptrdiff_t resolve(index_t n) const
        {return this->from_last ? static_cast<std::ptrdiff_t>(n) - this->value : this->value;}
    };

    // Type of range(...): from, from + step, ..., up to 'to'
    struct IndexRange
    {
        IndexBound     from;
        std::ptrdiff_t step;
        IndexBound     to;
    };

    // first, first + 1, ..., last (MATLAB first:last)
    inline IndexRange range(IndexBound first, IndexBound last)
    {return IndexRange{first, 1, last};}

    // first, first + step, ..., up to last (MATLAB first:step:last)
    template<std::integral S>
    IndexRange range(IndexBound first, S step, IndexBound last)
    {return IndexRange{first, static_cast<std::ptrdiff_t>(step), last};}


    namespace detail
    {
        // One index argument of a view: an integer, last (- k), all, or a range
        template<typename X>
        concept ViewIndex = (std::integral<X> && !std::same_as<X, bool>) || std::same_as<X, LastIndex> ||
                            std::same_as<X, AllIndices> || std::same_as<X, IndexRange>;

        // Selects several indices (all or a range)
        template<typename X>
        concept MultipleIndices = std::same_as<X, AllIndices> || std::same_as<X, IndexRange>;

        // Index arguments that select a view: at least one of them is all or a range
        // (A(1, all) is a view of row 1; A(1, 2) is an element)
        template<typename R, typename C>
        concept ViewIndices = ViewIndex<std::remove_cvref_t<R>> && ViewIndex<std::remove_cvref_t<C>> &&
                              (MultipleIndices<std::remove_cvref_t<R>> || MultipleIndices<std::remove_cvref_t<C>>);

        // Index arguments that select one element and use last: A(last, 1), A(2, last - 1)
        template<typename R, typename C>
        concept ElementIndices = ViewIndex<std::remove_cvref_t<R>> && ViewIndex<std::remove_cvref_t<C>> &&
                                 !MultipleIndices<std::remove_cvref_t<R>> && !MultipleIndices<std::remove_cvref_t<C>> &&
                                 (std::same_as<std::remove_cvref_t<R>, LastIndex> ||
                                  std::same_as<std::remove_cvref_t<C>, LastIndex>);

        // The index of an integer or last (- k) in a dimension of size n (an index outside the
        // dimension is reported by the element access)
        template<typename X>
        index_t element_index(X const& x, index_t n)
        {
            if constexpr (std::same_as<X, LastIndex>) return static_cast<index_t>(IndexBound(x).resolve(n));
            else                                      return static_cast<index_t>(x);
        }

        // The indices selected in one dimension: first, first + step, ... (count indices)
        struct Selection
        {
            std::ptrdiff_t first;
            std::ptrdiff_t step;
            index_t        count;
        };

        // Indices selected by an index argument in a dimension of size n. 'dimension' ("row" or
        // "column") names the dimension in error messages.
        template<typename X>
        Selection select_indices(X const& x, index_t n, [[maybe_unused]] char const* dimension)
        {
            Selection s{1, 1, 0};
            if constexpr (std::same_as<X, AllIndices>)
                s = Selection{1, 1, n};
            else if constexpr (std::same_as<X, IndexRange>)
            {
                if (x.step == 0) raise<index_error>("the step of a range must not be 0");
                const std::ptrdiff_t from = x.from.resolve(n);
                const std::ptrdiff_t to   = x.to.resolve(n);
                index_t count = 0;
                if (x.step > 0 && to >= from) count = static_cast<index_t>((to - from) / x.step) + 1;
                if (x.step < 0 && from >= to) count = static_cast<index_t>((from - to) / -x.step) + 1;
                s = Selection{from, x.step, count};
            }
            else
                s = Selection{static_cast<std::ptrdiff_t>(element_index(x, n)), 1, 1};

            #ifdef MATL_BOUNDS_CHECK
            if (s.count > 0)
            {
                // The first and the last selected index must be in the dimension (an empty
                // selection selects nothing, so it may start anywhere, e.g. range(last + 1, last))
                const std::ptrdiff_t size       = static_cast<std::ptrdiff_t>(n);
                const std::ptrdiff_t last_index = s.first + static_cast<std::ptrdiff_t>(s.count - 1) * s.step;
                const auto outside = [size](std::ptrdiff_t k) {return k < 1 || k > size;};
                if (outside(s.first) || outside(last_index))
                    raise<index_error>(std::string(dimension) + " index " +
                                       std::to_string(outside(s.first) ? s.first : last_index) +
                                       " is out of range for a matrix with " + std::to_string(n) + " " +
                                       dimension + "s (indices are 1-based)");
            }
            #endif
            return s;
        }
    }
}
#endif
