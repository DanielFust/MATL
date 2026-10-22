/**************************************************************************************************
 * File: _matl_errors_.hpp
 *
 * Description: Error reporting.
 *
 * Every error MATL detects is reported through detail::raise, so the policy is set in one place:
 *
 *   - By default, MATL throws an exception:
 *       matl::dimension_error  (derived from std::invalid_argument)  incompatible dimensions
 *       matl::index_error      (derived from std::out_of_range)      element index out of range
 *       std::logic_error                                             misuse of the TemporaryArena
 *   - If MATL_ABORT_ON_ERROR is defined, or exceptions are disabled (e.g. -fno-exceptions), MATL
 *     prints the message to stderr and calls std::abort(), which stops at the faulty statement
 *     under a debugger.
 *
 * MATL_THROWS_ON_ERROR is 1 if errors throw and 0 if they abort.
 *
 * Most checks are debug checks (MATL_BOUNDS_CHECK, MATL_DIMENSION_CHECK, see matl.hpp); they
 * report programming errors, so a program would normally not catch them except to log context or
 * shut down cleanly. A few checks are always on (e.g. allocate() of a fixed-size matrix with other
 * dimensions, or an invalid TemporaryArena::release()).
 *
 * An error detected inside a noexcept function (the move assignment of a matrix over external
 * memory) cannot propagate; the program terminates with the error message instead.
 **************************************************************************************************/
#ifndef MATL_ERRORS_HPP
#define MATL_ERRORS_HPP

#if !defined(MATL_ABORT_ON_ERROR) && (defined(__cpp_exceptions) || defined(__EXCEPTIONS) || defined(_CPPUNWIND))
#define MATL_THROWS_ON_ERROR 1
#else
#define MATL_THROWS_ON_ERROR 0
#endif

namespace matl
{
   /***********************************************************************************************
    * class dimension_error
    *
    * Thrown when the dimensions of operands, or of an assignment target, do not fit.
    ***********************************************************************************************/
    class dimension_error : public std::invalid_argument
    {
      public:
        using std::invalid_argument::invalid_argument;
    };

   /***********************************************************************************************
    * class index_error
    *
    * Thrown when an element index is out of range (indices are 1-based).
    ***********************************************************************************************/
    class index_error : public std::out_of_range
    {
      public:
        using std::out_of_range::out_of_range;
    };

    namespace detail
    {
       /*******************************************************************************************
        * void raise<EXCEPTION>(message)
        *
        * Reports an error: throws EXCEPTION(message), or prints the message and aborts (see the
        * top of this file).
        *******************************************************************************************/
        template<typename EXCEPTION>
        [[noreturn]] void raise(std::string const& message)
        {
            #if MATL_THROWS_ON_ERROR
            throw EXCEPTION(message);
            #else
            std::cerr << "MATL error: " << message << std::endl;
            std::abort();
            #endif
        }

        // "MxN"
        inline std::string dimensions(index_t rows, index_t cols)
        {return std::to_string(rows) + "x" + std::to_string(cols);}

        // Operands of an operation have incompatible dimensions
        [[noreturn]] inline void incompatible_dimensions(char const* operation,
                                                         index_t rows1, index_t cols1,
                                                         index_t rows2, index_t cols2)
        {
            raise<dimension_error>(std::string("incompatible matrix dimensions in ") + operation + ": " +
                                   dimensions(rows1, cols1) + " and " + dimensions(rows2, cols2));
        }

        // A value of one size is assigned (or copied) to a matrix of another size
        [[noreturn]] inline void assignment_dimensions(index_t target_rows, index_t target_cols,
                                                       index_t rows, index_t cols)
        {
            raise<dimension_error>("cannot assign a " + dimensions(rows, cols) + " value to a " +
                                   dimensions(target_rows, target_cols) + " matrix");
        }

        // Single (storage-order) index out of range
        [[noreturn]] inline void index_out_of_range(index_t i, index_t size)
        {
            raise<index_error>("index " + std::to_string(i) + " is out of range for a matrix with " +
                               std::to_string(size) + " elements (indices are 1-based)");
        }

        // Row and column index out of range
        [[noreturn]] inline void index_out_of_range(index_t i, index_t j, index_t rows, index_t cols)
        {
            raise<index_error>("index (" + std::to_string(i) + "," + std::to_string(j) +
                               ") is out of range for a " + dimensions(rows, cols) +
                               " matrix (indices are 1-based)");
        }
    }
}
#endif
