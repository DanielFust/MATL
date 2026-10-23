#ifndef MATL_HPP
#define MATL_HPP

// Version (MATL_VERSION is MAJOR * 10000 + MINOR * 100 + PATCH, e.g. 20200 for 2.2.0, so it can be
// compared in #if). CMakeLists.txt reads the version from these lines.
#define MATL_VERSION_MAJOR 2
#define MATL_VERSION_MINOR 2
#define MATL_VERSION_PATCH 0
#define MATL_VERSION (MATL_VERSION_MAJOR * 10000 + MATL_VERSION_MINOR * 100 + MATL_VERSION_PATCH)

// Debug checks
//
//   MATL_BOUNDS_CHECK       Index checks on element access
//   MATL_DIMENSION_CHECK    Run-time dimension checks of assignments and expressions
//
// Like assert(), the checks are enabled unless NDEBUG is defined (as in typical release builds).
// Define MATL_NO_DEBUG to disable them in any build, or define a check's macro before including
// matl.hpp to keep it in a release build. Dimensions known at compile time are always checked.
//
// Errors throw matl::dimension_error or matl::index_error; define MATL_ABORT_ON_ERROR to print the
// message and abort instead (see matl/errors/_matl_errors_.hpp).
//
// MATL_DEBUG (off by default) prints a trace of internal calls (assignments, allocations) to
// stdout. Define it before including matl.hpp to enable it.
#if !defined(NDEBUG) && !defined(MATL_NO_DEBUG)
    #ifndef MATL_DIMENSION_CHECK
    #define MATL_DIMENSION_CHECK
    #endif
    #ifndef MATL_BOUNDS_CHECK
    #define MATL_BOUNDS_CHECK
    #endif
#endif

// Standard Library dependencies
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <complex>
#include <vector>
#include <memory>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <initializer_list>
#include <limits>
#include <span>
#include <sstream>
#include <string>

// Namespace type definitions and settings
namespace matl
{
    // Index Type definitions
    typedef unsigned long index_t;

    // Default settings for MATL library. Settings for memory management can
    // be changed dynamically through the MemoryManager class (enable/disable and
    // CACHE_TEMPORARY_MATRICES).
    namespace defaults
    {
        // Default settings for memory manager (allocate evaluation temporaries from the
        // per-thread arena)
        const bool    ENABLE_MEMORY_MANAGER            = true;
    }

    // Set block size for cache blocking
    //const size_t BLOCK_SIZE = 1024;
}

// Error reporting (exceptions, or abort with MATL_ABORT_ON_ERROR or without exceptions)
#include "./matl/errors/_matl_errors_.hpp"

// Thread parallelism with OpenMP (when compiled with -fopenmp; see matl/parallel/_matl_parallel_.hpp)
#include "./matl/parallel/_matl_parallel_.hpp"

// Forward class and operator declarations
#include "./matl/_matl_declarations_.hpp" 

// Concepts and type traits
/*(Type traits MUST be included BEFORE class and function prototypes but AFTER class declarations) */
#include "./matl/_matl_type_traits_.hpp"

// Class Prototypes
#include "./matl/_matl_prototypes_.hpp"

// MemoryManager instance (memory_manager is a threadlocal singleton). The singleton pointer and
// static settings are defined as inline variables in the MemoryManager class prototype so that
// matl.hpp can be included in multiple translation units.

// Operators
#include "./matl/operators/_matl_operators_.hpp"

// Implementations
#include "./matl/_matl_implementations_.hpp"

// Additional API
#include "./matl/_matl_api_functions_.hpp"

// Linear systems (triangular solves, LU, Cholesky and QR factorizations, least squares)
#include "./matl/linear-algebra/_matl_linear_algebra_.hpp"

// Sparse matrices (CSR format, assembly, products with dense matrices)
#include "./matl/sparse/_matl_sparse_.hpp"

#endif