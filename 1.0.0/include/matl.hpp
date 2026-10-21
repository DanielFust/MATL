#ifndef MATL_HPP
#define MATL_HPP

// Debug Flags (enabled by default; define MATL_NO_DEBUG before including matl.hpp to disable the
// dimension and bounds checks and the MATL_DEBUG trace output, e.g. for benchmarks)
#ifndef MATL_NO_DEBUG
#define MATL_DIMENSION_CHECK
#define MATL_BOUNDS_CHECK
#define MATL_DEBUG
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
#include <cstddef>

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

#endif