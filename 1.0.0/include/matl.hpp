#ifndef MATL_HPP
#define MATL_HPP

// Debug Flags
#define MATL_DIMENSION_CHECK
#define MATL_BOUNDS_CHECK
#define MATL_DEBUG

// Standard Library dependencies
#include <type_traits>
#include <iostream>
#include <stdexcept>
#include <complex>
#include <vector>
#include <deque>
#include <memory>
#include <cstddef>

// Namespace type definitions and settings
namespace matl
{
    // Index Type definitions
    typedef unsigned long index_t;

    // Default settings for MATL library. Settings for memory management can
    // be changed dynamically by setting the static variables of the same
    // names in the MemoryManager class.
    namespace defaults
    {
        // Default settings for memory manager
        const bool    ENABLE_MEMORY_MANAGER            = true;
        const index_t INITIAL_SCRATCH_MEM_SIZE         = 1024;
        const index_t DEFAULT_NUM_SCRATCH_CONTAINERS   = 10;
        const index_t SCRATCH_BLOCK_SIZE               = 1024;
        const bool    ALLOCATE_MATRICES_IN_SCRATCH = false;
    }

    // Set block size for cache blocking
    //const size_t BLOCK_SIZE = 1024;
}

// Forward class and operator declarations
#include "./matl/_matl_declarations_.hpp" 

// Type Traits 
/*(Type traits MUST be included BEFORE class and function prototypes but AFTER class declarations) */
#include "./matl/_matl_type_traits_.hpp"

// Class Prototypes
#include "./matl/_matl_prototypes_.hpp"

// Utilities (Must be included after Type Traits but before implementations)
#include "./matl/utilities/_matl_utilities_.hpp"

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