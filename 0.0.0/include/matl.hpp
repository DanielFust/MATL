#ifndef MATL_HPP
#define MATL_HPP

// Debug Flags
#ifndef MATL_DIM_CHECK
#define MATL_DIM_CHECK
#endif
#ifndef MATL_BOUNDS_CHECK
#define MATL_BOUNDS_CHECK
#endif

// Namespace type definitions
namespace matl
{
    typedef unsigned long index_t;
}

// Forward class and operator declarations
#include "./matl/_matl_forward_declarations_.hpp" 

// Utilities & Type Traits 
/*(Type traits MUST be included BEFORE class and function prototypes but AFTER class declarations) */
#include "./matl/_matl_utilities_.hpp"

// Class Prototypes
#include "./matl/_matl_prototypes_.hpp"

// Implementations
#include "./matl/_matl_implementations_.hpp"


#endif