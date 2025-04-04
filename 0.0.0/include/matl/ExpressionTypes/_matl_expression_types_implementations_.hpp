#ifndef _MATL_EXPRESSION_TYPES_IMPLEMENTATIONS__HPP
#define _MATL_EXPRESSION_TYPES_IMPLEMENTATIONS__HPP

/**************************************************************************************************
 * Expression Template Supertypes
 * 
 * Polymorphic class interfaces for expression templates
 **************************************************************************************************/
 // Matrix Expression (Base class for all Matrix Expression templates)
#include "./MatExpression/_matl_matexpression_implementation_.hpp"

/**************************************************************************************************
 * Expression Templates
 **************************************************************************************************/
// Addition Expression Templates
#include "./Addition/_matl_addition_expressions_implementations_.hpp"
// Subtraction Expression Templates
#include "./Subtraction/_matl_subtraction_expressions_implementations_.hpp"
// Multiplication Expression Templates
#include "./Multiplication/_matl_multiplication_expressions_implementations_.hpp"
// Transposition Expression Templates
#include "./Transpose/_matl_transpose_expressions_implementations_.hpp"



#endif