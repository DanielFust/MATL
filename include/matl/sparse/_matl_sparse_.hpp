/**************************************************************************************************
 * File: _matl_sparse_.hpp
 *
 * Description: Sparse matrices: the CSR format, assembly with SparseBuilder, products with dense
 *              matrices and vectors, conversions, sparse arithmetic, tridiagonal and banded
 *              direct solvers, and the iterative solvers with their preconditioners. Included after the dense library, whose
 *              evaluation machinery the sparse products use.
 **************************************************************************************************/
#include "./csr/_matl_csr_matrix_.hpp"
#include "./builder/_matl_sparse_builder_.hpp"
#include "./product/_matl_sparse_product_.hpp"
#include "./functions/_matl_sparse_functions_.hpp"
#include "./arithmetic/_matl_sparse_arithmetic_.hpp"
#include "./banded/_matl_banded_.hpp"
#include "./solvers/_matl_preconditioners_.hpp"
#include "./solvers/_matl_iterative_solvers_.hpp"
