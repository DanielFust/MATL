/**************************************************************************************************
 * File: _matl_sparse_.hpp
 *
 * Description: Sparse matrices: the CSR format, assembly with SparseBuilder, products with dense
 *              matrices and vectors, and conversions. Included after the dense library, whose
 *              evaluation machinery the sparse products use.
 **************************************************************************************************/
#include "./csr/_matl_csr_matrix_.hpp"
#include "./builder/_matl_sparse_builder_.hpp"
#include "./product/_matl_sparse_product_.hpp"
#include "./functions/_matl_sparse_functions_.hpp"
