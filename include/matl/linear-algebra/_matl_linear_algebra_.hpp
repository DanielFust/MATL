/**************************************************************************************************
 * File: _matl_linear_algebra_.hpp
 *
 * Description: Dense linear algebra: triangular solves, the LU and Cholesky factorizations, the
 *              QR factorization with least squares, determinants and inverses, eigenvalues of
 *              symmetric matrices, and the singular value decomposition with the pseudo-inverse.
 *              Included after the API functions, which they use (norm_1, eye).
 **************************************************************************************************/
#include "./triangular/_matl_triangular_.hpp"
#include "./lu/_matl_lu_.hpp"
#include "./cholesky/_matl_cholesky_.hpp"
#include "./qr/_matl_qr_.hpp"
#include "./inverse/_matl_inverse_.hpp"
#include "./eigen/_matl_eigen_.hpp"
#include "./svd/_matl_svd_.hpp"
