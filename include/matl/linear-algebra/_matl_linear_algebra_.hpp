/**************************************************************************************************
 * File: _matl_linear_algebra_.hpp
 *
 * Description: Linear systems with dense matrices: triangular solves, the LU and Cholesky
 *              factorizations, and the QR factorization with least squares. Included after the API functions, which they use (norm_1, eye).
 **************************************************************************************************/
#include "./triangular/_matl_triangular_.hpp"
#include "./lu/_matl_lu_.hpp"
#include "./cholesky/_matl_cholesky_.hpp"
#include "./qr/_matl_qr_.hpp"
