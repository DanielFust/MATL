# Changelog

Notable changes to MATL, newest first. Version numbers follow [semantic versioning](https://semver.org):
a new major version may require changes to code written for the previous one.

## [Unreleased]

### Added

- **Sparse matrices:** `CSRMatrix<T, I>` (compressed sparse row, index type `I`), assembled with
  `SparseBuilder` (entries in any order, duplicates summed) or from the CSR arrays (0-based, as in
  SciPy, MKL and PETSc); element access, a fixed pattern with changeable values, products with
  dense matrices and vectors (`S * X`, `transpose(S) * X`, as expressions), scaling, `dense` and
  `sparse` conversions, and output. Dense factorizations reject sparse matrices with a compile-time
  message. Example 7 assembles 1D finite-element matrices.

## [2.1.0] - 2026-10-22

Adds the QR factorization and least squares, and makes the factorizations faster. Code written for
2.0.0 compiles unchanged.

### Added

- **QR factorization and least squares:** `qr(A)` (Householder, blocked) for matrices of any
  shape, with `solve` (least squares for m ≥ n), `Q()`, `full_Q()`, `R()`, `Qt_times`, `Q_times`
  and `rcond`; `qr_solve(A, B)` and `qr_inplace(A)`. Example 6 now fits with QR and compares the
  normal equations.

### Changed

- **Faster factorizations:** LU (Toledo) and Cholesky are now recursive, and QR factorizes its
  panels recursively (as LAPACK's `xGEQRT3`), so that almost all their arithmetic is matrix
  products. Block updates subtract the product inside the product kernel, without a temporary.
  At n = 1024: LU about 1.2–1.8× faster (15–17 GFLOP/s), Cholesky about 1.2–1.5× (11–12 GFLOP/s),
  QR about 1.3× (13–15 GFLOP/s), and QR of tall matrices about 1.6–1.9× (4000 x 400: 9–11 GFLOP/s).
- `norm_1` and `norm_inf` visit the elements in storage order (one running sum per column or row),
  which is faster for row-major `norm_1` and column-major `norm_inf`.

## [2.0.0] - 2026-10-22

Version 2.0.0 rewrites MATL's internals (expression evaluation, memory management, products) and
adds most of the library's features. Code written for 1.0.0 mostly compiles unchanged; see
[Breaking changes](#breaking-changes) for what does not.

### Added

- **Creating and printing matrices:** initializer lists (`Matrix<double> A = {{1, 2}, {3, 4}};`),
  `zeros` and `ones`, `operator<<` (aligned columns, using the stream's formatting), and `print`
  of expressions.
- **Elementwise operations:** negation; `abs`, `sqrt`, `cbrt`, `exp`, `log`, `log2`, `log10`, the
  trigonometric and hyperbolic functions, `floor`, `ceil` and `round`; `times`, `rdivide` and
  `pow` (MATLAB's `.*`, `./` and `.^`), with a scalar on either side.
- **Reductions:** `sum`, `prod`, `mean`, `min`, `max`, `trace`, `dot`, `norm_fro`, `norm_1` and
  `norm_inf`. Expressions are reduced without temporary matrices.
- **Column-major matrices:** `ColMajorMatrix<T>` and `FixedSizeColMajorMatrix<T, M, N>`, which mix
  freely with row-major matrices. They wrap Fortran and LAPACK arrays without copying.
- **Views:** `A(range(1, 3), all)`, `A(all, 2)`, `A(range(1, 2, last), all)`,
  `A(range(last, -1, 1), all)` select rows and columns without copying; assigning to a view writes
  the matrix. Element access with `last` (`A(last, 1)`), `swap` of two views, and `data()`.
- **Linear systems:** `lu` (partial pivoting) and `cholesky` factorizations, with `solve`,
  `solve_inplace`, `rcond` (condition estimate), `determinant` and the factors; `lu_solve`,
  `cholesky_solve`, `lu_inplace` and `cholesky_inplace`; triangular solves `solve_lower` and
  `solve_upper`.
- **Error reporting:** the exceptions `matl::dimension_error`, `matl::index_error`,
  `matl::singular_matrix_error` and `matl::not_positive_definite_error`. Define
  `MATL_ABORT_ON_ERROR` to print the message and abort instead.
- **Dimension checks:** fixed-size operands of the wrong dimensions do not compile, and every
  operation in an expression is checked at run time with `MATL_DIMENSION_CHECK`.
- **Examples:** six annotated programs in `examples/`, from the basics to linear systems.
- **Build and testing:** a CMake package (`find_package(matl)`, target `matl::matl`), a test suite
  (`make check`, or CTest), benchmarks, and CI on GCC 10–14, Clang 13–18 and Apple Clang.
- **Project:** the BSD-3-Clause license, and `CITATION.cff` for citing MATL.

### Changed

- **Expressions are formulas.** Building an expression computes nothing; assigning it to a matrix
  evaluates it in one pass. Stored expressions (`auto e = A + B;`) are safe and read their
  operands when evaluated, and temporary matrices inside expressions are kept alive.
- **Performance.** Matrix products use a cache-blocked kernel with SIMD vector extensions (10–20×
  faster for large matrices), transposes are tiled, and the temporaries of an evaluation come
  from a recycled per-thread arena, so repeated evaluations make no heap allocations.
- **Debug checks follow `NDEBUG`,** like `assert()`: `MATL_BOUNDS_CHECK` and `MATL_DIMENSION_CHECK`
  are on unless `NDEBUG` or `MATL_NO_DEBUG` is defined. The `MATL_DEBUG` trace output is opt-in.
- **Errors throw exceptions** instead of printing a message and aborting.
- **Types:** C++20 concepts (`MatrixType`, `ExpressionType`, `Scalar`, ...) replace `enable_if`;
  expression types are in `matl::expr`, and internal ones in `matl::detail`. `Scalar` means
  arithmetic types only.
- **Repository layout:** the library is at the repository root (`include/matl.hpp`) instead of
  under `1.0.0/`; versions are git tags.

### Breaking changes

| 1.0.0 | 2.0.0 |
| --- | --- |
| Scratch-memory allocation: `add_scratch_memory`, `ALLOCATE_MATRICES_IN_SCRATCH`, `ScratchContainer`, `MemoryManager::allocate` | Removed. For memory you own, use the external-memory constructor `Matrix<double> A(m, n, pointer)`; otherwise `allocate` and `deallocate` |
| Element access on an expression object, `(A + B)(1, 1)` | Assign the expression to a matrix first |
| Expression node types `_MATRIX_TIMES_MATRIX_`, `_MATRIX_TRANSPOSE_`, `_OWNED_MATRIX_` and the other `_..._` names | `expr::Product`, `expr::Transpose`, `expr::Owned`, `expr::Elementwise`; use `auto` rather than naming them |
| Interface parameters `_MATRIX_<MAT>&`, `_EXPRESSION_<EXPR> const&` | The concepts `MatrixType`, `ExpressionType` and `MatrixAlgebraObject` |
| Failed checks print a message and abort | They throw; define `MATL_ABORT_ON_ERROR` for the old behavior |
| Debug checks on in every build | Off when `NDEBUG` (or `MATL_NO_DEBUG`) is defined |
| The library under `1.0.0/` in the repository | At the repository root: put `include/` on the include path and `#include <matl.hpp>` |

### Fixed

- Wrong results for `A ± transpose(B)`, and for some expressions that mixed products and sums
  (e.g. `I*D + Z`).
- Stored expressions read freed memory once the statement that built them ended.
- Copying skipped the last row and column, and the generic copy constructor broke non-square
  matrices.
- Integer `/=` by a scalar gave 0.
- The dimension checks never ran (they tested a misspelled macro).
- Memory bugs: double destruction of temporaries, a leak when re-allocating a matrix, and
  several in the scratch allocator (since removed).
- `matl.hpp` could not be included from more than one translation unit.
- Assignment from `const` matrices did not compile, and scalar assignment returned `void`.

## [1.0.0] - 2025-04-07

The original version, kept as the git tag `v1.0.0-backup`.

[Unreleased]: https://github.com/DanielFust/MATL/compare/v2.1.0...HEAD
[2.1.0]: https://github.com/DanielFust/MATL/releases/tag/v2.1.0
[2.0.0]: https://github.com/DanielFust/MATL/releases/tag/v2.0.0
[1.0.0]: https://github.com/DanielFust/MATL/tree/v1.0.0-backup
