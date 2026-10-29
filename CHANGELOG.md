# Changelog

Notable changes to MATL, newest first. Version numbers follow [semantic versioning](https://semver.org):
a new major version may require changes to code written for the previous one.

## [2.3.0] - 2026-10-29

Adds multithreading with OpenMP, parallel finite-element assembly, fixed-size factorizations
without heap memory, determinants and inverses, symmetric eigenvalues, and the singular value
decomposition with the pseudo-inverse. Most code written for 2.2.0 compiles unchanged; see
Changed for fixed-size solutions and the rounding of large sums.

### Added

- **Threads (OpenMP):** compiled with `-fopenmp`, MATL runs elementwise evaluation, reductions,
  sparse matrix times dense products and the Jacobi preconditioner on several threads (above about
  65,000 elements of work, and never nested in a caller's parallel region). `MATL_NO_OPENMP` keeps
  MATL serial. Exceptions inside parallel loops are rethrown after them. `make check-openmp` and
  CTest (when OpenMP is found) run the test suite on 4 threads.
- **Parallel dense products:** the product kernel divides the row blocks (or, for products with
  few rows, the column panels) of large products (from about 2 million multiply-adds) among the
  threads; results are the same on any number of threads. The factorizations' block updates use
  it (1536 x 1536 on 8 threads: product 5.4x, LU 2.0x, Cholesky 1.5x, QR 2.4x).
- **`SerialScope`:** keeps MATL serial in the current thread (for threads of other frameworks,
  e.g. `std::thread` or TBB, which OpenMP cannot detect).
- **Parallel assembly:** `SparseBuilder::build(builders)` merges per-thread builders; with
  `schedule(static)` the result is the same as a serial assembly's.
- **Fixed-size factorizations:** `lu(A)` and `cholesky(A)` of a fixed-size matrix or expression
  return `FixedLUFactorization<T, N>` and `FixedCholeskyFactorization<T, N>`, which keep their
  factors in the object (pivots in a `std::array`) and factorize and solve without heap memory.
- **Determinants and inverses:** `det(A)` (0 for a singular matrix) and `inv(A)`, with closed
  forms for fixed-size matrices of order 1 to 3 (e.g. element Jacobians) and LU otherwise.
- **Symmetric eigenvalues:** `eig_sym(A)` (eigenvalues in ascending order and orthonormal
  eigenvectors, by cyclic Jacobi) and `eigenvalues_sym(A)`.
- **Singular value decomposition:** `svd(A)` for any shape (one-sided Jacobi), with `U()`,
  `singular_values()`, `V()`, `S()`, `rank()`, `cond()`, `pinv()` and the minimum-norm
  least-squares `solve`; `singular_values(A)`, `pinv(A)` and `rank(A)`, with an optional
  tolerance. Example 10 shows all of these on problems from mechanics and data fitting.
- Fixed-size matrices give fixed-size results without heap memory for `det`, `inv`, `eig_sym` and
  `svd`.
- **Fixed-size matrices from dynamic ones:** `Matrix<double, 4, 1> f = lu(K).solve(b);` copies a
  dynamic matrix into a fixed-size one (the dimensions are checked at run time).

### Changed

- **Reductions are reproducible:** `sum`, `prod`, `mean`, `min`, `max`, `dot` and `norm_fro` add
  their elements in fixed blocks of 4096, combined in order, so results are the same for any number
  of threads. For more than 4096 elements, sums round slightly differently than in 2.2.0. Large
  reductions keep their partial results in the temporary arena (no heap allocations).
- **Fixed-size solutions:** `solve` of LU and Cholesky factorizations, `lu_solve`, `cholesky_solve`,
  `solve_lower` and `solve_upper` return a fixed-size matrix when the right-hand side is fixed-size
  (before, a dynamic one). `norm_1` and `norm_inf` keep their column or row sums in the temporary
  arena instead of a `std::vector`.

## [2.2.0] - 2026-10-23

Adds sparse matrices: storage and assembly, products with dense matrices, iterative solvers with
preconditioners, sparse arithmetic, and direct solvers for tridiagonal and banded matrices. Code
written for 2.1.0 compiles unchanged.

### Added

- **Sparse matrices:** `CSRMatrix<T, I>` (compressed sparse row, index type `I`), assembled with
  `SparseBuilder` (entries in any order, duplicates summed) or from the CSR arrays (0-based, as in
  SciPy, MKL and PETSc); element access, a fixed pattern with changeable values, products with
  dense matrices and vectors (`S * X`, `transpose(S) * X`, as expressions), scaling, `dense` and
  `sparse` conversions, and output. Dense factorizations reject sparse matrices with a compile-time
  message. Example 7 assembles 1D finite-element matrices.
- **Iterative solvers:** `cg` (symmetric positive definite), `bicgstab` and restarted `gmres`
  (general), for sparse matrices, dense matrices, or operator functions `op(in, out)`; options for
  the tolerance, iteration limit and GMRES restart, and a result reporting convergence, iterations
  and the residual. Preconditioners: `JacobiPreconditioner` and `ILU0Preconditioner`. Example 8
  solves the 2D Poisson equation.
- **Sparse arithmetic:** `S + T`, `S - T`, `S * T` (Gustavson's algorithm), the explicit transpose
  (`CSRMatrix<double> St = transpose(S);`), and `diagonal(S)`.
- **Structured direct solvers:** `solve_tridiagonal` (Thomas algorithm), and `BandedMatrix` with
  banded LU (partial pivoting) and Cholesky factorizations (`lu(B)`, `cholesky(B)`, with `solve`,
  `rcond` and `determinant`), `banded(S)` conversion, and products. Example 9 steps the 1D heat
  equation and solves the 2D Poisson equation directly.

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

[2.3.0]: https://github.com/DanielFust/MATL/releases/tag/v2.3.0
[2.2.0]: https://github.com/DanielFust/MATL/releases/tag/v2.2.0
[2.1.0]: https://github.com/DanielFust/MATL/releases/tag/v2.1.0
[2.0.0]: https://github.com/DanielFust/MATL/releases/tag/v2.0.0
[1.0.0]: https://github.com/DanielFust/MATL/tree/v1.0.0-backup
