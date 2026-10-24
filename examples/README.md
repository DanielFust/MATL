# MATL examples

Short, annotated programs that show how to use MATL. Each one builds on the previous ones, so
they are best read in order:

| Example | Shows |
| --- | --- |
| [`01_basics.cpp`](01_basics.cpp) | Creating matrices (initializer lists, `zeros`, `ones`, `eye`, fixed sizes), 1-based element access and `last`, arithmetic, matrix products, `transpose`, elementwise functions (`times`, `rdivide`, `pow`, `sqrt`, ...), reductions (`sum`, `max`, norms, ...), and what `auto` does with an expression |
| [`02_views.cpp`](02_views.cpp) | Selecting rows, columns and blocks with `range`, `all` and `last`; writing through views; row operations and `swap`; views versus copies; Gaussian elimination written with views |
| [`03_linear_systems.cpp`](03_linear_systems.cpp) | `lu` and `cholesky` factorizations, solving for one or many right-hand sides, residuals, condition estimates (`rcond`), determinants, triangular solves, in-place factorization, singular matrices |
| [`04_storage_and_memory.cpp`](04_storage_and_memory.cpp) | Row-major and column-major matrices, matrices over memory you own (a buffer, a Fortran array), passing elements to C or Fortran code with `data()`, fixed-size matrices |
| [`05_errors_and_checks.cpp`](05_errors_and_checks.cpp) | Compile-time dimension checks, the debug checks (`MATL_DIMENSION_CHECK`, `MATL_BOUNDS_CHECK`) and when they run, the checks that always run, and catching MATL's exceptions |
| [`06_polynomial_fit.cpp`](06_polynomial_fit.cpp) | An application: a least-squares polynomial fit with QR, compared with the normal equations (Cholesky), combining views, elementwise powers, products, condition estimates and reductions |
| [`07_sparse_matrices.cpp`](07_sparse_matrices.cpp) | Sparse matrices: assembling 1D finite-element stiffness and mass matrices element by element, products with vectors, changing values in a fixed pattern, the CSR arrays, 32-bit indices |
| [`08_iterative_solvers.cpp`](08_iterative_solvers.cpp) | Iterative solvers: the 2D Poisson equation with CG, Jacobi and ILU(0) preconditioning, tolerances and iteration limits, a nonsymmetric problem with BiCGSTAB and GMRES, and a matrix-free stencil operator |
| [`09_direct_sparse_solvers.cpp`](09_direct_sparse_solvers.cpp) | Direct solvers for structured sparse matrices: the 1D heat equation with implicit time steps (tridiagonal), sparse arithmetic to build a 2D Laplacian, and the 2D Poisson equation by banded Cholesky, compared with CG |
| [`10_eigenvalues_and_svd.cpp`](10_eigenvalues_and_svd.cpp) | `det` and `inv` of an element Jacobian, principal stresses with `eig_sym`, natural frequencies of a spring-mass chain, and a rank-deficient least-squares fit with `svd`, `rank` and `pinv` |

## Building and running

From the repository root, with Make:

```sh
make examples                     # builds build/examples/01_basics, ...
./build/examples/01_basics
make run-examples                 # builds and runs all of them (as CI does)
make examples CXX=clang++         # choose a compiler (default: g++)
```

With CMake:

```sh
cmake -S . -B build/cmake -DMATL_BUILD_EXAMPLES=ON
cmake --build build/cmake
./build/cmake/matl_example_01_basics
```

Or compile one directly. MATL is header-only, so only the include path is needed:

```sh
g++ -std=c++20 -I include examples/01_basics.cpp -o basics
```

MATL's debug checks are on unless `NDEBUG` is defined. The Make build keeps them on, while CMake's
default Release build turns them off; `05_errors_and_checks.cpp` shows the difference.

## Using MATL in your own project

Add MATL's `include/` directory to your include path and compile with `-std=c++20`, or use the
`matl::matl` CMake target (see [Using MATL in your project](../README.md#using-matl-in-your-project)).
Then:

```cpp
#include <matl.hpp>
using namespace matl;     // or qualify the names: matl::Matrix<double>, matl::lu, ...
```
