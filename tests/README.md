# MATL test suite

## Running the tests

From the repository root:

```sh
make check                      # build and run the test suite
make check-sanitize             # same, with AddressSanitizer + UndefinedBehaviorSanitizer
make check-nodebug              # same, with MATL_NO_DEBUG (run-time debug checks compiled out)
make check-abort                # same, with MATL_ABORT_ON_ERROR (errors abort instead of throwing)
make check-noexceptions         # check that MATL compiles with -fno-exceptions and errors abort
make check-openmp               # with OpenMP, on 4 threads (GCC, or Clang with an OpenMP runtime)
make check FILTER=Transpose     # run only test cases whose name contains "Transpose"
make check CXX=clang++          # choose a compiler (default: g++)
./build/matl_tests --list       # list all test case names
```

With CMake, the suite is built when MATL is the top-level project (or with
`-DMATL_BUILD_TESTS=ON`), and `ctest` runs it three times: with the run-time checks enabled (even
in a Release build), with `MATL_ABORT_ON_ERROR`, and with `MATL_NO_DEBUG`; it also runs the
`-fno-exceptions` program. `tests/install/` is a small project that checks an
installed MATL can be found with `find_package(matl)`; CI builds it after `cmake --install`.

Test results are printed to stderr. The test builds enable MATL's run-time checks (they define
neither `NDEBUG` nor `MATL_NO_DEBUG`). The sanitizer build also enables the `MATL_DEBUG` trace
output, which goes to stdout; the Makefile writes it to `build/matl_tests_sanitize.stdout.log`.

The tests are compiled with `-Wall -Wextra -Werror`, so new compiler warnings in the library fail
the build.

## Layout

| File | Covers |
| --- | --- |
| `matl_test.hpp` | Test framework (macros below) and the naive reference implementation `matl_test::Ref` |
| `test_main.cpp` | Test runner |
| `test_matrix.cpp` | Construction, indexing, copy/move, scalar and compound assignment, `eye` |
| `test_parallel.cpp` | Thread parallelism: ranges cover every item once, exceptions inside parallel loops are rethrown, results are bit for bit the same on 1 and 4 threads (reductions, elementwise evaluation, dense and sparse products, LU, Cholesky and QR factors, CG), MATL runs serially inside a caller's parallel region and in a `SerialScope` (also from `std::thread`), and a finite-element style element loop with per-thread builders assembles the same matrix as a serial loop |
| `test_expressions.cpp` | Elementwise, scalar, product and transpose expressions; aliasing; fixed-size and mixed operands |
| `test_evaluation.cpp` | Deferred evaluation: nested products read current operand values, no allocation until assignment, temporaries released afterwards, aliasing with deferred products |
| `test_construction.cpp` | Initializer lists, `zeros` and `ones`, `operator<<` and `print` |
| `test_functions.cpp` | Negation, elementwise math functions, `times`, `rdivide` and `pow`; dimensions, lifetimes and aliasing |
| `test_reductions.cpp` | `sum`, `prod`, `mean`, `min`, `max`, `trace`, `dot` and the norms: expressions without temporaries, mixed storage orders, result types, empty matrices, dimension checks |
| `test_solvers.cpp` | Triangular solves, LU and Cholesky factorizations: factors, residuals of solutions (small and blocked, both storage orders), views, expressions, integer and `float` matrices, in-place variants, determinants, condition estimates against exact values, singular and non-positive-definite matrices, dimension errors |
| `test_sparse.cpp` | Sparse matrices: assembly (duplicates, explicit zeros, 32-bit indices), checked construction from CSR arrays, element access and changing values, conversions, products with dense vectors and matrices (transposed, nested, views, column-major, targets read by the product), stored expressions, scaling, output, dimension errors |
| `test_structured.cpp` | Sparse arithmetic (sums, differences, products, transposes, the diagonal, result patterns), tridiagonal solves, banded matrices (band access, conversions, products), banded LU and Cholesky against dense factorizations (residuals, determinants, condition estimates), errors |
| `test_column_major.cpp` | Column-major matrices: layout, Fortran-order external memory, every kind of operation across storage orders, conversions, aliasing |
| `test_types.cpp` | Concepts, the expression types the operators build, traits (fixed sizes, storage order, element-type promotion), assignment forms, compile-time dimension checks |
| `test_views.cpp` | Views with `range`, `all` and `last`: selecting and writing rows, columns, blocks and steps, views in expressions, products and reductions, views of views, column-major and fixed-size matrices, assignments that overlap what they read, views of temporaries, errors |
| `test_expression_lifetime.cpp` | Stored expressions (`auto e = ...`): ownership of sub-expressions and temporaries, live references to named matrices, initializing matrices from expressions |
| `test_products.cpp` | The cache-blocked product kernel at sizes around its block boundaries, element types, expression operands, aliasing and fixed-size operands |
| `test_qr.cpp` | QR factorization and least squares: `Q R = A`, orthonormal `Q` and triangular `R` for tall, square and wide matrices (small and blocked, both storage orders), least-squares solutions (exact fits, normal-equation residuals, square systems against LU), products with `Q`, views, expressions, integer and `float` matrices, in place, condition estimates, rank deficiency and errors |
| `test_temporary_arena.cpp` | TemporaryArena (stack order, growth, out-of-order release) and the exact temporaries each kind of assignment uses |
| `test_heap_allocations.cpp` | Counts heap allocations: none during evaluation once the arena has grown. Replaces the global `operator new` for the test program (counting only inside an explicit window) |
| `test_iterative.cpp` | Iterative solvers and preconditioners: CG, BiCGSTAB and GMRES on the 2D Laplacian and a nonsymmetric convection-diffusion matrix (against dense LU), Jacobi and ILU(0) (fewer iterations; exact for a tridiagonal matrix), restart lengths, dense matrices and operator functions, vector types, `float`, non-convergence, errors |
| `test_memory_manager.cpp` | Matrix memory ownership (heap and external buffers), `allocate`/`deallocate`, MemoryManager settings, per-thread arenas |
| `test_errors.cpp` | Error reporting: run-time bounds and dimension checks (including operands of different sizes), exception types and messages, recovery after a caught error |
| `no_exceptions/main.cpp` | A separate program built with `-fno-exceptions`: MATL compiles, and errors abort |
| `test_randomized.cpp` | Random sizes and contents (fixed seed) compared against the reference implementation |

Each test file is a separate translation unit, so the build also verifies that `matl.hpp` can be
included from more than one `.cpp` file.

## Writing tests

```cpp
#include "matl_test.hpp"
using namespace matl;
using matl_test::Ref;

TEST_CASE("Expressions: my new case")
{
    RowMajorMatrix<double> A(2, 3), B(3, 2), C(2, 2);
    matl_test::fill_sequence(A);              // 1, 2, 3, ... in row-major order
    matl_test::fill_random(B);                // small random integers (exact arithmetic)

    C = A * B;
    CHECK_MATRIX_EQ(C, Ref::from(A) * Ref::from(B));
}
```

| Macro | Meaning |
| --- | --- |
| `CHECK(cond)` | Record a failure and continue |
| `REQUIRE(cond)` | Record a failure and stop the test case |
| `CHECK_EQ(a, b)` | Equality, printing both values on failure |
| `CHECK_THROWS(expr)` / `CHECK_THROWS_AS(expr, type)` | Expect an exception |
| `CHECK_ABORTS(statements)` | Expect the process to terminate (runs in a forked child) |
| `CHECK_MATL_ERROR(type, statements)` | Expect a MATL error: the exception `type`, or termination with `MATL_ABORT_ON_ERROR` |
| `CHECK_MATRIX_EQ(mat, ref)` | Compare a MATL matrix element-by-element with a `Ref` |

Note that MATL only sizes a matrix automatically when it is empty: assigning a result of a
different size to an already-allocated matrix is a dimension error, so use a separate matrix for
each result size.

## Continuous integration

`.github/workflows/ci.yml` runs on pull requests, on pushes to `main`, and by hand: six compilers
(GCC 10, 13 and 14, Clang 13 and 18, Apple Clang) with every mode above, the examples, the
benchmarks, and the CMake build, install and `find_package`.

When GitHub Actions minutes are limited (macOS runners count ten times; a full run costs about
80-100 minutes), switch to the light CI, `.github/workflows/ci-light.yml`: two Linux jobs on the
minimum compilers (GCC 10 and Clang 13), on pull requests only, with every mode (the sanitizers
with Clang, OpenMP and CMake with GCC), about 15-20 minutes per run. The switch is a repository
variable, so no workflow file changes:

```sh
gh variable set CI_LIGHT --body true      # light CI from now on (the full CI's jobs are skipped)
gh variable delete CI_LIGHT               # back to the full CI
gh workflow run ci.yml --ref <branch>     # the full CI by hand, in either mode
```

Pushing a branch without opening a pull request runs no CI.

