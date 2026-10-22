# MATL test suite

## Running the tests

From the repository root:

```sh
make check                      # build and run the test suite
make check-sanitize             # same, with AddressSanitizer + UndefinedBehaviorSanitizer
make check-nodebug              # same, with MATL_NO_DEBUG (run-time debug checks compiled out)
make check-abort                # same, with MATL_ABORT_ON_ERROR (errors abort instead of throwing)
make check-noexceptions         # check that MATL compiles with -fno-exceptions and errors abort
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
| `test_expressions.cpp` | Elementwise, scalar, product and transpose expressions; aliasing; fixed-size and mixed operands |
| `test_evaluation.cpp` | Deferred evaluation: nested products read current operand values, no allocation until assignment, temporaries released afterwards, aliasing with deferred products |
| `test_construction.cpp` | Initializer lists, `zeros` and `ones`, `operator<<` and `print` |
| `test_functions.cpp` | Negation, elementwise math functions, `times`, `rdivide` and `pow`; dimensions, lifetimes and aliasing |
| `test_types.cpp` | Concepts, the expression types the operators build, traits (fixed sizes, storage order, element-type promotion), assignment forms, compile-time dimension checks |
| `test_expression_lifetime.cpp` | Stored expressions (`auto e = ...`): ownership of sub-expressions and temporaries, live references to named matrices, initializing matrices from expressions |
| `test_products.cpp` | The cache-blocked product kernel at sizes around its block boundaries, element types, expression operands, aliasing and fixed-size operands |
| `test_temporary_arena.cpp` | TemporaryArena (stack order, growth, out-of-order release) and the exact temporaries each kind of assignment uses |
| `test_heap_allocations.cpp` | Counts heap allocations: none during evaluation once the arena has grown. Replaces the global `operator new` for the test program (counting only inside an explicit window) |
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
