# MATL

MATL (Matrix Algebra Template Library) is a header-only C++20 library for matrix algebra with
syntax close to the mathematics: `C = A*B + transpose(D) - 2.0*E;`. Expressions are evaluated
lazily with expression templates, matrix products use a cache-blocked SIMD kernel, and
temporaries are recycled so that repeated evaluations make no heap allocations.

Version 2.0.0.

```cpp
#include <matl.hpp>
using namespace matl;

int main()
{
    Matrix<double> A(3, 3), B(3, 3), C(3, 3);     // dynamic size
    Matrix<double, 3, 3> F;                       // fixed size (no heap memory)
    A = 1.0;
    B = eye<double>(3);
    C = A*B + transpose(A) - 2.0*B;               // A*B into a recycled temporary, then one pass
    F = C * C;
    print(F);
}
```

## Using MATL in your project

MATL is header-only: add its `include/` directory to your include path and compile with
`-std=c++20`. With CMake, use the `matl::matl` target in either of these ways:

```cmake
# After installing MATL (cmake -S . -B build && cmake --install build --prefix <prefix>)
find_package(matl 2 REQUIRED)
target_link_libraries(my_program PRIVATE matl::matl)

# From a copy of the source tree (a git submodule, or FetchContent)
add_subdirectory(MATL)
target_link_libraries(my_program PRIVATE matl::matl)
```

The `matl::matl` target adds the include directory and requires C++20. In a version check,
`MATL_VERSION` is `MAJOR * 10000 + MINOR * 100 + PATCH` (`20000` for 2.0.0).

## Requirements

MATL uses **C++20 language features** (concepts and `requires` clauses) and is compiled with
`-std=c++20`. To stay portable across HPC systems it avoids C++20 features with uneven compiler
support: modules, coroutines, `<format>`, heavy use of `<ranges>`, and omitting `typename`
before dependent type names (P0634, which Clang only supports from version 16).

Minimum supported compilers, tested in CI:

| Compiler | Minimum version |
| --- | --- |
| GCC | 10 |
| Clang | 13 |
| Apple Clang | current Xcode |

Most vendor HPC compilers are based on LLVM (Intel oneAPI `icpx`, HPE/Cray CCE, AMD ROCm,
IBM Open XL); recent versions support the features MATL uses. The minimum versions above are the
ones verified by CI.

## Expressions, values and lifetimes

An arithmetic expression on matrices produces an *expression object* (a formula), not a matrix.
It is evaluated when it is assigned to a matrix:

```cpp
auto           e = A + B + C;   // formula: nothing is computed yet
Matrix<double> m = A + B + C;   // value: evaluated once, here
D = e;                          // evaluates the formula into D
```

Expression objects are safe to store, copy and return from functions:

- **Named matrices are held by reference.** A stored expression is a *live* formula: it reads the
  matrices' current values each time it is evaluated. A stored expression must not outlive the
  named matrices it refers to (the same rule as for any reference).
- **Sub-expressions are held by value**, so the temporaries created while building an expression
  (such as `A + B` in `A + B + C`) can safely go out of scope.
- **Temporary matrices are moved into the expression**, which then owns them, e.g.
  `auto e = A + eye<double>(3);`.

To compute a result once and reuse it, assign it to a matrix.

### When work happens

Building an expression computes nothing and allocates nothing. All work happens when the
expression is assigned to a matrix:

```cpp
D = A*B + transpose(C*A) - E;
// 1. A*B and C*A are computed into temporary matrices
// 2. the elementwise part is evaluated in a single pass: D(i,j) = T1(i,j) + T2(j,i) - E(i,j)
// 3. the temporaries are released (in reverse order of allocation)
```

- A product assigned directly (`C = A*B`) is computed straight into the target, with no
  temporary.
- Operands of a product that are expressions (e.g. `A + B` in `(A + B)*C`) are first evaluated
  into temporaries, because each of their elements is read many times.
- Temporaries for fixed-size matrices are plain arrays and need no allocation. Temporaries for
  dynamic matrices come from a recycled arena (see below).
- Because products are computed before any element of the target is written, the target may
  appear anywhere inside them (`C = A*C + C`, `C = transpose(C*A)`). MATL also detects when the
  target is read in a way that would be overwritten too early (`C = A*C`, `C = C + transpose(C)`)
  and evaluates into a temporary first.

A stored expression is re-evaluated each time it is assigned, so assigning an expression that
contains a product to several matrices computes the product each time.

### Memory for temporaries

Temporaries are allocated from a per-thread `TemporaryArena`, a stack of preallocated memory
that is reused by every assignment. Once the arena has grown to fit the temporaries of an
expression, evaluating it makes **no heap allocations**. The arena grows automatically, and you
can reserve memory up front so that even the first evaluation allocates nothing:

```cpp
MemoryManager::temporary_arena().reserve(64 * 1024 * 1024);    // bytes

// Statistics
MemoryManager::temporary_arena().bytes_in_use();      // 0 between assignments
MemoryManager::temporary_arena().high_water_mark();   // peak bytes used by temporaries
MemoryManager::temporary_arena().release_memory();    // return the arena's memory to the system
```

Set `MemoryManager::CACHE_TEMPORARY_MATRICES = false` to allocate temporaries on the heap
instead. Each thread has its own arena, which is freed when the thread exits.

`make bench` reports the time and heap allocations per assignment for several expressions,
with and without the arena.

## Types and namespaces

| Namespace | Contents |
| --- | --- |
| `matl` | Everything you use directly: matrix types (`Matrix`, `RowMajorMatrix`, `FixedSizeRowMajorMatrix`), operators, `transpose`, `eye`, `allocate`/`deallocate`, `MemoryManager`, and the concepts below |
| `matl::expr` | Expression types, i.e. what `auto e = A + B;` holds: `Elementwise<Op, L, R>` (every `+`, `-`, and `*` or `/` by a scalar), `Product<L, R>`, `Transpose<M>`, `Owned<M>` |
| `matl::detail` | Implementation details (evaluators, traits); not part of the interface |

You rarely need to name an expression type, but they appear in compiler messages and debuggers,
e.g. `A + B*2.0` has type
`matl::expr::Elementwise<matl::expr::Add, RowMajorMatrix<double>, matl::expr::Elementwise<matl::expr::Multiply, RowMajorMatrix<double>, double>>`.

Write generic functions with the concepts:

| Concept | Satisfied by |
| --- | --- |
| `MatrixType` | Matrices (any MATL matrix type) |
| `ExpressionType` | Expressions |
| `MatrixAlgebraObject` | Matrices and expressions |
| `Scalar` | Arithmetic types (`int`, `double`, ...) |
| `FixedSize` | Matrices and expressions with compile-time dimensions |

```cpp
template<matl::MatrixType MAT>
double trace(MAT const& A)
{
    double sum = 0;
    for (matl::index_t i = 1; i <= A.rows(); ++i) sum += A(i,i);
    return sum;
}
```

The operators only accept MATL matrices, expressions and arithmetic scalars, so a misuse such as
`A + "text"` is reported as an unsatisfied constraint rather than an error inside the library.

## Matrix memory

A `RowMajorMatrix` owns its data, which it allocates on the heap. Copies are deep, and moves
transfer the data without copying it.

A matrix can also use memory that you own, e.g. to place many matrices in one contiguous block
or to wrap a buffer received through MPI. The matrix reads and writes the buffer in place (in
row-major order) and never frees it, so the buffer must outlive the matrix:

```cpp
std::vector<double> block(2 * 16);
RowMajorMatrix<double> A(4, 4, block.data());        // A uses block[0..15]
RowMajorMatrix<double> B(4, 4, block.data() + 16);   // B uses block[16..31]
A = B * 2.0;                                         // result written into the block
```

Assigning to such a matrix (including move assignment) copies the values into the buffer. A copy
of it owns its own data. `allocate(A, rows, cols)` and `deallocate(A)` allocate a matrix on the
heap and release its memory; for a matrix over external memory, `deallocate` only stops using
the buffer.

## Dimension checks

Dimensions known at compile time are checked at compile time. Adding a 2x3 and a 3x2 fixed-size
matrix, multiplying matrices whose inner dimensions differ, or assigning a result to a fixed-size
matrix of another size does not compile:

```cpp
Matrix<double,2,3> A;  Matrix<double,3,2> B;
auto e = A + B;    // error: ... because 'SizesMatch<...<double,2,3>, ...<double,3,2>>' evaluated to false
auto p = A * B;    // OK: 2x3 * 3x2
```

The concepts `SizesMatch<A, B>` and `InnerSizesMatch<A, B>` express these rules. When a dimension
is only known at run time (dynamic matrices, or a mix of fixed-size and dynamic operands), the
check happens when the expression is assigned, if dimension checks are enabled (see below):
operands of different sizes, or a product with different inner dimensions, are reported as a
`matl::dimension_error` (see Errors). Because the check runs at assignment, it also catches a
stored expression whose matrices were reallocated after it was built.

## Debug checks

Like `assert`, MATL's run-time checks are enabled unless `NDEBUG` is defined, so debug builds
check and release builds (`-DNDEBUG`, the default of CMake's Release configuration) do not:

| Macro | Checks |
| --- | --- |
| `MATL_BOUNDS_CHECK` | Indices of element access |
| `MATL_DIMENSION_CHECK` | Dimensions of assignments and of the operations in an expression |

Define `MATL_NO_DEBUG` to disable the checks in any build, or define one of the macros above
before including `matl.hpp` to keep that check in a release build. Compile-time dimension checks
are always on.

Define `MATL_DEBUG` to print a trace of internal calls (assignments, allocations) to stdout. It is
off by default.

## Errors

MATL reports the errors it detects by throwing an exception:

| Exception | Derived from | Reported when |
| --- | --- | --- |
| `matl::dimension_error` | `std::invalid_argument` | Operands, or an assignment target, have incompatible dimensions |
| `matl::index_error` | `std::out_of_range` | An element index is out of range (indices are 1-based) |
| `std::logic_error` | | The `TemporaryArena` is misused (e.g. releasing a block twice) |

```cpp
try { R = A + B; }
catch (matl::dimension_error const& e) { std::cerr << e.what() << '\n'; }
// incompatible matrix dimensions in elementwise operation: 2x2 and 3x3
```

Define `MATL_ABORT_ON_ERROR` to print the message and call `std::abort()` instead, which stops at
the faulty statement under a debugger. MATL also aborts when it is compiled without exceptions
(e.g. `-fno-exceptions`). `MATL_THROWS_ON_ERROR` is `1` when errors throw and `0` when they abort.

Most of these errors are found by the debug checks above, so they indicate a programming error;
a program would normally catch them only to log context or shut down cleanly. A caught error
leaves MATL in a usable state: the temporaries of the failed evaluation are released (the target
of a failed assignment may be partly written).

## Performance

Matrix products use a cache-blocked kernel in the structure of GotoBLAS/BLIS: blocks of the
operands are packed into contiguous buffers (taken from the temporary arena, so they make no heap
allocations) and multiplied in small tiles held in SIMD registers. Small products (below about
10x10x10) use a simple loop, which is faster there. Each element of a product is still the sum of
its products in increasing order, so without fused multiply-add (FMA) instructions blocking
gives exactly the result of a simple loop. When FMA is available (e.g. `-march=native` on recent
x86, or ARM), compilers fuse multiplications and additions by default; results may then differ
from a simple loop in the last bits. Compile with `-ffp-contract=off` to disable fusion.

Elementwise expressions are evaluated in a single pass, at the speed of a hand-written loop.
Transposes, and expressions mixing storage orders, are evaluated in small tiles so that both
the rows and the columns being accessed stay in cache.

Compile with optimization and for your processor's vector instructions for the best performance,
e.g. `-O3 -march=native` (GCC, Clang) — products use 32-byte AVX vectors when they are enabled,
and 16-byte vectors (SSE2, NEON) otherwise. On the development machine (Intel i9-9880H), a
1024x1024 double-precision product runs at about 15 GFLOP/s with `-O2` and about 27 GFLOP/s with
`-O3 -march=native`, compared with 1.4 GFLOP/s for the simple loop. `make bench` builds and runs
the benchmarks in `benchmarks/`.

The kernel's block sizes are in `detail::product_blocking` (`_matl_product_kernel_.hpp`).

## Indexing

Matrices use **1-based indexing**: `A(1,1)` is the first element, and `A(i)` indexes elements in
storage (row-major) order.

## Tests

See [tests/README.md](tests/README.md). In short, from the repository root, with Make:

```sh
make check            # build and run the test suite
make check-sanitize   # with AddressSanitizer and UndefinedBehaviorSanitizer
make check-nodebug    # with MATL_NO_DEBUG
make bench            # benchmarks
```

or with CMake:

```sh
cmake -S . -B build/cmake -DMATL_BUILD_BENCHMARKS=ON
cmake --build build/cmake
cd build/cmake && ctest
```

## License and citation

MATL is released under the [BSD 3-Clause License](LICENSE). If you use MATL in your research,
please cite it: [CITATION.cff](CITATION.cff) has the citation metadata, and GitHub's
"Cite this repository" button formats it (e.g. as BibTeX).
