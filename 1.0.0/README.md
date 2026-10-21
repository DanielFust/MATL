# MATL 1.0.0

MATL (Matrix Algebra Template Library) is a header-only C++ library for matrix algebra with
syntax close to the mathematics: `C = A*B + transpose(D) - 2.0*E;`. Expressions are evaluated
lazily with expression templates.

```cpp
#include "matl.hpp"
```

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

## Debug checks

`matl.hpp` enables dimension checks, bounds checks and `MATL_DEBUG` trace output by default.
Define `MATL_NO_DEBUG` before including it (or compile with `-DMATL_NO_DEBUG`) to disable them,
e.g. for benchmarks and production builds.

## Indexing

Matrices use **1-based indexing**: `A(1,1)` is the first element, and `A(i)` indexes elements in
storage (row-major) order.

## Tests

See [tests/README.md](tests/README.md). In short, from this directory:

```sh
make check            # build and run the test suite
make check-sanitize   # with AddressSanitizer and UndefinedBehaviorSanitizer
make check-nodebug    # with MATL_NO_DEBUG
make bench            # benchmarks
```
