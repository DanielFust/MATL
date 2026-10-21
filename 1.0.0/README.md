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
- Temporaries for dynamic matrices are allocated according to the `MemoryManager` settings;
  temporaries for fixed-size matrices need no allocation.
- Because products are computed before any element of the target is written, the target may
  appear anywhere inside them (`C = A*C + C`, `C = transpose(C*A)`). MATL also detects when the
  target is read in a way that would be overwritten too early (`C = A*C`, `C = C + transpose(C)`)
  and evaluates into a temporary first.

A stored expression is re-evaluated each time it is assigned, so assigning an expression that
contains a product to several matrices computes the product each time.

Matrices use **1-based indexing**: `A(1,1)` is the first element, and `A(i)` indexes elements in
storage (row-major) order.

## Tests

See [tests/README.md](tests/README.md). In short, from this directory:

```sh
make check            # build and run the test suite
make check-sanitize   # with AddressSanitizer and UndefinedBehaviorSanitizer
```
