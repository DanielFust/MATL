# MATL

MATL (Matrix Algebra Template Library) is a header-only C++20 library for matrix algebra with
syntax close to the mathematics: `C = A*B + transpose(D) - 2.0*E;`. Expressions are evaluated
lazily with expression templates, matrix products use a cache-blocked SIMD kernel, and
temporaries are recycled so that repeated evaluations make no heap allocations.

Version 2.1.0; see [CHANGELOG.md](CHANGELOG.md) for what changed.

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

The [`examples/`](examples) folder has short, annotated programs for each feature, from the basics
to linear systems; `make run-examples` builds and runs them.

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
`MATL_VERSION` is `MAJOR * 10000 + MINOR * 100 + PATCH` (`20100` for 2.1.0).

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

## Creating and printing matrices

```cpp
Matrix<double> A(3, 4);                   // 3x4 dynamic matrix (elements uninitialized)
Matrix<double, 2, 2> F;                   // 2x2 fixed-size matrix

Matrix<double> B = {{1, 2, 3},            // from an initializer list, one inner list per row
                    {4, 5, 6}};
F = {{1, 0}, {0, 1}};                     // assignment works the same way

auto Z = zeros<double>(3, 4);             // also ones<T>, eye<T>; (n) gives n x n
auto I = eye<double, 3>();                // fixed-size forms: <T, M, N>() and <T, N>()
A = 2.0;                                  // set every element

std::cout << B << '\n';                   // rows on separate lines, columns aligned:
                                          // 1  2  3
                                          // 4  5  6
std::cout << std::setprecision(3) << B / 3.0 << '\n';   // expressions too; stream formatting applies
print(B);                                 // the same, to std::cout, with a newline
```

The rows of an initializer list must have the same length, and a fixed-size matrix needs a list
of its exact dimensions; both are always checked (`matl::dimension_error`). Assigning a list to
an allocated matrix requires the same dimensions, like any assignment; `A = {}` empties `A`.

## Elementwise operations

| Operation | Meaning |
| --- | --- |
| `A + B`, `A - B`, `-A` | Sum, difference, negation |
| `A * s`, `s * A`, `A / s`, `A + s`, `s - A`, ... | With a scalar `s`, elementwise |
| `times(A, B)`, `rdivide(A, B)`, `pow(A, B)` | Elementwise product, quotient, power (MATLAB `.*`, `./`, `.^`); either operand may be a scalar, e.g. `pow(A, 2)`, `rdivide(1.0, A)` |
| `abs`, `sqrt`, `cbrt`, `exp`, `log`, `log2`, `log10` | Elementwise math functions, as in `<cmath>` |
| `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh` | |
| `floor`, `ceil`, `round` | |
| `A * B`, `transpose(A)` | Matrix product, transpose |

All of them build expressions, so `R = sqrt(abs(A)) + 2.0 * times(A, B);` is evaluated in a single
pass over `R`. The element type of the result follows the function, e.g. `sqrt` of an integer
matrix has `double` elements.

## Reductions

| Function | Result |
| --- | --- |
| `sum(A)`, `prod(A)`, `mean(A)` | Sum, product and mean of all the elements |
| `min(A)`, `max(A)` | Smallest and largest element |
| `trace(A)` | Sum of the diagonal (`A` square) |
| `dot(A, B)` | Sum of the elementwise products; for vectors, the dot product |
| `norm_fro(A)` | Frobenius norm (for vectors, the Euclidean 2-norm) |
| `norm_1(A)`, `norm_inf(A)` | Largest absolute column sum, largest absolute row sum |

Unlike MATLAB, where `sum(A)` of a matrix gives the column sums, each of these gives a single
scalar for the whole matrix. Expressions are reduced element by element without a temporary
matrix, e.g. `max(abs(A - B))` reads `A` and `B` once. `mean` and `norm_fro` of an integer matrix
are `double`; `mean`, `min` and `max` of an empty matrix are reported as `matl::dimension_error`.

## Views

A view selects rows and columns of a matrix without copying them. Reading a view reads the
matrix, and assigning to a view writes into it:

```cpp
A(range(1, 3), all)          // rows 1 to 3                    MATLAB: A(1:3, :)
A(all, 2)                    // column 2                       MATLAB: A(:, 2)
A(4, range(2, last))         // row 4, columns 2 to the end    MATLAB: A(4, 2:end)
A(range(1, 2, last), all)    // every other row                MATLAB: A(1:2:end, :)
A(range(last, -1, 1), all)   // the rows in reverse order      MATLAB: A(end:-1:1, :)

A(all, 1) = 0.0;                         // sets column 1 of A to 0
A(range(1, 3), range(1, 3)) = B * C;     // a product written into a block of A
A(1, all) += 2.0 * A(2, all);            // a row operation
double s = sum(A(all, 2));               // a view works anywhere a matrix does
double x = A(last, last);                // last also works for single elements
```

| Index | Selects | MATLAB |
| --- | --- | --- |
| `i` | row or column `i` | `i` |
| `all` | every row or column | `:` |
| `range(first, last)` | `first`, `first+1`, ..., `last` | `first:last` |
| `range(first, step, last)` | `first`, `first+step`, ..., up to `last` (the step may be negative) | `first:step:last` |
| `last`, `last - k` | the last row or column, the `k`-th before it | `end`, `end-k` |

Ranges include their last index, as in MATLAB, and a range whose step points away from its last
index is empty: `range(k + 1, last)` selects nothing when `k` is the last index. A view can be
indexed again (`v(range(1, 2), all)` selects from the view), and a view of a row or column
vector can be indexed with one index, `A(all, j)(i)`.

**A view refers to the matrix.** Copying a view gives another view of the same elements. To copy
the elements, make a matrix:

```cpp
auto v = A(all, 1);              // v is a view of column 1 of A
v = 0.0;                         // writes column 1 of A
Matrix<double> c = A(all, 1);    // c is a new matrix holding a copy of column 1
```

Like a stored expression, a view is valid only while the matrix keeps its memory: the matrix must
not be destroyed, moved from, or reallocated (e.g. by `A = {}`). A view of a temporary matrix,
such as `ones<double>(3, 3)(all, 1)`, keeps the temporary alive. Views of const matrices and of
temporaries are read-only. Views have the type `MatrixView<T, StorageOrder>` (with `const T` for
read-only views); name them with `auto`.

`swap(v, w)` exchanges the elements of two views with the same dimensions, e.g.
`swap(A(1, all), A(3, all))` swaps rows 1 and 3. Call it unqualified or as `matl::swap`:
`std::swap(v, w)` would exchange copies of the views instead (copying a view does not copy
elements), leaving both with `w`'s elements.

**Assignments.** A view cannot be resized, so a value assigned to it must have its dimensions
(checked with `MATL_DIMENSION_CHECK`). Indices outside the matrix are reported with
`MATL_BOUNDS_CHECK`, and a step of 0 is always reported. An assignment may read the elements it
writes, e.g. `A(range(2, last), all) = A(range(1, last - 1), all)` shifts the rows down: when the
target shares elements with a value it reads at other positions, the value is evaluated into a
temporary first. Assignments that read each element only where they write it, such as
`v = 2.0 * v + B`, need no temporary.

**Performance.** Products read views directly, without copying them. Elementwise operations visit
a view's elements in the matrix's storage order, so a block of rows of a row-major matrix (or of
columns of a column-major one) is read sequentially.

## Linear systems

Dense linear systems are solved with factorizations named for their algorithms, so that the
choice, and its requirements, are explicit:

| Function | Solves | Requires |
| --- | --- | --- |
| `lu(A)`, `lu_solve(A, B)` | `A X = B` by LU factorization with partial pivoting, `P A = L U` | `A` square |
| `cholesky(A)`, `cholesky_solve(A, B)` | `A X = B` by Cholesky factorization, `A = L Lᵀ` | `A` symmetric positive definite (only its lower triangle is read) |
| `qr(A)`, `qr_solve(A, B)` | `A X ≈ B` in the least-squares sense (minimizes `‖A X − B‖`) by Householder QR, `A = Q R`; for a square `A`, `A X = B` | `A` m x n with m ≥ n and full column rank (`qr(A)` itself accepts any shape) |
| `solve_lower(L, B)`, `solve_upper(U, B)` | `L X = B`, `U X = B` by forward or back substitution | `L`, `U` triangular (only that triangle is read); `Diagonal::Unit` takes the diagonal to be 1 |

```cpp
auto F = lu(A);                    // factorize once ...
Matrix<double> x = F.solve(b);     // ... and solve for as many right-hand sides as needed
Matrix<double> X = F.solve(B);     // each column of B is one system

if (F.rcond() < 1e-12)             // estimated reciprocal condition number
    std::cerr << "A is nearly singular\n";
double d = F.determinant();

auto C = cholesky(S);              // S symmetric positive definite: half the work of LU
C.solve_inplace(X);                // overwrites X (a matrix or a view) with the solution

auto Q = qr(V);                    // V tall (more equations than unknowns), e.g. a data fit
Matrix<double> c = Q.solve(y);     // least squares: c minimizes ||V c - y||
```

A factorization object holds the factors packed in one matrix, as LAPACK does: for LU, `U` on and
above the diagonal and `L` below it (its unit diagonal is not stored), plus the row interchanges
(`pivots()`). `F.L()`, `F.U()` and `F.P()` return the factors as separate matrices, with
`P A = L U`. Matrices of integers are factorized in `double`; `float` stays `float`.

**QR and least squares.** `qr(A)` writes `A` as `Q R`, with `Q` orthogonal (a product of
Householder reflections) and `R` upper triangular. `F.Q()` and `F.R()` are the thin factors (`Q` is
m x min(m, n), and `F.Q() * F.R()` is `A`), `F.full_Q()` is the m x m `Q`, and `F.Qt_times(B)` and
`F.Q_times(B)` multiply by `Qᵀ` and `Q` without forming `Q`. For least squares, QR is the accurate
choice: solving the normal equations `AᵀA x = Aᵀb` (e.g. with Cholesky) is cheaper, but squares the
condition number of `A`. QR exists for every matrix, so `qr(A)` never fails; `solve` reports a matrix
that is exactly rank deficient (a zero on the diagonal of `R`) as `matl::singular_matrix_error`,
and `rcond()` estimates the condition number of `R`.

**In place.** `lu_inplace(A)`, `cholesky_inplace(A)` and `qr_inplace(A)` factorize a matrix in
its own memory instead of a copy: `A` then holds the factors (for Cholesky, `L` with zeros above
the diagonal), and must outlive the factorization object, which uses its memory.

**Singular and ill-conditioned matrices.** A pivot that is exactly zero means the matrix is
singular, and is reported as `matl::singular_matrix_error` (`matl::not_positive_definite_error`
for a Cholesky pivot that is not positive). A matrix that is merely close to singular is
factorized, and `rcond()` tells how close: it estimates `1 / (‖A‖₁ ‖A⁻¹‖₁)` from a few triangular
solves (Hager and Higham's method, as in LAPACK). It is near 1 for a well-conditioned matrix, and a
solution loses about `log10(1 / rcond)` of the 16 significant digits of a `double`. The
dimensions are always checked (`matl::dimension_error`), since the check costs nothing next to the
factorization.

**Performance.** The factorizations are recursive (Toledo's LU, a recursive Cholesky, and
LAPACK's recursive QR panel, `xGEQRT3`, inside a blocked QR): they split the matrix in halves down
to 16 columns, so that almost all the arithmetic is matrix products of every size, which run on
the product kernel through views and subtract from the matrix without temporaries. At n = 1024
(Intel i9-9880H, GCC 15 and Clang 21, `-O2 -march=native`), LU runs at about 15–17 GFLOP/s,
Cholesky at about 11–12 GFLOP/s and QR at about 13–15 GFLOP/s (9–11 GFLOP/s for a 4000 x 400
matrix), against 21–25 GFLOP/s for the product.

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
| `matl` | Everything you use directly: matrix types (see Storage order), operators and functions, `allocate`/`deallocate`, `MemoryManager`, and the concepts below |
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

## Storage order

Matrices store their elements row by row (row-major, the default) or column by column
(column-major, as in Fortran, LAPACK and MATLAB):

| Type | Dimensions | Storage order |
| --- | --- | --- |
| `Matrix<T>` = `RowMajorMatrix<T>` | run time | row-major |
| `ColMajorMatrix<T>` | run time | column-major |
| `Matrix<T, M, N>` = `FixedSizeRowMajorMatrix<T, M, N>` | compile time | row-major |
| `FixedSizeColMajorMatrix<T, M, N>` | compile time | column-major |

These are aliases of two class templates, `DynamicMatrix<T, StorageOrder>` and
`FixedSizeMatrix<T, M, N, StorageOrder>` (the names compiler messages show). Both storage orders
are used the same way: `A(i,j)` is row `i`, column `j`, initializer lists give the rows, and
matrices of either order can be mixed freely in expressions and assignments. Only the single
index `A(i)` and external memory follow the storage order.

Column-major matrices are mainly useful for exchanging data with Fortran or LAPACK without
copying, e.g. wrapping an array a Fortran routine filled:

```cpp
ColMajorMatrix<double> A(m, n, fortran_array);   // A(i,j) is fortran_array(i,j)
```

## Matrix memory

A matrix owns its data, which it allocates on the heap. Copies are deep, and moves
transfer the data without copying it.

A matrix can also use memory that you own, e.g. to place many matrices in one contiguous block
or to wrap a buffer received through MPI. The matrix reads and writes the buffer in place (in
its storage order) and never frees it, so the buffer must outlive the matrix:

```cpp
std::vector<double> block(2 * 16);
RowMajorMatrix<double> A(4, 4, block.data());        // A uses block[0..15]
RowMajorMatrix<double> B(4, 4, block.data() + 16);   // B uses block[16..31]
A = B * 2.0;                                         // result written into the block
```

Assigning to such a matrix (including move assignment) copies the values into the buffer. A copy
of it owns its own data. `allocate(A, rows, cols)` and `deallocate(A)` allocate a matrix on the
heap and release its memory; for a matrix over external memory, `deallocate` only stops using
the buffer. `A.data()` points to a matrix's elements in storage order, e.g. to pass the matrix to
a C or Fortran library.

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
| `matl::singular_matrix_error` | `std::runtime_error` | A factorization or solve meets an exactly zero pivot (a singular, or for QR rank-deficient, matrix) |
| `matl::not_positive_definite_error` | `std::runtime_error` | A Cholesky factorization meets a pivot that is not positive |
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
a program would normally catch them only to log context or shut down cleanly. The errors of the
linear solvers are always checked, and a singular or indefinite matrix may well be a property of
the data rather than a bug, so these are worth catching where they can occur. A caught error
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
storage order (along the rows for a row-major matrix, down the columns for a column-major one).
`last` stands for the last row or column, as MATLAB's `end` does: `A(last, 1)`, `A(2, last - 1)`.
Rows, columns and blocks are selected with views (see [Views](#views)).

## Tests

See [tests/README.md](tests/README.md). In short, from the repository root, with Make:

```sh
make check            # build and run the test suite
make run-examples     # build and run the example programs
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
