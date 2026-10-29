/**************************************************************************************************
 * File: _matl_parallel_.hpp
 *
 * Description: Thread parallelism with OpenMP.
 *
 * MATL runs its large loops on several threads when it is compiled with OpenMP (-fopenmp, which
 * defines _OPENMP), unless MATL_NO_OPENMP is defined. Without OpenMP, the same code runs serially.
 * The number of threads is OpenMP's (OMP_NUM_THREADS, omp_set_num_threads). MATL_OPENMP is 1 if
 * MATL uses OpenMP and 0 otherwise.
 *
 * What runs in parallel: the elementwise evaluation of expressions (e.g. the vector updates of
 * the iterative solvers), reductions (sum, dot, norm_fro, ...), sparse matrix times dense
 * products S * X, and the Jacobi preconditioner. A loop runs serially when its work is below
 * PARALLEL_THRESHOLD (threads cost more than they save on small problems) and inside a parallel
 * region of the caller (e.g. one matrix per thread), so MATL never nests parallelism.
 *
 * Reproducibility: reductions add their elements in blocks of REDUCTION_BLOCK elements, and add
 * the blocks' results in order, whatever the number of threads (and without OpenMP). A result is
 * therefore the same, bit for bit, on 1 thread or 16, and an iterative solver takes the same
 * iterations. Elementwise results never depend on the number of threads.
 *
 * Errors: an exception thrown inside a parallel loop (e.g. by a bounds check) is caught there and
 * rethrown after the loop (OpenMP would otherwise terminate the program).
 **************************************************************************************************/
#ifndef MATL_PARALLEL_HPP
#define MATL_PARALLEL_HPP

#if defined(_OPENMP) && !defined(MATL_NO_OPENMP)
#include <omp.h>
#define MATL_OPENMP 1
#else
#define MATL_OPENMP 0
#endif

#if defined(__cpp_exceptions) || defined(__EXCEPTIONS) || defined(_CPPUNWIND)
#define MATL_PARALLEL_EXCEPTIONS 1
#else
#define MATL_PARALLEL_EXCEPTIONS 0
#endif

namespace matl::detail
{
    // Work (elements, or stored entries times columns) below which a loop runs serially. Starting
    // the threads costs a few microseconds: on an 8-core laptop, 8 threads were slower than one
    // for vector updates of 40,000 elements and faster from about 60,000.
    inline constexpr index_t PARALLEL_THRESHOLD = index_t(1) << 16;

    // Elements per block of a reduction (see the description above)
    inline constexpr index_t REDUCTION_BLOCK = 4096;

    // 'true' if a loop with the given work should use several threads
    inline bool use_threads([[maybe_unused]] index_t work)
    {
        #if MATL_OPENMP
        return work >= PARALLEL_THRESHOLD && !omp_in_parallel() && omp_get_max_threads() > 1;
        #else
        return false;
        #endif
    }

   /***********************************************************************************************
    * void parallel_ranges(n, work, f, split)
    *
    * Calls f(first, last) for consecutive ranges [first, last) that cover [0, n), one range per
    * thread, or f(0, n) when the loop runs serially (see use_threads). f processes its range with
    * an ordinary loop, so the compiler can vectorize it.
    *
    * Parameters:
    *     n        Number of items (e.g. rows, or elements)
    *     work     Estimated work of the whole loop (decides whether to use threads)
    *     f        f(first, last) for a range of items (0-based, last excluded)
    *     split    split(t, threads): the first item of thread t's range (default: equal ranges;
    *              e.g. ranges with equal numbers of stored entries for a sparse matrix)
    ***********************************************************************************************/
    struct equal_split
    {
        index_t n;
        index_t operator() (index_t t, index_t threads) const {return n / threads * t + std::min(t, n % threads);}
    };

    template<typename F, typename SPLIT>
    void parallel_ranges(index_t n, [[maybe_unused]] index_t work, F&& f, [[maybe_unused]] SPLIT const& split)
    {
        #if MATL_OPENMP
        if (n > 1 && use_threads(work))
        {
            #if MATL_PARALLEL_EXCEPTIONS
            std::exception_ptr error = nullptr;
            #endif
            #pragma omp parallel
            {
                const index_t t = static_cast<index_t>(omp_get_thread_num());
                const index_t threads = static_cast<index_t>(omp_get_num_threads());
                const index_t first = (t == 0) ? 0 : split(t, threads);
                const index_t last  = (t + 1 == threads) ? n : split(t + 1, threads);
                #if MATL_PARALLEL_EXCEPTIONS
                try
                {
                    if (first < last) f(first, last);
                }
                catch (...)
                {
                    #pragma omp critical(matl_parallel_error)
                    if (!error) error = std::current_exception();
                }
                #else
                if (first < last) f(first, last);
                #endif
            }
            #if MATL_PARALLEL_EXCEPTIONS
            if (error) std::rethrow_exception(error);
            #endif
            return;
        }
        #endif
        if (n > 0) f(index_t(0), n);
    }

    template<typename F>
    void parallel_ranges(index_t n, index_t work, F&& f)
    {parallel_ranges(n, work, std::forward<F>(f), equal_split{n});}

}
#endif
