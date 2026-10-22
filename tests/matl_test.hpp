/**************************************************************************************************
 * File: matl_test.hpp
 *
 * Description: Minimal, dependency-free test framework for the MATL test suite.
 *
 *   TEST_CASE("name") { ... }       Registers a test case (any number per translation unit)
 *   CHECK(cond)                     Records a failure if cond is false, continues the test
 *   REQUIRE(cond)                   Records a failure and aborts the current test case
 *   CHECK_EQ(a, b)                  Equality check that prints both values on failure
 *   CHECK_THROWS(expr)              Checks that expr throws any exception
 *   CHECK_THROWS_AS(expr, type)     Checks that expr throws the given exception type
 *   CHECK_ABORTS(statements)        Checks that the statements terminate the process (run in a
 *                                   forked child process, so the test runner survives)
 *   CHECK_MATRIX_EQ(mat, ref)       Compares a MATL matrix element-by-element against a
 *                                   matl_test::Ref reference matrix
 *
 * The matl_test::Ref type is a deliberately naive reference implementation of matrix
 * arithmetic used to check MATL results.
 **************************************************************************************************/
#ifndef MATL_TEST_HPP
#define MATL_TEST_HPP

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

#include "matl.hpp"

namespace matl_test
{
    /*************************************** Registry ***************************************/

    struct TestCase
    {
        const char*           name;
        const char*           file;
        int                   line;
        std::function<void()> body;
    };

    inline std::vector<TestCase>& registry()
    {
        static std::vector<TestCase> tests;
        return tests;
    }

    struct Registrar
    {
        Registrar(const char* name, const char* file, int line, void (*body)())
        {registry().push_back({name, file, line, body});}
    };

    // Thrown by REQUIRE to abort the current test case
    struct RequireFailure {};

    // Failure count for the currently running test case
    inline int& current_failures()
    {
        static int failures = 0;
        return failures;
    }

    inline void report_failure(const char* file, int line, const std::string& message)
    {
        std::fprintf(stderr, "    %s:%d: FAILED: %s\n", file, line, message.c_str());
        ++current_failures();
    }

    template<typename T>
    std::string to_string(T const& value)
    {
        std::ostringstream out;
        out << value;
        return out.str();
    }

    /********************************** Death test support **********************************/

    // Runs 'body' in a forked child process. Returns true if the child terminated abnormally
    // (signal such as SIGABRT, or a non-zero exit status) and false if 'body' returned normally.
    inline bool terminates(const std::function<void()>& body)
    {
        std::fflush(stdout);
        std::fflush(stderr);
        pid_t pid = fork();
        if (pid < 0) throw std::runtime_error("fork() failed");
        if (pid == 0)
        {
            // Child: silence output, run the body, exit cleanly if it returns
            if (std::freopen("/dev/null", "w", stdout) == nullptr) _exit(2);
            if (std::freopen("/dev/null", "w", stderr) == nullptr) _exit(2);
            try { body(); } catch (...) { _exit(3); }
            _exit(0);
        }
        int status = 0;
        waitpid(pid, &status, 0);
        return !(WIFEXITED(status) && WEXITSTATUS(status) == 0);
    }

    /********************************* Reference matrices ***********************************/

    // Naive row-major reference matrix with 1-based indexing, used to verify MATL results
    struct Ref
    {
        matl::index_t       rows = 0;
        matl::index_t       cols = 0;
        std::vector<double> data;

        Ref() = default;
        Ref(matl::index_t m, matl::index_t n, double value = 0.0) : rows(m), cols(n), data(m*n, value) {}

        double& operator()(matl::index_t i, matl::index_t j)       {return data[(i-1)*cols + (j-1)];}
        double  operator()(matl::index_t i, matl::index_t j) const {return data[(i-1)*cols + (j-1)];}

        // Build a reference from any MATL matrix
        template<typename MAT>
        static Ref from(MAT const& mat)
        {
            Ref r(mat.rows(), mat.cols());
            for (matl::index_t i = 1; i <= r.rows; ++i)
                for (matl::index_t j = 1; j <= r.cols; ++j)
                    r(i,j) = static_cast<double>(mat(i,j));
            return r;
        }
    };

    inline Ref operator+(Ref const& a, Ref const& b)
    {Ref r(a.rows, a.cols); for (size_t k = 0; k < r.data.size(); ++k) r.data[k] = a.data[k] + b.data[k]; return r;}

    inline Ref operator-(Ref const& a, Ref const& b)
    {Ref r(a.rows, a.cols); for (size_t k = 0; k < r.data.size(); ++k) r.data[k] = a.data[k] - b.data[k]; return r;}

    inline Ref operator*(Ref const& a, Ref const& b)
    {
        Ref r(a.rows, b.cols);
        for (matl::index_t i = 1; i <= a.rows; ++i)
            for (matl::index_t j = 1; j <= b.cols; ++j)
            {
                double sum = 0.0;
                for (matl::index_t k = 1; k <= a.cols; ++k) sum += a(i,k) * b(k,j);
                r(i,j) = sum;
            }
        return r;
    }

    inline Ref operator*(Ref const& a, double s) {Ref r = a; for (double& x : r.data) x *= s; return r;}
    inline Ref operator*(double s, Ref const& a) {return a * s;}
    inline Ref operator+(Ref const& a, double s) {Ref r = a; for (double& x : r.data) x += s; return r;}
    inline Ref operator+(double s, Ref const& a) {return a + s;}
    inline Ref operator-(Ref const& a, double s) {Ref r = a; for (double& x : r.data) x -= s; return r;}
    inline Ref operator-(double s, Ref const& a) {Ref r = a; for (double& x : r.data) x = s - x; return r;}
    inline Ref operator/(Ref const& a, double s) {Ref r = a; for (double& x : r.data) x /= s; return r;}

    inline Ref transpose(Ref const& a)
    {
        Ref r(a.cols, a.rows);
        for (matl::index_t i = 1; i <= a.rows; ++i)
            for (matl::index_t j = 1; j <= a.cols; ++j)
                r(j,i) = a(i,j);
        return r;
    }

    // Compare a MATL matrix against a reference. Small integers are used throughout the suite so
    // results are exact, but a tolerance is applied to allow for division.
    template<typename MAT>
    bool matrix_matches(MAT const& mat, Ref const& ref, std::string& why, double tol = 1e-12)
    {
        if (mat.rows() != ref.rows || mat.cols() != ref.cols)
        {
            why = "dimensions " + to_string(mat.rows()) + "x" + to_string(mat.cols())
                + " != expected " + to_string(ref.rows) + "x" + to_string(ref.cols);
            return false;
        }
        for (matl::index_t i = 1; i <= ref.rows; ++i)
            for (matl::index_t j = 1; j <= ref.cols; ++j)
                if (std::fabs(static_cast<double>(mat(i,j)) - ref(i,j)) > tol)
                {
                    why = "element (" + to_string(i) + "," + to_string(j) + ") = "
                        + to_string(static_cast<double>(mat(i,j))) + ", expected " + to_string(ref(i,j));
                    return false;
                }
        return true;
    }

    /********************************** Data generation *************************************/

    inline std::mt19937& rng()
    {
        static std::mt19937 engine(20261020u); // fixed seed: failures are reproducible
        return engine;
    }

    // Fill a MATL matrix with small random integers in [lo, hi] so arithmetic is exact
    template<typename MAT>
    void fill_random(MAT& mat, int lo = -5, int hi = 5)
    {
        std::uniform_int_distribution<int> dist(lo, hi);
        for (matl::index_t i = 1; i <= mat.size(); ++i) mat(i) = static_cast<typename MAT::dtype>(dist(rng()));
    }

    // Fill a MATL matrix with 1, 2, 3, ... in row-major order
    template<typename MAT>
    void fill_sequence(MAT& mat, double start = 1.0)
    {
        for (matl::index_t i = 1; i <= mat.size(); ++i) mat(i) = static_cast<typename MAT::dtype>(start + double(i - 1));
    }

    inline matl::index_t random_dim(matl::index_t lo = 1, matl::index_t hi = 12)
    {
        std::uniform_int_distribution<matl::index_t> dist(lo, hi);
        return dist(rng());
    }
}

/***************************************** Macros *****************************************/

#define MATL_TEST_CONCAT_(a, b) a##b
#define MATL_TEST_CONCAT(a, b)  MATL_TEST_CONCAT_(a, b)

#define TEST_CASE_IMPL_(name, fn)                                                               \
    static void fn();                                                                           \
    static matl_test::Registrar MATL_TEST_CONCAT(fn, _registrar)(name, __FILE__, __LINE__, fn); \
    static void fn()

#define TEST_CASE(name) TEST_CASE_IMPL_(name, MATL_TEST_CONCAT(matl_test_case_, __LINE__))

#define CHECK(cond)                                                                             \
    do { if (!(cond)) matl_test::report_failure(__FILE__, __LINE__, "CHECK(" #cond ")"); } while (0)

#define REQUIRE(cond)                                                                           \
    do { if (!(cond)) { matl_test::report_failure(__FILE__, __LINE__, "REQUIRE(" #cond ")");    \
                        throw matl_test::RequireFailure(); } } while (0)

#define CHECK_EQ(a, b)                                                                          \
    do { auto&& matl_va_ = (a); auto&& matl_vb_ = (b);                                          \
         if (!(matl_va_ == matl_vb_))                                                           \
             matl_test::report_failure(__FILE__, __LINE__, "CHECK_EQ(" #a ", " #b "): "         \
                 + matl_test::to_string(matl_va_) + " != " + matl_test::to_string(matl_vb_));   \
    } while (0)

#define CHECK_THROWS(expr)                                                                      \
    do { bool matl_threw_ = false; try { (void)(expr); } catch (...) { matl_threw_ = true; }    \
         if (!matl_threw_) matl_test::report_failure(__FILE__, __LINE__, "CHECK_THROWS(" #expr ")"); \
    } while (0)

#define CHECK_THROWS_AS(expr, type)                                                             \
    do { bool matl_threw_ = false;                                                              \
         try { (void)(expr); } catch (type const&) { matl_threw_ = true; } catch (...) {}       \
         if (!matl_threw_)                                                                      \
             matl_test::report_failure(__FILE__, __LINE__, "CHECK_THROWS_AS(" #expr ", " #type ")"); \
    } while (0)

#define CHECK_ABORTS(...)                                                                       \
    do { if (!matl_test::terminates([&]() { __VA_ARGS__; }))                                    \
             matl_test::report_failure(__FILE__, __LINE__, "CHECK_ABORTS(" #__VA_ARGS__ ")");   \
    } while (0)

#define CHECK_MATRIX_EQ(mat, ref)                                                               \
    do { std::string matl_why_;                                                                 \
         if (!matl_test::matrix_matches((mat), (ref), matl_why_))                               \
             matl_test::report_failure(__FILE__, __LINE__,                                      \
                 "CHECK_MATRIX_EQ(" #mat ", " #ref "): " + matl_why_);                          \
    } while (0)

#endif
