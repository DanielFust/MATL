/**************************************************************************************************
 * File: tests/no_exceptions/main.cpp
 *
 * Description: Checks that MATL compiles without exceptions (-fno-exceptions) and that an error
 *              then aborts (see _matl_errors_.hpp). Built and run by 'make check-noexceptions'
 *              and by CTest. Exits with 0 if all checks pass.
 **************************************************************************************************/
#include <cstdio>
#include <sys/wait.h>
#include <unistd.h>

#include "matl.hpp"

using namespace matl;

static_assert(MATL_THROWS_ON_ERROR == 0, "built without exceptions, MATL errors abort");

int main()
{
    // Ordinary use works
    Matrix<double> A = eye<double>(20), B(20, 20), C(20, 20);
    B = 2.0;
    C = A * B + transpose(B) - 1.0;
    if (C(1, 1) != 3.0 || C(20, 1) != 3.0) {std::puts("wrong result"); return 1;}      // 2 + 2 - 1

    // An error aborts: incompatible operands, in a child process
    std::fflush(nullptr);
    const pid_t pid = fork();
    if (pid == 0)
    {
        Matrix<double> D(3, 3), R(20, 20);
        D = 1.0;
        R = A + D;                       // 20x20 + 3x3
        _exit(0);                        // not reached if the error aborts
    }
    int status = 0;
    waitpid(pid, &status, 0);
    const bool aborted = WIFSIGNALED(status);
    std::puts(aborted ? "no-exceptions build: ok" : "no-exceptions build: the error did not abort");
    return aborted ? 0 : 1;
}
