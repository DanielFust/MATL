/**************************************************************************************************
 * File: test_main.cpp
 *
 * Description: Test runner for the MATL test suite.
 *
 * Usage:
 *     matl_tests                 Run all test cases
 *     matl_tests <substring>     Run test cases whose name contains <substring>
 *     matl_tests --list          List test case names
 *
 * Results are written to stderr. MATL's MATL_DEBUG trace output (enabled in the sanitizer build)
 * goes to stdout, so it can be redirected separately (the Makefile writes it to a log file in
 * build/).
 **************************************************************************************************/
#include <cstring>

#include "matl_test.hpp"

int main(int argc, char** argv)
{
    const char* filter = (argc > 1) ? argv[1] : nullptr;

    if (filter != nullptr && std::strcmp(filter, "--list") == 0)
    {
        for (const matl_test::TestCase& test : matl_test::registry()) std::printf("%s\n", test.name);
        return 0;
    }

    int run = 0, failed = 0;
    for (const matl_test::TestCase& test : matl_test::registry())
    {
        if (filter != nullptr && std::strstr(test.name, filter) == nullptr) continue;

        ++run;
        matl_test::current_failures() = 0;

        // Announce the test first, so a crash (e.g. a segfault) identifies the test case
        std::fprintf(stderr, "[ RUN    ] %s\n", test.name);
        try
        {
            test.body();
        }
        catch (matl_test::RequireFailure const&)
        {
            // Failure already reported by REQUIRE
        }
        catch (std::exception const& e)
        {
            matl_test::report_failure(test.file, test.line, std::string("unexpected exception: ") + e.what());
        }
        catch (...)
        {
            matl_test::report_failure(test.file, test.line, "unexpected non-standard exception");
        }

        if (matl_test::current_failures() == 0)
            std::fprintf(stderr, "[  PASS  ] %s\n", test.name);
        else
        {
            std::fprintf(stderr, "[  FAIL  ] %s (%d failure%s)\n", test.name,
                         matl_test::current_failures(), matl_test::current_failures() == 1 ? "" : "s");
            ++failed;
        }
    }

    std::fprintf(stderr, "\n%d test case%s run, %d passed, %d failed\n",
                 run, run == 1 ? "" : "s", run - failed, failed);
    if (run == 0)
    {
        std::fprintf(stderr, "No test cases matched the filter\n");
        return 1;
    }
    return failed == 0 ? 0 : 1;
}
