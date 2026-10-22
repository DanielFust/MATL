# MATL test suite (see tests/README.md). Run with 'ctest' from the build directory.
file(GLOB MATL_TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# The suite with MATL's run-time checks enabled, whatever the build type (a Release build defines
# NDEBUG, which would otherwise disable them)
add_executable(matl_tests ${MATL_TEST_SOURCES})
target_link_libraries(matl_tests PRIVATE matl::matl matl_warnings)
target_include_directories(matl_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(matl_tests PRIVATE MATL_BOUNDS_CHECK MATL_DIMENSION_CHECK)
add_test(NAME matl_tests COMMAND matl_tests)

# The same suite with errors that abort instead of throwing
add_executable(matl_tests_abort ${MATL_TEST_SOURCES})
target_link_libraries(matl_tests_abort PRIVATE matl::matl matl_warnings)
target_include_directories(matl_tests_abort PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(matl_tests_abort PRIVATE MATL_BOUNDS_CHECK MATL_DIMENSION_CHECK MATL_ABORT_ON_ERROR)
add_test(NAME matl_tests_abort COMMAND matl_tests_abort)

# MATL compiled without exceptions (errors abort)
if(NOT MSVC)
    add_executable(matl_no_exceptions no_exceptions/main.cpp)
    target_link_libraries(matl_no_exceptions PRIVATE matl::matl matl_warnings)
    target_compile_definitions(matl_no_exceptions PRIVATE MATL_BOUNDS_CHECK MATL_DIMENSION_CHECK)
    target_compile_options(matl_no_exceptions PRIVATE -fno-exceptions)
    add_test(NAME matl_no_exceptions COMMAND matl_no_exceptions)
endif()

# The same suite with the run-time checks compiled out
add_executable(matl_tests_nodebug ${MATL_TEST_SOURCES})
target_link_libraries(matl_tests_nodebug PRIVATE matl::matl matl_warnings)
target_include_directories(matl_tests_nodebug PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(matl_tests_nodebug PRIVATE MATL_NO_DEBUG)
add_test(NAME matl_tests_nodebug COMMAND matl_tests_nodebug)
