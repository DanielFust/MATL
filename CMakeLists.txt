# MATL (Matrix Algebra Template Library): header-only, C++20.
#
# Use from another CMake project, either after installing MATL:
#     find_package(matl 2 REQUIRED)
#     target_link_libraries(my_program PRIVATE matl::matl)
# or from a copy of the source tree (e.g. a git submodule or FetchContent):
#     add_subdirectory(MATL)
#     target_link_libraries(my_program PRIVATE matl::matl)
cmake_minimum_required(VERSION 3.14)

# The version is defined once, in include/matl.hpp
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/matl.hpp" matl_version_lines
     REGEX "^#define MATL_VERSION_(MAJOR|MINOR|PATCH) [0-9]+$")
foreach(line IN LISTS matl_version_lines)
    string(REGEX MATCH "^#define MATL_VERSION_(MAJOR|MINOR|PATCH) ([0-9]+)$" _ "${line}")
    set(MATL_VERSION_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
endforeach()

project(matl
    VERSION ${MATL_VERSION_MAJOR}.${MATL_VERSION_MINOR}.${MATL_VERSION_PATCH}
    DESCRIPTION "Matrix Algebra Template Library"
    LANGUAGES CXX)

# Tests, examples and install rules are on by default only when MATL is the top-level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(MATL_IS_TOP_LEVEL ON)
else()
    set(MATL_IS_TOP_LEVEL OFF)
endif()

option(MATL_BUILD_TESTS          "Build the MATL test suite"                ${MATL_IS_TOP_LEVEL})
option(MATL_BUILD_EXAMPLES       "Build the MATL example program"           OFF)
option(MATL_BUILD_BENCHMARKS     "Build the MATL benchmarks"                OFF)
option(MATL_INSTALL              "Generate MATL's install rules"            ${MATL_IS_TOP_LEVEL})
option(MATL_WARNINGS_AS_ERRORS   "Treat compiler warnings as errors"        OFF)

# Optimized builds by default (Release defines NDEBUG, which disables MATL's run-time checks)
if(MATL_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(GNUInstallDirs)

# The library
add_library(matl INTERFACE)
add_library(matl::matl ALIAS matl)
target_include_directories(matl INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(matl INTERFACE cxx_std_20)

# Warning flags for MATL's own programs (tests, examples, benchmarks)
add_library(matl_warnings INTERFACE)
if(MSVC)
    target_compile_options(matl_warnings INTERFACE /W4 $<$<BOOL:${MATL_WARNINGS_AS_ERRORS}>:/WX>)
else()
    target_compile_options(matl_warnings INTERFACE -Wall -Wextra $<$<BOOL:${MATL_WARNINGS_AS_ERRORS}>:-Werror>)
endif()

if(MATL_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(MATL_BUILD_EXAMPLES)
    add_executable(matl_example examples/test.cpp)
    target_link_libraries(matl_example PRIVATE matl::matl matl_warnings)
endif()

if(MATL_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# Installation: headers, and a CMake package for find_package(matl)
if(MATL_INSTALL)
    include(CMakePackageConfigHelpers)
    set(MATL_CMAKE_DIR "${CMAKE_INSTALL_DATADIR}/cmake/matl")

    install(TARGETS matl EXPORT matlTargets)
    install(DIRECTORY include/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
    install(EXPORT matlTargets NAMESPACE matl:: DESTINATION "${MATL_CMAKE_DIR}")

    configure_package_config_file(cmake/matlConfig.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/matlConfig.cmake"
        INSTALL_DESTINATION "${MATL_CMAKE_DIR}")
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/matlConfigVersion.cmake"
        COMPATIBILITY SameMajorVersion
        ARCH_INDEPENDENT)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/matlConfig.cmake"
                  "${CMAKE_CURRENT_BINARY_DIR}/matlConfigVersion.cmake"
            DESTINATION "${MATL_CMAKE_DIR}")
endif()
