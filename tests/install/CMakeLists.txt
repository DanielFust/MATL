# Checks that an installed MATL can be used with find_package (run by CI):
#     cmake --install <build> --prefix <prefix>
#     cmake -S tests/install -B <dir> -DCMAKE_PREFIX_PATH=<prefix>
#     cmake --build <dir> && <dir>/matl_install_test
cmake_minimum_required(VERSION 3.14)
project(matl_install_test LANGUAGES CXX)

find_package(matl 2 REQUIRED)

add_executable(matl_install_test main.cpp)
target_link_libraries(matl_install_test PRIVATE matl::matl)
