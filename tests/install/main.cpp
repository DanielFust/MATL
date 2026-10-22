// Uses an installed MATL through find_package(matl) (see CMakeLists.txt)
#include <matl.hpp>

static_assert(MATL_VERSION_MAJOR == 2, "found the expected MATL major version");

int main()
{
    matl::Matrix<double> A = matl::eye<double>(3), B(3, 3);
    B = 2.0 * A + A * A;
    return (B(1,1) == 3.0 && B(1,2) == 0.0) ? 0 : 1;
}
