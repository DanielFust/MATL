/**************************************************************************************************
 * Example 10: determinants, inverses, eigenvalues and the SVD
 *
 * Four small problems from mechanics and data fitting:
 *
 *   1. The Jacobian of a finite-element mapping: det and inv of a fixed-size 2 x 2 matrix, as at
 *      every integration point of an element (closed forms, no heap memory).
 *   2. Principal stresses: the eigenvalues and eigenvectors of a symmetric 3 x 3 stress tensor
 *      (eig_sym).
 *   3. Natural frequencies of a spring-mass chain: the generalized problem K v = w^2 M v, turned
 *      into a symmetric one with M^(-1/2).
 *   4. A rank-deficient least-squares fit: the SVD reveals the rank, and the pseudo-inverse gives
 *      the smallest solution (pinv, svd(A).solve).
 *
 * Build and run from the repository root:
 *
 *     make examples && ./build/examples/10_eigenvalues_and_svd
 **************************************************************************************************/
#include <cmath>
#include <iostream>
#include <matl.hpp>

using namespace matl;

int main()
{
    /******************************* 1. An element Jacobian *************************************/

    // A bilinear quadrilateral maps the reference square (xi, eta) in [-1, 1]^2 to the element.
    // Corner coordinates (one row per node), and the shape-function derivatives at the center
    // (xi = eta = 0): row k holds (dN_k/dxi, dN_k/deta).
    const Matrix<double, 4, 2> nodes = {{0.0, 0.0}, {2.0, 0.0}, {2.5, 1.5}, {0.5, 1.0}};
    const Matrix<double, 4, 2> dN_dxi = {{-0.25, -0.25}, {0.25, -0.25}, {0.25, 0.25}, {-0.25, 0.25}};

    // J = dx/dxi (2 x 2). det(J) scales areas (it must be positive for a valid element), and
    // inv(J) maps the derivatives to physical coordinates: dN/dx = dN/dxi * J^-1.
    const Matrix<double, 2, 2> J = transpose(nodes) * dN_dxi;
    const double detJ = det(J);
    const Matrix<double, 2, 2> Jinv = inv(J);                  // fixed-size: closed form
    const Matrix<double, 4, 2> dN_dx = dN_dxi * Jinv;
    std::cout << "J =\n" << J << "\ndet(J) = " << detJ
              << " (det(J) is linear in xi and eta, so the element's area is 4 det(J) at the center: "
              << 4 * detJ << ")\n";
    std::cout << "dN/dx at the center:\n" << dN_dx << "\n";
    std::cout << "sum of dN/dx over the nodes (zero: the shape functions sum to 1):\n"
              << sum(dN_dx(all, 1)) << " " << sum(dN_dx(all, 2)) << "\n\n";

    /******************************* 2. Principal stresses **************************************/

    // A stress tensor (MPa) is symmetric; its eigenvalues are the principal stresses, and its
    // eigenvectors the principal directions (orthonormal)
    const Matrix<double, 3, 3> sigma = {{50, 30, 20}, {30, -20, -10}, {20, -10, 10}};
    const auto P = eig_sym(sigma);                            // ascending: sigma_3 <= sigma_2 <= sigma_1
    const double s3 = P.values()(1), s2 = P.values()(2), s1 = P.values()(3);
    std::cout << "principal stresses: " << s1 << ", " << s2 << ", " << s3 << " MPa\n";
    std::cout << "direction of the largest:\n" << P.vectors()(all, 3) << "\n";

    // Invariants: the trace and the determinant are the sum and product of the eigenvalues
    std::cout << "trace " << trace(sigma) << " = " << s1 + s2 + s3
              << ",  det " << det(sigma) << " = " << s1 * s2 * s3 << "\n";

    // The von Mises stress from the principal stresses, and the largest shear stress
    const double von_mises = std::sqrt(0.5 * ((s1 - s2) * (s1 - s2) + (s2 - s3) * (s2 - s3) + (s3 - s1) * (s3 - s1)));
    std::cout << "von Mises stress " << von_mises << " MPa, largest shear " << 0.5 * (s1 - s3) << " MPa\n\n";

    /*************************** 3. Natural frequencies of a chain *******************************/

    // Three masses between two walls, connected by four springs: K v = w^2 M v. With M diagonal,
    // the substitution v = M^(-1/2) u gives the symmetric problem (M^(-1/2) K M^(-1/2)) u = w^2 u.
    const double k = 1000.0;                                  // N/m
    const Matrix<double, 3, 1> m = {{1.0}, {2.0}, {1.0}};     // kg
    const Matrix<double, 3, 3> K = {{2 * k, -k, 0}, {-k, 2 * k, -k}, {0, -k, 2 * k}};
    Matrix<double, 3, 3> A;
    for (index_t i = 1; i <= 3; ++i)
        for (index_t j = 1; j <= 3; ++j) A(i, j) = K(i, j) / std::sqrt(m(i) * m(j));
    const auto modes = eig_sym(A);
    std::cout << "natural frequencies (rad/s):";
    for (index_t i = 1; i <= 3; ++i) std::cout << " " << std::sqrt(modes.values()(i));
    Matrix<double, 3, 1> v = modes.vectors()(all, 1);          // the lowest mode, v = M^(-1/2) u
    for (index_t i = 1; i <= 3; ++i) v(i) /= std::sqrt(m(i));
    std::cout << "\nlowest mode shape (all masses move together):\n" << v / v(2) << "\n\n";

    /*************************** 4. A rank-deficient least-squares fit ***************************/

    // Fit y ~ c1 + c2 x + c3 (2x - 1): the third column is a combination of the first two, so the
    // design matrix has rank 2 and the coefficients are not unique. QR would find R singular;
    // the SVD shows the rank, and the pseudo-inverse picks the solution of smallest norm.
    const index_t n = 6;
    Matrix<double> V(n, 3), y(n, 1);
    for (index_t i = 1; i <= n; ++i)
    {
        const double x = double(i - 1) / double(n - 1);
        V(i, 1) = 1.0;  V(i, 2) = x;  V(i, 3) = 2.0 * x - 1.0;
        y(i) = 3.0 + 4.0 * x;
    }
    const auto S = svd(V);
    std::cout << "singular values of V:\n" << S.singular_values() << "\n";
    std::cout << "rank " << S.rank() << " (the default tolerance " << S.tolerance() << " counts the last as zero)\n";

    const Matrix<double> c = S.solve(y);                      // = pinv(V) * y
    std::cout << "smallest coefficients:\n" << c << "\n";
    std::cout << "residual ||V c - y|| = " << norm_fro(V * c - y) << "\n";

    // Every fit c + t (1, -2, 1) is as good (V maps that direction to zero); the pseudo-inverse's
    // solution has no component along it
    const Matrix<double> null_direction = S.V()(all, 3);
    std::cout << "component of c along the null space: " << dot(c, null_direction) << "\n";
    std::cout << "pinv(V) =\n" << pinv(V) << "\n";

    return 0;
}
