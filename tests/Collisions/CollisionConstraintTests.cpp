#include <catch2/catch.hpp>

#include "CollisionConstraints.hpp"
#include "CTCD.h"

#include <Eigen/Dense>
#include <finitediff.hpp>

const std::vector<std::string> constraintNames = {
    "VOLUME", "GRAPHICS", "NONSMOOTH_NEWMARK", "GAP_FUNCTION", "CMR",
    "VERSCHOOR"
};

TEST_CASE("Test Point-Triangle Collision Constraints",
    "[collision-constraints][sqp]")
{
    // point
    Eigen::Vector3d v0(0, 1, -0.5);
    // triangle = (v1, v2, v3) in counter-clockwise order
    Eigen::Vector3d v1(-1, 0, 1);
    Eigen::Vector3d v2(1, 0, 1);
    Eigen::Vector3d v3(0, 0, -1);

    // displacements
    double u0y = GENERATE(-1.1, 0.0, 1.1);
    Eigen::Vector3d u0(0, u0y, 0);
    double u1y = GENERATE(-1.1, 0.0, 1.1);
    Eigen::Vector3d u1(0, u1y, 0);

    bool intersects = (-u0y + u1y >= 1);

    double toi;
    bool ccd_intersects = CTCD::vertexFaceCTCD(
        v0, v1, v3, v2, v0 + u0, v1 + u1, v3 + u1, v2 + u1, 0, toi);
    REQUIRE(ccd_intersects == intersects);

    IPC::CollisionConstraintType constraintType = GENERATE(
        IPC::CollisionConstraintType::VOLUME,
        IPC::CollisionConstraintType::GRAPHICS,
        IPC::CollisionConstraintType::VERSCHOOR);
    double c;
    IPC::compute_collision_constraint(
        v0, v1, v3, v2,
        v0 + u0, v1 + u1, v3 + u1, v2 + u1,
        constraintType, /*is_edge_edge=*/false, toi,
        c);

    CAPTURE(u0y, u1y, c, constraintNames[constraintType]);
    CHECK((c < 0) == intersects);

    Eigen::Vector12d grad_c;
    IPC::compute_collision_constraint_gradient(
        v0, v1, v3, v2,
        v0 + u0, v1 + u1, v3 + u1, v2 + u1,
        constraintType, /*is_edge_edge=*/false, toi,
        grad_c);

    Eigen::VectorXd x(12);
    x.segment<3>(0) = v0 + u0;
    x.segment<3>(3) = v1 + u1;
    x.segment<3>(6) = v3 + u1; // Swap order to match compute_collision_constraint
    x.segment<3>(9) = v2 + u1; // Swap order to match compute_collision_constraint
    std::function<double(const Eigen::VectorXd&)> f = [&](const Eigen::VectorXd& x) {
        assert(x.size() == 12);
        double c;
        IPC::compute_collision_constraint(
            v0, v1, v3, v2,
            x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9),
            constraintType, /*is_edge_edge=*/false, toi,
            c);
        return c;
    };
    Eigen::VectorXd finite_grad_c(12);
    fd::finite_gradient(x, f, finite_grad_c);
    CAPTURE(grad_c.transpose(), finite_grad_c.transpose());
    CHECK(fd::compare_gradient(grad_c, finite_grad_c));
}

TEST_CASE("Test Edge-Edge Collision Constraints",
    "[collision-constraints][sqp]")
{
    // e0 = (v0, v1)
    Eigen::Vector3d v0(-1, -1, 0);
    Eigen::Vector3d v1(1, -1, 0);
    // e2 = (v2, v3)
    Eigen::Vector3d v2(0, 1, -1);
    Eigen::Vector3d v3(0, 1, 1);

    // displacements
    double y_displacement = GENERATE(-2, 0, 2);
    Eigen::Vector3d u0(0, y_displacement, 0);
    Eigen::Vector3d u1(0, -y_displacement, 0);

    bool intersects = y_displacement >= 1.0;
    double toi;
    bool ccd_intersects = CTCD::edgeEdgeCTCD(
        v0, v1, v2, v3, v0 + u0, v1 + u0, v2 + u1, v3 + u1, 0, toi);
    CAPTURE(y_displacement);
    REQUIRE(ccd_intersects == intersects);

    IPC::CollisionConstraintType constraintType = GENERATE(
        IPC::CollisionConstraintType::VOLUME,
        IPC::CollisionConstraintType::GRAPHICS,
        IPC::CollisionConstraintType::VERSCHOOR);

    double c;
    IPC::compute_collision_constraint(
        v0, v1, v2, v3,
        v0 + u0, v1 + u0, v2 + u1, v3 + u1,
        constraintType, /*is_edge_edge=*/true, toi,
        c);

    CAPTURE(y_displacement, c, constraintNames[constraintType]);
    CHECK((c < 0) == intersects);

    typedef Eigen::Matrix<double, 12, 1> Vector12d;
    Vector12d grad_c;
    IPC::compute_collision_constraint_gradient(
        v0, v1, v2, v3,
        v0 + u0, v1 + u0, v2 + u1, v3 + u1,
        constraintType, /*is_edge_edge=*/true, toi,
        grad_c);

    Eigen::VectorXd x(12);
    x.segment<3>(0) = v0 + u0;
    x.segment<3>(3) = v1 + u0;
    x.segment<3>(6) = v2 + u1;
    x.segment<3>(9) = v3 + u1;
    std::function<double(const Eigen::VectorXd&)> f = [&](const Eigen::VectorXd& x) {
        assert(x.size() == 12);
        double c;
        IPC::compute_collision_constraint(
            v0, v1, v2, v3,
            x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9),
            constraintType, /*is_edge_edge=*/true, toi,
            c);
        return c;
    };
    Eigen::VectorXd finite_grad_c(12);
    fd::finite_gradient(x, f, finite_grad_c);
    CAPTURE(grad_c.transpose(), finite_grad_c.transpose());
    CHECK(fd::compare_gradient(grad_c, finite_grad_c));
}

TEST_CASE("Barycentric coordinates", "[barycentric]")
{
    // Create random triangle
    Eigen::Vector3d a = Eigen::Vector3d::Random();
    Eigen::Vector3d b = Eigen::Vector3d::Random();
    Eigen::Vector3d c = Eigen::Vector3d::Random();
    // Create random point
    Eigen::Vector3d barycentric_coords = Eigen::Vector3d::Random();
    barycentric_coords(0) = 1 - barycentric_coords(1) - barycentric_coords(2);
    Eigen::Vector3d p = barycentric_coords(0) * a + barycentric_coords(1) * b + barycentric_coords(2) * c;

    Eigen::Vector3d computed_barycentric_coords;
    IPC::barycentric_coordinates(p.transpose(), a.transpose(), b.transpose(), c.transpose(), computed_barycentric_coords);
    Eigen::Vector3d computed_p = computed_barycentric_coords(0) * a + computed_barycentric_coords(1) * b + computed_barycentric_coords(2) * c;

    CAPTURE(barycentric_coords, computed_barycentric_coords);
    CHECK((computed_barycentric_coords - barycentric_coords).norm() < 1e-12);
    CAPTURE(p, computed_p);
    CHECK((computed_p - p).norm() < 1e-12);
}
