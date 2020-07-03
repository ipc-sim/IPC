#include <catch2/catch.hpp>

#ifdef USE_EXACT_CCD

#include "ExactCCD.hpp"

static const double EPSILON = std::numeric_limits<float>::epsilon();

TEST_CASE("Test Edge-Edge Exact Continous Collision Detection",
    "[ccd][exact-ccd][edge-edge]")
{
    // EE Exact CCD unit test
    // e0 = (v0, v1)
    Eigen::Vector3d v0(-1, -1, 0);
    Eigen::Vector3d v1(1, -1, 0);
    // e2 = (v2, v3)
    double e1x = GENERATE(
        -1 - EPSILON, -1, -1 + EPSILON, -0.5, 0, 0.5, 1 - EPSILON, 1,
        1 + EPSILON);
    Eigen::Vector3d v2(e1x, 1, -1);
    Eigen::Vector3d v3(e1x, 1, 1);

    // displacements
    double y_displacement = GENERATE(
        -1.0, 0.0, 1 - EPSILON, 1.0, 1 + EPSILON, 2.0);
    Eigen::Vector3d u0(0, y_displacement, 0);
    Eigen::Vector3d u1(0, -y_displacement, 0);

    IPC::ExactCCD::Method method = GENERATE(
        IPC::ExactCCD::Method::ROOT_PARITY,
        IPC::ExactCCD::Method::BSC);

    bool hit = IPC::ExactCCD::edgeEdgeCCD(
        v0, v1, v2, v3,
        v0 + u0, v1 + u0, v2 + u1, v3 + u1,
        method);

    const char* methodNames[3] = { "None", "RootParity", "BSC" };
    CAPTURE(y_displacement, e1x, methodNames[method]);
    CHECK(hit == (y_displacement >= 1.0 && e1x >= -1 && e1x <= 1));
}

TEST_CASE("Test Point-Triangle Exact Continous Collision Detection",
    "[ccd][exact-ccd][point-triangle]")
{
    // PT Exact CCD unit test
    // point
    double v0z = GENERATE(0.0, -1.0);
    Eigen::Vector3d v0(0, 1, v0z);
    // triangle = (v1, v2, v3)
    Eigen::Vector3d v1(-1, 0, 1);
    Eigen::Vector3d v2(1, 0, 1);
    Eigen::Vector3d v3(0, 0, -1);

    // displacements
    double u0y = -GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
    double u0z = GENERATE(-EPSILON, 0.0, EPSILON);
    Eigen::Vector3d u0(0, u0y, u0z);
    double u1y = GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
    Eigen::Vector3d u1(0, u1y, 0);

    IPC::ExactCCD::Method method = GENERATE(
        IPC::ExactCCD::Method::ROOT_PARITY,
        IPC::ExactCCD::Method::BSC);

    SECTION("Clockwise triangle")
    {
        // already in clockwise order
    }
    SECTION("Counter-clockwise triangle")
    {
        std::swap(v1, v2);
    }

    bool hit = IPC::ExactCCD::vertexFaceCCD(
        v0, v1, v2, v3,
        v0 + u0, v1 + u1, v2 + u1, v3 + u1,
        method);

    const char* methodNames[3] = { "None", "RootParity", "BSC" };
    CAPTURE(v0z, u0y, u1y, u0z, EPSILON, methodNames[method]);
    CHECK(hit == ((-u0y + u1y >= 1) && (v0z + u0z >= v3.z())));
}

#endif
