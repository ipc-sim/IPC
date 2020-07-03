// Eigen wrappers for the exact CCD methods.

#include "ExactCCD.hpp"

#ifdef USE_EXACT_CCD
#include <bsc.h>
#include <rootparitycollisiontest.h>
#include "ECCD.hpp"
#endif

namespace IPC {
namespace ExactCCD {

// Detect collisions between two edges as they move.
bool edgeEdgeCCD(
    const Eigen::Vector3d& q0start,
    const Eigen::Vector3d& p0start,
    const Eigen::Vector3d& q1start,
    const Eigen::Vector3d& p1start,
    const Eigen::Vector3d& q0end,
    const Eigen::Vector3d& p0end,
    const Eigen::Vector3d& q1end,
    const Eigen::Vector3d& p1end,
    const Method method)
{
    switch (method) {
#ifdef USE_EXACT_CCD
    case Method::ROOT_PARITY:
        return rootparity::RootParityCollisionTest(
            // Edge 1 at t=0
            Vec3d(q0start.data()), Vec3d(p0start.data()),
            // Edge 2 at t=0
            Vec3d(q1start.data()), Vec3d(p1start.data()),
            // Edge 1 at t=1
            Vec3d(q0end.data()), Vec3d(p0end.data()),
            // Edge 2 at t=1
            Vec3d(q1end.data()), Vec3d(p1end.data()),
            /* is_edge_edge = */ true)
            .run_test();
    case Method::BSC:
        return bsc::Intersect_EE_robust(
            // Edge 1 at t=0
            Vec3d(q0start.data()), Vec3d(p0start.data()),
            // Edge 2 at t=0
            Vec3d(q1start.data()), Vec3d(p1start.data()),
            // Edge 1 at t=1
            Vec3d(q0end.data()), Vec3d(p0end.data()),
            // Edge 2 at t=1
            Vec3d(q1end.data()), Vec3d(p1end.data()));
    case Method::RATIONAL_ROOT_PARITY:
        return eccd::edgeEdgeCCD(q0start, p0start, q1start, p1start,
            q0end, p0end, q1end, p1end);
#endif
    default:
        return false;
    }
}

// Detect collisions between a vertex and a triangular face.
bool vertexFaceCCD(
    const Eigen::Vector3d& q0start,
    const Eigen::Vector3d& q1start,
    const Eigen::Vector3d& q2start,
    const Eigen::Vector3d& q3start,
    const Eigen::Vector3d& q0end,
    const Eigen::Vector3d& q1end,
    const Eigen::Vector3d& q2end,
    const Eigen::Vector3d& q3end,
    const Method method)
{
    switch (method) {
#ifdef USE_EXACT_CCD
    case Method::ROOT_PARITY:
        return rootparity::RootParityCollisionTest(
            // Point at t=0
            Vec3d(q0start.data()),
            // Triangle at t = 0
            Vec3d(q2start.data()), Vec3d(q1start.data()), Vec3d(q3start.data()),
            // Point at t=1
            Vec3d(q0end.data()),
            // Triangle at t = 1
            Vec3d(q2end.data()), Vec3d(q1end.data()), Vec3d(q3end.data()),
            /* is_edge_edge = */ false)
            .run_test();
        // rootparity::RootParityCollisionTest(
        // // Point at t=0
        // Vec3d(q0start.data()),
        // // Triangle at t = 0
        // Vec3d(q1start.data()), Vec3d(q2start.data()), Vec3d(q3start.data()),
        // // Point at t=1
        // Vec3d(q0end.data()),
        // // Triangle at t = 1
        // Vec3d(q1end.data()), Vec3d(q2end.data()), Vec3d(q3end.data()),
        // /* is_edge_edge = */ false)
        // .run_test();
    case Method::BSC:
        return bsc::Intersect_VF_robust(
            // Point at t=0
            Vec3d(q0start.data()),
            // Triangle at t = 0
            Vec3d(q1start.data()), Vec3d(q2start.data()), Vec3d(q3start.data()),
            // Point at t=1
            Vec3d(q0end.data()),
            // Triangle at t = 1
            Vec3d(q1end.data()), Vec3d(q2end.data()), Vec3d(q3end.data()));
    case Method::RATIONAL_ROOT_PARITY:
        return eccd::vertexFaceCCD(q0start, q1start, q2start, q3start,
            q0end, q1end, q2end, q3end);
#endif
    default:
        return false;
    }
}

}
} // namespace IPC::ExactCCD
