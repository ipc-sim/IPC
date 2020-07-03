/// @brief Constraint functions for use with constraint based collisions.

#include "CollisionConstraints.hpp"
#include "MeshCollisionUtils.hpp"

#include "CTCD.h"

#include <spdlog/spdlog.h>

namespace IPC {

void compute_collision_constraint(
    const Eigen::Vector3d& v0_t0, // First vertex at the start of time-step
    const Eigen::Vector3d& v1_t0, // Second vertex at the start of time-step
    const Eigen::Vector3d& v2_t0, // Third vertex at the start of time-step
    const Eigen::Vector3d& v3_t0, // Fourth vertex at the start of time-step
    const Eigen::Vector3d& v0_t1, // First vertex at the end of time-step
    const Eigen::Vector3d& v1_t1, // Second vertex at the end of time-step
    const Eigen::Vector3d& v2_t1, // Third vertex at the end of time-step
    const Eigen::Vector3d& v3_t1, // Fourth vertex at the end of time-step
    const CollisionConstraintType constraintType, bool is_edge_edge,
    double toi, double& c)
{
    switch (constraintType) {
    case CollisionConstraintType::NONSMOOTH_NEWMARK:
    case CollisionConstraintType::VOLUME:
        compute_collision_volume_constraint(
            v0_t1, v1_t1, v2_t1, v3_t1, c);
        break;
    case CollisionConstraintType::GAP_FUNCTION:
    case CollisionConstraintType::GRAPHICS:
        if (is_edge_edge) {
            compute_graphics_edge_edge_constraint(
                v0_t1, v1_t1, v2_t1, v3_t1, c);
        }
        else {
            // Swap order to counter-clockwise
            compute_graphics_point_triangle_constraint(
                v0_t1, v1_t1, v3_t1, v2_t1, c);
        }
        break;
    case CollisionConstraintType::CMR:
    case CollisionConstraintType::VERSCHOOR:
        if (is_edge_edge) {
            compute_Verschoor_edge_edge_constraint(
                v0_t0, v1_t0, v2_t0, v3_t0, v0_t1, v1_t1, v2_t1, v3_t1, toi, c);
        }
        else {
            // Swap order to counter-clockwise
            compute_Verschoor_point_triangle_constraint(
                v0_t0, v1_t0, v3_t0, v2_t0, v0_t1, v1_t1, v3_t1, v2_t1, toi, c);
        }
        break;
    case CollisionConstraintType::STIV:
        if (is_edge_edge) {
            compute_Verschoor_edge_edge_constraint(
                v0_t0, v1_t0, v2_t0, v3_t0, v0_t1, v1_t1, v2_t1, v3_t1, toi, c);
        }
        else {
            // Swap order to counter-clockwise
            compute_STIV_point_triangle_constraint(
                v0_t0, v1_t0, v3_t0, v2_t0, v0_t1, v1_t1, v3_t1, v2_t1, toi, c);
        }
        break;
    }
}

void compute_collision_constraint_gradient(
    const Eigen::Vector3d& v0_t0, // First vertex at the start of time-step
    const Eigen::Vector3d& v1_t0, // Second vertex at the start of time-step
    const Eigen::Vector3d& v2_t0, // Third vertex at the start of time-step
    const Eigen::Vector3d& v3_t0, // Fourth vertex at the start of time-step
    const Eigen::Vector3d& v0_t1, // First vertex at the end of time-step
    const Eigen::Vector3d& v1_t1, // Second vertex at the end of time-step
    const Eigen::Vector3d& v2_t1, // Third vertex at the end of time-step
    const Eigen::Vector3d& v3_t1, // Fourth vertex at the end of time-step
    const CollisionConstraintType constraintType, bool is_edge_edge,
    double toi, Eigen::Vector12d& grad_c)
{
    switch (constraintType) {
    case CollisionConstraintType::NONSMOOTH_NEWMARK:
    case CollisionConstraintType::VOLUME:
        compute_collision_volume_constraint_gradient(
            v0_t1, v1_t1, v2_t1, v3_t1, grad_c);
        break;
    case CollisionConstraintType::GAP_FUNCTION:
    case CollisionConstraintType::GRAPHICS:
        if (is_edge_edge) {
            compute_graphics_edge_edge_constraint_gradient(
                v0_t1, v1_t1, v2_t1, v3_t1, grad_c);
        }
        else {
            // Expects triangle vertices in counter-clockwise order.
            // Swap order to counter-clockwise.
            compute_graphics_point_triangle_constraint_gradient(
                v0_t1, v1_t1, v3_t1, v2_t1, grad_c);
            // Swap the gradient segments to the expected clockwise order
            Eigen::Vector3d tmp = grad_c.segment<3>(2 * 3);
            grad_c.segment<3>(2 * 3) = grad_c.segment<3>(3 * 3);
            grad_c.segment<3>(3 * 3) = tmp;
        }
        break;
    case CollisionConstraintType::CMR:
    case CollisionConstraintType::VERSCHOOR:
        if (is_edge_edge) {
            compute_Verschoor_edge_edge_constraint_gradient(
                v0_t0, v1_t0, v2_t0, v3_t0, v0_t1, v1_t1, v2_t1, v3_t1, toi,
                grad_c);
        }
        else {
            // Expects triangle vertices in counter-clockwise order.
            // Swap order to counter-clockwise.
            compute_Verschoor_point_triangle_constraint_gradient(
                v0_t0, v1_t0, v3_t0, v2_t0, v0_t1, v1_t1, v3_t1, v2_t1, toi,
                grad_c);
            // Swap the gradient segments to the expected clockwise order
            Eigen::Vector3d tmp = grad_c.segment<3>(2 * 3);
            grad_c.segment<3>(2 * 3) = grad_c.segment<3>(3 * 3);
            grad_c.segment<3>(3 * 3) = tmp;
        }
        break;
    case CollisionConstraintType::STIV:
        if (is_edge_edge) {
            compute_Verschoor_edge_edge_constraint_gradient(
                v0_t0, v1_t0, v2_t0, v3_t0, v0_t1, v1_t1, v2_t1, v3_t1, toi,
                grad_c);
        }
        else {
            // Expects triangle vertices in counter-clockwise order.
            // Swap order to counter-clockwise.
            compute_STIV_point_triangle_constraint_gradient(
                v0_t0, v1_t0, v3_t0, v2_t0, v0_t1, v1_t1, v3_t1, v2_t1, toi,
                grad_c);
            // Swap the gradient segments to the expected clockwise order
            Eigen::Vector3d tmp = grad_c.segment<3>(2 * 3);
            grad_c.segment<3>(2 * 3) = grad_c.segment<3>(3 * 3);
            grad_c.segment<3>(3 * 3) = tmp;
        }
        break;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Tretrahedron volume constraint

void compute_collision_volume_constraint(
    const Eigen::Vector3d& v0, // point
    const Eigen::Vector3d& v1, // triangle point 0
    const Eigen::Vector3d& v2, // triangle point 2
    const Eigen::Vector3d& v3, // triangle point 1
    double& c)
{
    c = (v3 - v0).dot((v1 - v0).cross(v2 - v0));
}

void compute_collision_volume_constraint_gradient(
    const Eigen::Vector3d& v0, // point
    const Eigen::Vector3d& v1, // triangle point 0
    const Eigen::Vector3d& v2, // triangle point 2
    const Eigen::Vector3d& v3, // triangle point 1
    Eigen::Vector12d& grad_c)
{
    grad_c.segment<3>(3) = (v2 - v0).cross(v3 - v0);
    grad_c.segment<3>(6) = (v3 - v0).cross(v1 - v0);
    grad_c.segment<3>(9) = (v1 - v0).cross(v2 - v0);
    grad_c.segment<3>(0) = -grad_c.segment<3>(3) - grad_c.segment<3>(6) - grad_c.segment<3>(9);
}

///////////////////////////////////////////////////////////////////////////////
// Robust Treatment of Simultaneous Collisions [Harmon et al. 2008]

static std::string toString(const Eigen::MatrixXd& mat)
{
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

// Computes the barycentric coordinates of any point in the plane containing
// the triangle (a, b, c).
void barycentric_coordinates(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c,
    Eigen::Vector3d& coords)
{
    const Eigen::Array3d v0 = b.array() - a.array();
    const Eigen::Array3d v1 = c.array() - a.array();
    const Eigen::Array3d v2 = p.array() - a.array();
    double d00 = (v0 * v0).sum();
    double d01 = (v0 * v1).sum();
    double d11 = (v1 * v1).sum();
    double d20 = (v2 * v0).sum();
    double d21 = (v2 * v1).sum();
    double denom = d00 * d11 - d01 * d01;
    if (denom == 0) {
        spdlog::error("Unable to compute barycentric coordinates!");
        spdlog::debug("p=({}) a=({}) b=({}) c=({})",
            toString(p.transpose()), toString(a.transpose()),
            toString(b.transpose()), toString(c.transpose()));
    }
    coords(1) = (d11 * d20 - d01 * d21) / denom;
    coords(2) = (d00 * d21 - d01 * d20) / denom;
    coords(0) = 1.0 - (coords(1) + coords(2));
}

void compute_graphics_point_triangle_constraint(
    const Eigen::Vector3d& v0, // point
    const Eigen::Vector3d& v1, // triangle point 0
    const Eigen::Vector3d& v2, // triangle point 1
    const Eigen::Vector3d& v3, // triangle point 2
    double& c)
{
    // 1. Project point onto the plane containing the triangle.
    Eigen::Vector3d normal = (v2 - v1).cross(v3 - v1).normalized();
    c = normal.dot(v0 - v1); // distance to plane
}

void compute_graphics_edge_edge_constraint(
    const Eigen::Vector3d& v0, // first edge's first vertex
    const Eigen::Vector3d& v1, // first edge's second vertex
    const Eigen::Vector3d& v2, // second edge's first vertex
    const Eigen::Vector3d& v3, // second edge's second vertex
    double& c)
{
    // 1. Find closest points between e0 and e1 (see: https://bit.ly/2qwI33t)
    Eigen::Vector3d dir0 = v1 - v0;
    Eigen::Vector3d dir1 = v3 - v2;
    Eigen::Vector3d dir2 = dir1.cross(dir0);
    // v0 + t0 * dir0 + t2 * dir2 = v2 + t1 * dir1
    // [dir0, -dir1, dir2] * [t0; t1; t2] = v2 - v0
    // TODO: Make this more efficient.
    Eigen::Matrix3d A;
    A.col(0) = dir0;
    A.col(1) = -dir1;
    A.col(2) = dir2;
    Eigen::Vector3d params = A.lu().solve(v2 - v0);
    // Check the solution is valid
    if (!(params.array().isFinite().all())) {
        spdlog::error("Unable to compute edge-edge graphics constraint!");
        spdlog::debug("Linear solve for parameters along edges resulted in ({}); e₁×e₂=({})",
            toString(params.transpose()), toString(dir2.transpose()));
        c = 1e28;
        return;
    }
    assert((A * params - (v2 - v0)).squaredNorm() < 1e-12);
    // Project the points back onto the edges.
    Eigen::Vector3d closest_point0 = dir0 * std::clamp(params(0), 0.0, 1.0) + v0;
    Eigen::Vector3d closest_point1 = dir1 * std::clamp(params(1), 0.0, 1.0) + v2;
    // 2. Compute the distance from the points doted with the normal.
    Eigen::Vector3d normal = (v3 - v2).cross(v1 - v0).normalized();
    c = normal.dot(closest_point1 - closest_point0);
}

void compute_graphics_point_triangle_constraint_gradient(
    const Eigen::Vector3d& v0, // point
    const Eigen::Vector3d& v1, // triangle point 0
    const Eigen::Vector3d& v2, // triangle point 1
    const Eigen::Vector3d& v3, // triangle point 2
    Eigen::Vector12d& grad_c)
{
    // 1. Project point on to plane
    Eigen::Vector3d normal = (v2 - v1).cross(v3 - v1).normalized();
    double dist = normal.dot(v0 - v1);
    Eigen::Vector3d projected_point = v0 - dist * normal;
    // 2. Compute barycentric coordinates of projected point ([α₁, α₂, α₃])
    Eigen::Vector3d barycentric_coords;
    barycentric_coordinates(projected_point, v1, v2, v3, barycentric_coords);
    // 3. ∇C = [N, -α₁N, -α₂N, -α₃N]
    for (int i = 0; i < 4; i++) {
        double coef = i == 0 ? 1 : (-barycentric_coords(i - 1));
        grad_c.segment<3>(3 * i) = coef * normal;
    }
}

void compute_graphics_edge_edge_constraint_gradient(
    const Eigen::Vector3d& v0, // first edge's first vertex
    const Eigen::Vector3d& v1, // first edge's second vertex
    const Eigen::Vector3d& v2, // second edge's first vertex
    const Eigen::Vector3d& v3, // second edge's second vertex
    Eigen::Vector12d& grad_c)
{
    // 1. Find closest points between e0 and e1 (see: https://bit.ly/2qwI33t)
    Eigen::Vector3d dir0 = v1 - v0;
    Eigen::Vector3d dir1 = v3 - v2;
    Eigen::Vector3d dir2 = dir1.cross(dir0);
    // v0 + t0 * dir0 + t2 * dir2 = v2 + t1 * dir1
    // [dir0, -dir1, dir2] * [t0; t1; t2] = v2 - v0
    Eigen::Matrix3d A;
    A.col(0) = dir0;
    A.col(1) = -dir1;
    A.col(2) = dir2;
    Eigen::Vector3d params = A.lu().solve(v2 - v0);
    if (!(params.array().isFinite().all())) {
        spdlog::error("Unable to compute edge-edge graphics constraint gradient!");
        spdlog::debug("Linear solve for parameters along edges resulted in ({}); e₁×e₂=({})",
            toString(params.transpose()), toString(dir2.transpose()));
        grad_c.setZero();
        return;
    }
    assert((A * params - (v2 - v0)).squaredNorm() < 1e-12);
    // 3. ∇C = [-(1 - α₀)N, -α₀N, (1 - α₁)N, α₁N]
    params(0) = std::clamp(params(0), 0.0, 1.0);
    params(1) = std::clamp(params(1), 0.0, 1.0);
    Eigen::Vector3d normal = (v3 - v2).cross(v1 - v0).normalized();
    grad_c.segment<3>(0) = -(1 - params(0)) * normal;
    grad_c.segment<3>(3) = -params(0) * normal;
    grad_c.segment<3>(6) = (1 - params(1)) * normal;
    grad_c.segment<3>(9) = params(1) * normal;
}

///////////////////////////////////////////////////////////////////////////////
// Efficient and Accurate Collision Response for Elastically Deformable Models
// [Verschoor et al. 2019]

void compute_Verschoor_point_triangle_constraint(
    const Eigen::Vector3d& v0_t0, // point at start of the timestep
    const Eigen::Vector3d& v1_t0, // triangle point 0 at start of the timestep
    const Eigen::Vector3d& v2_t0, // triangle point 1 at start of the timestep
    const Eigen::Vector3d& v3_t0, // triangle point 2 at start of the timestep
    const Eigen::Vector3d& v0_t1, // point at end of the timestep
    const Eigen::Vector3d& v1_t1, // triangle point 0 at end of the timestep
    const Eigen::Vector3d& v2_t1, // triangle point 1 at end of the timestep
    const Eigen::Vector3d& v3_t1, // triangle point 2 at end of the timestep
    double toi, double& c)
{
    // 1. Find the point of contact using CCD
    if (toi < 0 || toi > 1) {
        c = 1e28;
        return;
    }
    const Eigen::Vector3d vc_toi = (v0_t1 - v0_t0) * toi + v0_t0;
    // 2. Compute barycentric coordinates of contact point ([α₁, α₂, α₃])
    Eigen::Vector3d barycentric_coords;
    const Eigen::Vector3d v1_toi = (v1_t1 - v1_t0) * toi + v1_t0;
    const Eigen::Vector3d v2_toi = (v2_t1 - v2_t0) * toi + v2_t0;
    const Eigen::Vector3d v3_toi = (v3_t1 - v3_t0) * toi + v3_t0;
    barycentric_coordinates(vc_toi, v1_toi, v2_toi, v3_toi, barycentric_coords);
    // 3. Compute contact point's position at the start of the iteration
    const Eigen::Vector3d vc_t1 = barycentric_coords(0) * v1_t1
        + barycentric_coords(1) * v2_t1 + barycentric_coords(2) * v3_t1;
    // 4. Compute the contact normal
    // Is the normal at time of impact?
    // Eigen::Vector3d normal = (v2_toi - v1_toi).cross(v3_toi - v1_toi).normalized();
    // or is the normal at start of iteration?
    Eigen::Vector3d normal = (v2_t1 - v1_t1).cross(v3_t1 - v1_t1).normalized();
    // 5. Compute the distance constraint
    c = normal.dot(v0_t1 - vc_t1); // distance to plane
}

void compute_Verschoor_edge_edge_constraint(
    const Eigen::Vector3d& v0_t0, // first edge's first vertex at t = 0
    const Eigen::Vector3d& v1_t0, // first edge's second vertex at t = 0
    const Eigen::Vector3d& v2_t0, // second edge's first vertex at t = 0
    const Eigen::Vector3d& v3_t0, // second edge's second vertex at t = 0
    const Eigen::Vector3d& v0_t1, // first edge's first vertex at t = 1
    const Eigen::Vector3d& v1_t1, // first edge's second vertex at t = 1
    const Eigen::Vector3d& v2_t1, // second edge's first vertex at t = 1
    const Eigen::Vector3d& v3_t1, // second edge's second vertex at t = 1
    double toi, double& c)
{
    // 1. Find the time of impact using CCD
    if (toi < 0 || toi > 1 || !std::isfinite(toi)) {
        c = 1e28;
        return;
    }
    const Eigen::Vector3d v0_toi = (v0_t1 - v0_t0) * toi + v0_t0;
    const Eigen::Vector3d v1_toi = (v1_t1 - v1_t0) * toi + v1_t0;
    const Eigen::Vector3d v2_toi = (v2_t1 - v2_t0) * toi + v2_t0;
    const Eigen::Vector3d v3_toi = (v3_t1 - v3_t0) * toi + v3_t0;
    // 2. Find points of contact between e0 and e1 (see: https://bit.ly/2qwI33t)
    Eigen::Vector3d dir0 = v1_toi - v0_toi;
    Eigen::Vector3d dir1 = v3_toi - v2_toi;
    Eigen::Vector3d dir2 = dir1.cross(dir0);
    // v0 + t0 * dir0 + t2 * dir2 = v2 + t1 * dir1
    // [dir0, -dir1, dir2] * [t0; t1; t2] = v2 - v0
    // TODO: Make this more efficient.
    Eigen::Matrix3d A;
    A.col(0) = dir0;
    A.col(1) = -dir1;
    A.col(2) = dir2;
    Eigen::Vector3d params = A.lu().solve(v2_toi - v0_toi);
    // Check the solution is valid
    if (!(params.array().isFinite().all())) {
        spdlog::error("Unable to compute edge-edge Verschoor constraint!");
        spdlog::debug("Linear solve for parameters along edges resulted in ({}); e₁×e₂=({})",
            toString(params.transpose()), toString(dir2.transpose()));
        c = 1e28;
        return;
    }
    assert((A * params - (v2_toi - v0_toi)).squaredNorm() < 1e-12);
    // 3. Project the points back onto the edges.
    Eigen::Vector3d vc0_t1 = (v1_t1 - v0_t1) * std::clamp(params(0), 0.0, 1.0) + v0_t1;
    Eigen::Vector3d vc1_t1 = (v3_t1 - v2_t1) * std::clamp(params(1), 0.0, 1.0) + v2_t1;
    // 4. Compute the contact normal
    // Is the normal at time of impact?
    // Eigen::Vector3d normal = (v3_toi - v2_toi).cross(v1_toi - v0_toi).normalized();
    // or is the normal at start of iteration?
    Eigen::Vector3d normal = (v3_t1 - v2_t1).cross(v1_t1 - v0_t1).normalized();
    // 5. Compute the distance from the points doted with the normal.
    c = normal.dot(vc1_t1 - vc0_t1);
}

void compute_Verschoor_point_triangle_constraint_gradient(
    const Eigen::Vector3d& v0_t0, // point at start of the timestep
    const Eigen::Vector3d& v1_t0, // triangle point 0 at start of the timestep
    const Eigen::Vector3d& v2_t0, // triangle point 1 at start of the timestep
    const Eigen::Vector3d& v3_t0, // triangle point 2 at start of the timestep
    const Eigen::Vector3d& v0_t1, // point at end of the timestep
    const Eigen::Vector3d& v1_t1, // triangle point 0 at end of the timestep
    const Eigen::Vector3d& v2_t1, // triangle point 1 at end of the timestep
    const Eigen::Vector3d& v3_t1, // triangle point 2 at end of the timestep
    double toi, Eigen::Vector12d& grad_c)
{
    // 1. Find the point of contact using CCD
    if (toi < 0 || toi > 1 || !std::isfinite(toi)) {
        grad_c.setZero();
        return;
    }
    const Eigen::Vector3d vc_toi = (v0_t1 - v0_t0) * toi + v0_t0;
    // 2. Compute barycentric coordinates of contact point ([α₁, α₂, α₃])
    Eigen::Vector3d barycentric_coords;
    const Eigen::Vector3d v1_toi = (v1_t1 - v1_t0) * toi + v1_t0;
    const Eigen::Vector3d v2_toi = (v2_t1 - v2_t0) * toi + v2_t0;
    const Eigen::Vector3d v3_toi = (v3_t1 - v3_t0) * toi + v3_t0;
    barycentric_coordinates(vc_toi, v1_toi, v2_toi, v3_toi, barycentric_coords);
    // 3. Compute the contact normal
    // Is the normal at time of impact?
    // Eigen::Vector3d normal = (v2_toi - v1_toi).cross(v3_toi - v1_toi).normalized();
    // or is the normal at start of iteration?
    Eigen::Vector3d normal = (v2_t1 - v1_t1).cross(v3_t1 - v1_t1).normalized();
    // 4. ∇C = [N, -α₁N, -α₂N, -α₃N]
    for (int i = 0; i < 4; i++) {
        double coef = i == 0 ? 1 : (-barycentric_coords(i - 1));
        grad_c.segment<3>(3 * i) = coef * normal;
    }
}

void compute_Verschoor_edge_edge_constraint_gradient(
    const Eigen::Vector3d& v0_t0, // first edge's first vertex at t = 0
    const Eigen::Vector3d& v1_t0, // first edge's second vertex at t = 0
    const Eigen::Vector3d& v2_t0, // second edge's first vertex at t = 0
    const Eigen::Vector3d& v3_t0, // second edge's second vertex at t = 0
    const Eigen::Vector3d& v0_t1, // first edge's first vertex at t = 1
    const Eigen::Vector3d& v1_t1, // first edge's second vertex at t = 1
    const Eigen::Vector3d& v2_t1, // second edge's first vertex at t = 1
    const Eigen::Vector3d& v3_t1, // second edge's second vertex at t = 1
    double toi, Eigen::Vector12d& grad_c)
{
    // 1. Find the time of impact using CCD
    if (toi < 0 || toi > 1 || !std::isfinite(toi)) {
        grad_c.setZero();
        return;
    }
    const Eigen::Vector3d v0_toi = (v0_t1 - v0_t0) * toi + v0_t0;
    const Eigen::Vector3d v1_toi = (v1_t1 - v1_t0) * toi + v1_t0;
    const Eigen::Vector3d v2_toi = (v2_t1 - v2_t0) * toi + v2_t0;
    const Eigen::Vector3d v3_toi = (v3_t1 - v3_t0) * toi + v3_t0;
    // 2. Find points of contact between e0 and e1 (see: https://bit.ly/2qwI33t)
    Eigen::Vector3d dir0 = v1_toi - v0_toi;
    Eigen::Vector3d dir1 = v3_toi - v2_toi;
    Eigen::Vector3d dir2 = dir1.cross(dir0);
    // v0 + t0 * dir0 + t2 * dir2 = v2 + t1 * dir1
    // [dir0, -dir1, dir2] * [t0; t1; t2] = v2 - v0
    // TODO: Make this more efficient.
    Eigen::Matrix3d A;
    A.col(0) = dir0;
    A.col(1) = -dir1;
    A.col(2) = dir2;
    Eigen::Vector3d params = A.lu().solve(v2_toi - v0_toi);
    // Check the solution is valid
    if (!(params.array().isFinite().all())) {
        spdlog::error("Unable to compute edge-edge Verschoor constraint gradient!");
        spdlog::debug("Linear solve for parameters along edges resulted in ({}); e₁×e₂=({})",
            toString(params.transpose()), toString(dir2.transpose()));
        grad_c.setZero();
        return;
    }
    assert((A * params - (v2_toi - v0_toi)).squaredNorm() < 1e-12);
    // 3. Compute the contact normal
    // Is the normal at time of impact?
    // Eigen::Vector3d normal = (v3_toi - v2_toi).cross(v1_toi - v0_toi).normalized();
    // or is the normal at start of iteration?
    Eigen::Vector3d normal = (v3_t1 - v2_t1).cross(v1_t1 - v0_t1).normalized();
    // 4. ∇C = [-(1 - α₀)N, -α₀N, (1 - α₁)N, α₁N]
    params(0) = std::clamp(params(0), 0.0, 1.0);
    params(1) = std::clamp(params(1), 0.0, 1.0);
    grad_c.segment<3>(0) = -(1 - params(0)) * normal;
    grad_c.segment<3>(3) = -params(0) * normal;
    grad_c.segment<3>(6) = (1 - params(1)) * normal;
    grad_c.segment<3>(9) = params(1) * normal;
}

////////////////////////////////////////////////////////////////////////////////
// Parallel contact-aware simulations of deformable particles in 3D Stokes flow
// [Lu et al. 2018]

void compute_STIV_point_triangle_constraint(
    const Eigen::Vector3d& v0_t0, // point at start of the timestep
    const Eigen::Vector3d& v1_t0, // triangle point 0 at start of the timestep
    const Eigen::Vector3d& v2_t0, // triangle point 1 at start of the timestep
    const Eigen::Vector3d& v3_t0, // triangle point 2 at start of the timestep
    const Eigen::Vector3d& v0_t1, // point at end of the timestep
    const Eigen::Vector3d& v1_t1, // triangle point 0 at end of the timestep
    const Eigen::Vector3d& v2_t1, // triangle point 1 at end of the timestep
    const Eigen::Vector3d& v3_t1, // triangle point 2 at end of the timestep
    double toi, double& c)
{
    if (toi < 0 || toi > 1 || !std::isfinite(toi)) {
        c = 1e28;
        return;
    }
    // V(t,X) = (1 - τ) * √(ϵ² + (X_k(t_{n+1}) - X_k(t_{n+1}))⋅n(τ))²) * |t|
    // 0. Avoid zero volume from parallel displacment and normal.
    const double EPS = 1e-3;
    // 1. Compute the displacment of the vertex
    const Eigen::Vector3d u0 = v0_t1 - v0_t0;
    // 2. Compute normal of the triangle at the time of impact
    const Eigen::Vector3d v1_toi = (v1_t1 - v1_t0) * toi + v1_t0;
    const Eigen::Vector3d v2_toi = (v2_t1 - v2_t0) * toi + v2_t0;
    const Eigen::Vector3d v3_toi = (v3_t1 - v3_t0) * toi + v3_t0;
    Eigen::Vector3d n_toi = (v2_toi - v1_toi).cross(v3_toi - v1_toi);
    // 4. Compute the area of triangle
    const double tri_area = n_toi.norm() / 2;
    n_toi.normalize();
    // 5. Compute the STIV
    c = -(1 - toi) * sqrt(EPS * EPS + pow(u0.dot(n_toi), 2)) * tri_area;
    spdlog::debug("STIV={:g}", c);
}

void compute_STIV_edge_edge_constraint(
    const Eigen::Vector3d& v0_t0, // first edge's first vertex at t = 0
    const Eigen::Vector3d& v1_t0, // first edge's second vertex at t = 0
    const Eigen::Vector3d& v2_t0, // second edge's first vertex at t = 0
    const Eigen::Vector3d& v3_t0, // second edge's second vertex at t = 0
    const Eigen::Vector3d& v0_t1, // first edge's first vertex at t = 1
    const Eigen::Vector3d& v1_t1, // first edge's second vertex at t = 1
    const Eigen::Vector3d& v2_t1, // second edge's first vertex at t = 1
    const Eigen::Vector3d& v3_t1, // second edge's second vertex at t = 1
    double toi, double& c)
{
    throw "not implemented";
}

void compute_STIV_point_triangle_constraint_gradient(
    const Eigen::Vector3d& v0_t0, // point at start of the timestep
    const Eigen::Vector3d& v1_t0, // triangle point 0 at start of the timestep
    const Eigen::Vector3d& v2_t0, // triangle point 1 at start of the timestep
    const Eigen::Vector3d& v3_t0, // triangle point 2 at start of the timestep
    const Eigen::Vector3d& v0_t1, // point at end of the timestep
    const Eigen::Vector3d& v1_t1, // triangle point 0 at end of the timestep
    const Eigen::Vector3d& v2_t1, // triangle point 1 at end of the timestep
    const Eigen::Vector3d& v3_t1, // triangle point 2 at end of the timestep
    double toi, Eigen::Vector12d& grad_c)
{
    // 1. Find the point of contact using CCD
    if (toi < 0 || toi > 1 || !std::isfinite(toi)) {
        grad_c.setZero();
        return;
    }
    const Eigen::Vector3d vc_toi = (v0_t1 - v0_t0) * toi + v0_t0;
    // 2. Compute barycentric coordinates of contact point ([α₁, α₂, α₃])
    Eigen::Vector3d barycentric_coords;
    const Eigen::Vector3d v1_toi = (v1_t1 - v1_t0) * toi + v1_t0;
    const Eigen::Vector3d v2_toi = (v2_t1 - v2_t0) * toi + v2_t0;
    const Eigen::Vector3d v3_toi = (v3_t1 - v3_t0) * toi + v3_t0;
    barycentric_coordinates(vc_toi, v1_toi, v2_toi, v3_toi, barycentric_coords);
    // 3. Compute the contact normal
    // Is the normal at time of impact?
    // Eigen::Vector3d normal = (v2_toi - v1_toi).cross(v3_toi - v1_toi).normalized();
    // or is the normal at start of iteration?
    Eigen::Vector3d normal = (v2_t1 - v1_t1).cross(v3_t1 - v1_t1).normalized();
    // 4. ∇C = [N, -α₁N, -α₂N, -α₃N]
    for (int i = 0; i < 4; i++) {
        double coef = i == 0 ? 1 : (-barycentric_coords(i - 1));
        grad_c.segment<3>(3 * i) = coef * normal;
    }
}

void compute_STIV_edge_edge_constraint_gradient(
    const Eigen::Vector3d& v0_t0, // first edge's first vertex at t = 0
    const Eigen::Vector3d& v1_t0, // first edge's second vertex at t = 0
    const Eigen::Vector3d& v2_t0, // second edge's first vertex at t = 0
    const Eigen::Vector3d& v3_t0, // second edge's second vertex at t = 0
    const Eigen::Vector3d& v0_t1, // first edge's first vertex at t = 1
    const Eigen::Vector3d& v1_t1, // first edge's second vertex at t = 1
    const Eigen::Vector3d& v2_t1, // second edge's first vertex at t = 1
    const Eigen::Vector3d& v3_t1, // second edge's second vertex at t = 1
    double toi, Eigen::Vector12d& grad_c)
{
    throw "not implemented";
}

} // namespace IPC
