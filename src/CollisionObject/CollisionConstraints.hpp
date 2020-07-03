/// @brief Constraint functions for use with constraint based collisions.

#pragma once

#include <Eigen/Core>

namespace Eigen {
typedef Matrix<double, 12, 1> Vector12d;
}

namespace IPC {

/// @brief Methods for computing the constraints of the optimization problem.
enum CollisionConstraintType {
    VOLUME, ///< @brief Volume of a tetrahedron formed by the vertices.
    GRAPHICS, ///< @brief Distance constraint from Harmon et al. [2008].
    NONSMOOTH_NEWMARK, ///< @brief Variation of volume constraints proposed by Kane et al. [1999].
    GAP_FUNCTION, ///< @brief Common method from computational mechanics [Wriggers 1995]
    CMR, ///< @brief Constraint manifold refinement from Otaduy et al. [2009].
    VERSCHOOR, ///< @brief Variant of standard graphics approach by Verschoor and Jalba [2019].
    STIV ///< @brief Space-Time Interference Volume of Harmon et al. [2011] and Lu et al. [2018].
};

/**
 * @brief Compute the collision constraint value.
 *
 * Triangle vertices expected in clockwise order.
 *
 * @param[in]  v0_t0           First vertex at the start of time-step
 * @param[in]  v1_t0           Second vertex at the start of time-step
 * @param[in]  v2_t0           Third vertex at the start of time-step
 * @param[in]  v3_t0           Fourth vertex at the start of time-step
 * @param[in]  v0_t1           First vertex at the end of time-step
 * @param[in]  v1_t1           Second vertex at the end of time-step
 * @param[in]  v2_t1           Third vertex at the end of time-step
 * @param[in]  v3_t1           Fourth vertex at the end of time-step
 * @param[in]  constraintType  Type of collision constraint to compute.
 * @param[in]  is_edge_edge    Are the vertices of two edge:
 *                             (v0, v1) and (v2, v3)?
 * @param[in]  toi             Normalized time of impact between the two
 *                             primitives.
 * @param[out] c               Computed constraint value.
 */
void compute_collision_constraint(
    const Eigen::Vector3d& v0_t0, const Eigen::Vector3d& v1_t0,
    const Eigen::Vector3d& v2_t0, const Eigen::Vector3d& v3_t0,
    const Eigen::Vector3d& v0_t1, const Eigen::Vector3d& v1_t1,
    const Eigen::Vector3d& v2_t1, const Eigen::Vector3d& v3_t1,
    const CollisionConstraintType constraintType, bool is_edge_edge,
    double toi, double& c);

/**
 * @brief Compute the collision constraint gradient.
 *
 * Triangle vertices expected in clockwise order.
 *
 * @param[in]  v0_t0           First vertex at the start of time-step
 * @param[in]  v1_t0           Second vertex at the start of time-step
 * @param[in]  v2_t0           Third vertex at the start of time-step
 * @param[in]  v3_t0           Fourth vertex at the start of time-step
 * @param[in]  v0_t1           First vertex at the end of time-step
 * @param[in]  v1_t1           Second vertex at the end of time-step
 * @param[in]  v2_t1           Third vertex at the end of time-step
 * @param[in]  v3_t1           Fourth vertex at the end of time-step
 * @param[in]  constraintType  Type of collision constraint to compute.
 * @param[in]  is_edge_edge    Are the vertices of two edge:
 *                             (v0, v1) and (v2, v3)?
 * @param[in]  toi             Normalized time of impact between the two
 *                             primitives.
 * @param[out] grad_c          Computed gradient of constraint value.
 */
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
    const double toi, Eigen::Vector12d& grad_c);

///////////////////////////////////////////////////////////////////////////////
// Volume Constraint

/**
 * @brief Compute the collision constraint volume.
 *
 * Triangle vertices expected in clockwise order.
 */
void compute_collision_volume_constraint(
    const Eigen::Vector3d& v0,
    const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
    double& c);

/**
 * @brief Compute the collision constraint volume gradient.
 *
 * Triangle vertices expected in clockwise order.
 */
void compute_collision_volume_constraint_gradient(
    const Eigen::Vector3d& v0,
    const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
    Eigen::Vector12d& grad_c);

///////////////////////////////////////////////////////////////////////////////
// Graphics (Distance) Constraint

/**
* @brief Compute the barycentric coordinates of a point in triangle (a, b, c).
*
* Computes the barycentric coordinates of any point in the plane containing
* the triangle (a, b, c).
*
* @param[in]  p       Compute the barycentric coordinates of this point.
* @param[in]  a       First vertex in the triangle.
* @param[in]  b       Second vertex in the triangle.
* @param[in]  b       Third vertex in the triangle.
* @param[out] coords  Computed barycentric coordinates.
*/
void barycentric_coordinates(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c,
    Eigen::Vector3d& coords);

/**
 * @brief Compute the standard graphics point-triangle constraint value.
 *
 * Triangle vertices expected in counter-clockwise order.
 */
void compute_graphics_point_triangle_constraint(
    const Eigen::Vector3d& v0, // point
    const Eigen::Vector3d& v1, // triangle point 0
    const Eigen::Vector3d& v2, // triangle point 2
    const Eigen::Vector3d& v3, // triangle point 1
    double& c);

/**
 * @brief Compute the standard graphics edge-edge collision constraint value.
 */
void compute_graphics_edge_edge_constraint(
    const Eigen::Vector3d& v0, // first edge's first vertex
    const Eigen::Vector3d& v1, // first edge's second vertex
    const Eigen::Vector3d& v2, // second edge's first vertex
    const Eigen::Vector3d& v3, // second edge's second vertex
    double& c);

/**
 * @brief Compute the standard graphics point-triangle constraint gradient.
 *
 * Triangle vertices expected in counter-clockwise order.
 */
void compute_graphics_point_triangle_constraint_gradient(
    const Eigen::Vector3d& v0, // point
    const Eigen::Vector3d& v1, // triangle point 0
    const Eigen::Vector3d& v2, // triangle point 2
    const Eigen::Vector3d& v3, // triangle point 1
    Eigen::Vector12d& grad_c);

/**
 * @brief Compute the standard graphics edge-edge collision constraint gradient.
 */
void compute_graphics_edge_edge_constraint_gradient(
    const Eigen::Vector3d& v0, // first edge's first vertex
    const Eigen::Vector3d& v1, // first edge's second vertex
    const Eigen::Vector3d& v2, // second edge's first vertex
    const Eigen::Vector3d& v3, // second edge's second vertex
    Eigen::Vector12d& grad_c);

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
    double toi, double& c);

void compute_Verschoor_edge_edge_constraint(
    const Eigen::Vector3d& v0_t0, // first edge's first vertex at t = 0
    const Eigen::Vector3d& v1_t0, // first edge's second vertex at t = 0
    const Eigen::Vector3d& v2_t0, // second edge's first vertex at t = 0
    const Eigen::Vector3d& v3_t0, // second edge's second vertex at t = 0
    const Eigen::Vector3d& v0_t1, // first edge's first vertex at t = 1
    const Eigen::Vector3d& v1_t1, // first edge's second vertex at t = 1
    const Eigen::Vector3d& v2_t1, // second edge's first vertex at t = 1
    const Eigen::Vector3d& v3_t1, // second edge's second vertex at t = 1
    double toi, double& c);

void compute_Verschoor_point_triangle_constraint_gradient(
    const Eigen::Vector3d& v0_t0, // point at start of the timestep
    const Eigen::Vector3d& v1_t0, // triangle point 0 at start of the timestep
    const Eigen::Vector3d& v2_t0, // triangle point 1 at start of the timestep
    const Eigen::Vector3d& v3_t0, // triangle point 2 at start of the timestep
    const Eigen::Vector3d& v0_t1, // point at end of the timestep
    const Eigen::Vector3d& v1_t1, // triangle point 0 at end of the timestep
    const Eigen::Vector3d& v2_t1, // triangle point 1 at end of the timestep
    const Eigen::Vector3d& v3_t1, // triangle point 2 at end of the timestep
    double toi, Eigen::Vector12d& grad_c);

void compute_Verschoor_edge_edge_constraint_gradient(
    const Eigen::Vector3d& v0_t0, // first edge's first vertex at t = 0
    const Eigen::Vector3d& v1_t0, // first edge's second vertex at t = 0
    const Eigen::Vector3d& v2_t0, // second edge's first vertex at t = 0
    const Eigen::Vector3d& v3_t0, // second edge's second vertex at t = 0
    const Eigen::Vector3d& v0_t1, // first edge's first vertex at t = 1
    const Eigen::Vector3d& v1_t1, // first edge's second vertex at t = 1
    const Eigen::Vector3d& v2_t1, // second edge's first vertex at t = 1
    const Eigen::Vector3d& v3_t1, // second edge's second vertex at t = 1
    double toi, Eigen::Vector12d& grad_c);

///////////////////////////////////////////////////////////////////////////////
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
    double toi, double& c);

void compute_STIV_edge_edge_constraint(
    const Eigen::Vector3d& v0_t0, // first edge's first vertex at t = 0
    const Eigen::Vector3d& v1_t0, // first edge's second vertex at t = 0
    const Eigen::Vector3d& v2_t0, // second edge's first vertex at t = 0
    const Eigen::Vector3d& v3_t0, // second edge's second vertex at t = 0
    const Eigen::Vector3d& v0_t1, // first edge's first vertex at t = 1
    const Eigen::Vector3d& v1_t1, // first edge's second vertex at t = 1
    const Eigen::Vector3d& v2_t1, // second edge's first vertex at t = 1
    const Eigen::Vector3d& v3_t1, // second edge's second vertex at t = 1
    double toi, double& c);

void compute_STIV_point_triangle_constraint_gradient(
    const Eigen::Vector3d& v0_t0, // point at start of the timestep
    const Eigen::Vector3d& v1_t0, // triangle point 0 at start of the timestep
    const Eigen::Vector3d& v2_t0, // triangle point 1 at start of the timestep
    const Eigen::Vector3d& v3_t0, // triangle point 2 at start of the timestep
    const Eigen::Vector3d& v0_t1, // point at end of the timestep
    const Eigen::Vector3d& v1_t1, // triangle point 0 at end of the timestep
    const Eigen::Vector3d& v2_t1, // triangle point 1 at end of the timestep
    const Eigen::Vector3d& v3_t1, // triangle point 2 at end of the timestep
    double toi, Eigen::Vector12d& grad_c);

void compute_STIV_edge_edge_constraint_gradient(
    const Eigen::Vector3d& v0_t0, // first edge's first vertex at t = 0
    const Eigen::Vector3d& v1_t0, // first edge's second vertex at t = 0
    const Eigen::Vector3d& v2_t0, // second edge's first vertex at t = 0
    const Eigen::Vector3d& v3_t0, // second edge's second vertex at t = 0
    const Eigen::Vector3d& v0_t1, // first edge's first vertex at t = 1
    const Eigen::Vector3d& v1_t1, // first edge's second vertex at t = 1
    const Eigen::Vector3d& v2_t1, // second edge's first vertex at t = 1
    const Eigen::Vector3d& v3_t1, // second edge's second vertex at t = 1
    double toi, Eigen::Vector12d& grad_c);

} // namespace IPC
