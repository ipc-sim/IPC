/// @brief Eigen wrappers for the exact CCD methods.

#pragma once

#include <Eigen/Core>

namespace IPC {
namespace ExactCCD {

/// @brief Methods of geometrically exact continous collision detection.
enum Method {
    NONE, ///< @brief Do not use exact CCD.
    ROOT_PARITY, ///< @brief Root parity method of Brochu et al. [2012].
    BSC, ///< @brief Bernstein sign classification method of Tang et al. [2014].
    RATIONAL_ROOT_PARITY ///< @brief Teseo's reimplementation of Brochu et al. [2012] using rationals
};

/**
 * @brief Detect collisions between two edges as they move.
 *
 * Looks for collisions between edges (q0start, p0start) and (q1start, p1start)
 * as they move towards (q0end, p0end) and (q1end, p1end). Returns true if the
 * edges collide.
 *
 * @param[in]  q0start  Start position of the first edge's first vertex.
 * @param[in]  p0start  Start position of the first edge's second vertex.
 * @param[in]  q1start  Start position of the second edge's first vertex.
 * @param[in]  p1start  Start position of the second edge's second vertex.
 * @param[in]  q0end    End position of the first edge's first vertex.
 * @param[in]  p0end    End position of the first edge's second vertex.
 * @param[in]  q1end    End position of the second edge's first vertex.
 * @param[in]  p1end    End position of the second edge's second vertex.
 * @param[in]  method   Method of exact CCD.
 *
 * @returns True if the edges collide.
 */
bool edgeEdgeCCD(
    const Eigen::Vector3d& q0start,
    const Eigen::Vector3d& p0start,
    const Eigen::Vector3d& q1start,
    const Eigen::Vector3d& p1start,
    const Eigen::Vector3d& q0end,
    const Eigen::Vector3d& p0end,
    const Eigen::Vector3d& q1end,
    const Eigen::Vector3d& p1end,
    const Method method);

/**
 * @brief Detect collisions between a vertex and a triangular face.
 *
 * Looks for collisions between the vertex q0start and the face
 * (q1start, q2start, q3start) as they move towards q0end and
 * (q1end, q2end, q3end). Returns true if the vertex and face collide.
 *
 * @param[in]  q0start  Start position of the vertex.
 * @param[in]  q1start  Start position of the first vertex of the face.
 * @param[in]  q2start  Start position of the second vertex of the face.
 * @param[in]  q3start  Start position of the third vertex of the face.
 * @param[in]  q0end    End position of the vertex.
 * @param[in]  q1end    End position of the first vertex of the face.
 * @param[in]  q2end    End position of the second vertex of the face.
 * @param[in]  q3end    End position of the third vertex of the face.
 * @param[in]  method   Method of exact CCD.
 *
 * @returns  True if the vertex and face collide.
 */
bool vertexFaceCCD(
    const Eigen::Vector3d& q0start,
    const Eigen::Vector3d& q1start,
    const Eigen::Vector3d& q2start,
    const Eigen::Vector3d& q3start,
    const Eigen::Vector3d& q0end,
    const Eigen::Vector3d& q1end,
    const Eigen::Vector3d& q2end,
    const Eigen::Vector3d& q3end,
    const Method method);

}
} // namespace IPC::ExactCCD
