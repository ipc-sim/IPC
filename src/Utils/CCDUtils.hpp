#pragma once

#include <vector>

#include <Eigen/Core>

#include <ccd.hpp>

#include <CollisionObject/CollisionObject.h>
#include <Mesh.hpp>

namespace IPC {

static constexpr double DEFAULT_CCD_TOLERANCE = 1e-6;
#ifdef IPC_WITH_TIGHT_INCLUSION
static constexpr int TIGHT_INCLUSION_MAX_ITER = 1e6;
static constexpr int TIGHT_INCLUSION_CCD_TYPE = 1;
static constexpr bool TIGHT_INCLUSION_NO_ZERO_TOI = true;
static constexpr double TIGHT_INCLUSION_DIST_P = 0.2;
static constexpr double TIGHT_INCLUSION_MIN_DIST = DEFAULT_CCD_TOLERANCE;
#endif

#ifdef IPC_WITH_TIGHT_INCLUSION
extern std::array<double, 3> tight_inclusion_vf_err, tight_inclusion_ee_err;
#endif

extern Eigen::Vector3d invShift;

#ifdef IPC_WITH_TIGHT_INCLUSION
template <int dim>
void computeTightInclusionError(
    Mesh<dim>& mesh,
    std::vector<std::shared_ptr<CollisionObject<dim>>>& meshCollisionObjects);

template <int dim>
void computeTightInclusionError(
    Mesh<dim>& mesh,
    std::vector<std::shared_ptr<CollisionObject<dim>>>& meshCollisionObjects,
    const Eigen::Vector3d& world_min,
    const Eigen::Vector3d& world_max);
#endif

#ifdef IPC_WITH_FPRP_CCD
template <int dim>
Eigen::Vector3d shiftWorld(
    Mesh<dim>& mesh,
    std::vector<std::shared_ptr<CollisionObject<dim>>>& meshCollisionObjects);

template <int dim>
Eigen::Vector3d shiftWorld(
    Mesh<dim>& mesh,
    std::vector<std::shared_ptr<CollisionObject<dim>>>& meshCollisionObjects,
    const Eigen::Vector3d& world_min,
    const Eigen::Vector3d& world_max);
#endif

bool vertexFaceToIBisection(
    const Eigen::Vector3d& v_t0,
    const Eigen::Vector3d& f0_t0,
    const Eigen::Vector3d& f1_t0,
    const Eigen::Vector3d& f2_t0,
    const Eigen::Vector3d& v_t1,
    const Eigen::Vector3d& f0_t1,
    const Eigen::Vector3d& f1_t1,
    const Eigen::Vector3d& f2_t1,
    const ccd::CCDMethod ccdMetod,
    double& toi,
    double eps_t = 1e-3);

bool edgeEdgeToIBisection(
    const Eigen::Vector3d& ea0_t0,
    const Eigen::Vector3d& ea1_t0,
    const Eigen::Vector3d& eb0_t0,
    const Eigen::Vector3d& eb1_t0,
    const Eigen::Vector3d& ea0_t1,
    const Eigen::Vector3d& ea1_t1,
    const Eigen::Vector3d& eb0_t1,
    const Eigen::Vector3d& eb1_t1,
    const ccd::CCDMethod ccdMetod,
    double& toi,
    double eps_t = 1e-3);

} // namespace IPC
