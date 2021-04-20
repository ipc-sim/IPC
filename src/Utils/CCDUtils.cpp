#include <CCDUtils.hpp>

#ifdef USE_FPRP_CCD
#include <doubleCCD/doubleccd.hpp>
#endif

#include <spdlog/spdlog.h>

namespace IPC {

Eigen::Vector3d invShift = Eigen::Vector3d::Zero();

template <int dim>
Eigen::Vector3d shiftVertices(
    Mesh<dim>& mesh,
    std::vector<CollisionObject<dim>*>& meshCollisionObjects)
{
    assert(dim == 3);

    // compute a bounding box for the entire scene
    Eigen::Vector3d world_min = mesh.V.colwise().minCoeff().transpose();
    Eigen::Vector3d world_max = mesh.V.colwise().maxCoeff().transpose();
    for (CollisionObject<dim>* meshCO : meshCollisionObjects) {
        world_min = world_min.cwiseMin(meshCO->V.colwise().minCoeff().transpose());
        world_max = world_max.cwiseMax(meshCO->V.colwise().maxCoeff().transpose());
    }

    double world_radius = (world_max - world_min).norm() / 2;
    Eigen::Vector3d world_center = (world_max + world_min) / 2;
    spdlog::info(
        "world_min=[{:g}, {:g}, {:g}] world_max=[{:g}, {:g}, {:g}] world_radius={:g}",
        world_min.x(), world_min.y(), world_min.z(),
        world_max.x(), world_max.y(), world_max.z(),
        world_radius);

    world_min = world_center - 10 * world_radius * Eigen::Vector3d(1, 1, 1).normalized();
    world_max = world_center + 10 * world_radius * Eigen::Vector3d(1, 1, 1).normalized();

    spdlog::info(
        "using a world bbox of [{:g}, {:g}, {:g}]×[{:g}, {:g}, {:g}]",
        world_min.x(), world_min.y(), world_min.z(),
        world_max.x(), world_max.y(), world_max.z());

    return shiftVertices(mesh, meshCollisionObjects, world_min, world_max);
}

template <int dim>
Eigen::Vector3d shiftVertices(
    Mesh<dim>& mesh,
    std::vector<CollisionObject<dim>*>& meshCollisionObjects,
    const Eigen::Vector3d& world_min,
    const Eigen::Vector3d& world_max)
{
    assert(dim == 3);

#ifdef USE_FPRP_CCD
    // Stack all vertices to shift by the same amount
    size_t num_vertices = mesh.V.rows();
    for (CollisionObject<dim>* meshCO : meshCollisionObjects) {
        num_vertices += meshCO->V.rows();
    }
    Eigen::MatrixX3d all_vertices(num_vertices, 3);
    size_t start_vi = 0;
    for (CollisionObject<dim>* meshCO : meshCollisionObjects) {
        all_vertices.middleRows(start_vi, meshCO->V.rows()) = meshCO->V;
        start_vi += meshCO->V.rows();
    }
    all_vertices.bottomRows(mesh.V.rows()) = mesh.V;

    Eigen::Vector3d invShift;
    double shift_err = doubleccd::get_whole_mesh_shifted(
        all_vertices, world_min, world_max, invShift);
    spdlog::info("FPRP shift_err={:.17g}", shift_err);
    spdlog::info(
        "shifting world by [{:g}, {:g}, {:g}]",
        -invShift.x(), -invShift.y(), -invShift.z());

    // Split all vertices
    start_vi = 0;
    for (CollisionObject<dim>* meshCO : meshCollisionObjects) {
        meshCO->V = all_vertices.middleRows(start_vi, meshCO->V.rows());
        start_vi += meshCO->V.rows();
    }
    mesh.V = all_vertices.bottomRows(mesh.V.rows());

    return invShift;
#else
    throw "FPRP disabled, unable to shift vertices";
#endif
}

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
    double eps_t)
{
    if (ccd::is_time_of_impact_computed(ccdMetod)) {
        spdlog::warn(
            "Using vertexFaceToIBisection with a method that directly computes "
            "the time of impact.");
    }

    bool colliding = ccd::vertexFaceCCD(
        v_t0, f0_t0, f1_t0, f2_t0,
        v_t1, f0_t1, f1_t1, f2_t1,
        ccdMetod);
    if (!colliding) {
        return false;
    }

    double t0 = 0, t1 = 1;
    while (abs(t1 - t0) > eps_t) {
        double t = (t0 + t1) / 2;
        Eigen::VectorXd v_t = (v_t1 - v_t0) * t + v_t0;
        Eigen::VectorXd f0_t = (f0_t1 - f0_t0) * t + f0_t0;
        Eigen::VectorXd f1_t = (f1_t1 - f1_t0) * t + f1_t0;
        Eigen::VectorXd f2_t = (f2_t1 - f2_t0) * t + f2_t0;

        colliding = ccd::vertexFaceCCD(
            v_t0, f0_t0, f1_t0, f2_t0,
            v_t, f0_t, f1_t, f2_t,
            ccdMetod);

        if (colliding) {
            t1 = t;
        }
        else {
            t0 = t;
        }
    }
    toi = t0; // conservative answer
    return true;
}

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
    double eps_t)
{
    // Assumes there is a collision
    if (ccd::is_time_of_impact_computed(ccdMetod)) {
        spdlog::warn(
            "Using edgeEdgeToIBisection with a method that directly computes "
            "the time of impact.");
    }

    bool colliding = ccd::edgeEdgeCCD(
        ea0_t0, ea1_t0, eb0_t0, eb1_t0,
        ea0_t1, ea1_t1, eb0_t1, eb1_t1,
        ccdMetod);
    if (!colliding) {
        return false;
    }

    double t0 = 0, t1 = 1;
    while (abs(t1 - t0) > eps_t) {
        double t = (t0 + t1) / 2;
        Eigen::VectorXd ea0_t = (ea0_t1 - ea0_t0) * t + ea0_t0;
        Eigen::VectorXd ea1_t = (ea1_t1 - ea1_t0) * t + ea1_t0;
        Eigen::VectorXd eb0_t = (eb0_t1 - eb0_t0) * t + eb0_t0;
        Eigen::VectorXd eb1_t = (eb1_t1 - eb1_t0) * t + eb1_t0;

        colliding = ccd::edgeEdgeCCD(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0,
            ea0_t, ea1_t, eb0_t, eb1_t,
            ccdMetod);

        if (colliding) {
            t1 = t;
        }
        else {
            t0 = t;
        }
    }
    toi = t0; // conservative answer
    return true;
}

///////////////////////////////////////////////////////////////////////////////

template Eigen::Vector3d shiftVertices(
    Mesh<3>& mesh,
    std::vector<CollisionObject<3>*>& meshCollisionObjects);
template Eigen::Vector3d shiftVertices(
    Mesh<3>& mesh,
    std::vector<CollisionObject<3>*>& meshCollisionObjects,
    const Eigen::Vector3d& world_min,
    const Eigen::Vector3d& world_max);

} // namespace IPC
