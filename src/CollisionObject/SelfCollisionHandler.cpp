//
//  SelfCollisionHandler.cpp
//  IPC
//
//  Created by Minchen Li on 6/03/19.
//

#include "SelfCollisionHandler.hpp"
#include "FrictionUtils.hpp"

#include <stdexcept>

#ifdef USE_TBB
#include <mutex>
#include <tbb/parallel_for.h>
#endif

// Etienne Vouga's CCD using a root finder in floating points
#include <CTCD.h>
#ifdef IPC_WITH_TIGHT_INCLUSION
#include <tight_inclusion/inclusion_ccd.hpp>
#endif
#include "CCDUtils.hpp"

#include "get_feasible_steps.hpp"
#include "IglUtils.hpp"
#include "BarrierFunctions.hpp"
#include "Timer.hpp"

#include <spdlog/spdlog.h>

extern std::ofstream logFile;
extern Timer timer_temp3, timer_mt;

namespace IPC {

template <int dim>
void SelfCollisionHandler<dim>::evaluateConstraint(const Mesh<dim>& mesh,
    const MMCVID& MMCVIDI, double& val, double coef)
{
    if (MMCVIDI[0] >= 0) {
        // edge-edge
        d_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), val);
    }
    else {
        // point-triangle and degenerate edge-edge
        assert(MMCVIDI[1] >= 0);
        if (MMCVIDI[2] < 0) {
            // PP
            d_PP(mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(MMCVIDI[1]), val);
        }
        else if (MMCVIDI[3] < 0) {
            // PE
            d_PE(mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), val);
        }
        else {
            // PT
            d_PT(mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), val);
        }
    }
}

template <int dim>
void SelfCollisionHandler<dim>::evaluateConstraints(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& activeSet,
    Eigen::VectorXd& val, double coef)
{
    const int constraintStartInd = val.size();
    val.conservativeResize(constraintStartInd + activeSet.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)activeSet.size(), 1, [&](int cI)
#else
    for (int cI = 0; cI < activeSet.size(); ++cI)
#endif
        {
            evaluateConstraint(mesh, activeSet[cI], val[constraintStartInd + cI], coef);
        }
#ifdef USE_TBB
    );
#endif
}

template <int dim>
void SelfCollisionHandler<dim>::leftMultiplyConstraintJacobianT(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& activeSet,
    const Eigen::VectorXd& input,
    Eigen::VectorXd& output_incremental,
    double coef)
{
    // TODO: parallelize
    if constexpr (dim == 3) {
        int constraintI = 0;
        for (const auto& MMCVIDI : activeSet) {
            if (MMCVIDI[0] >= 0) {
                // edge-edge
                Eigen::Matrix<double, dim * 4, 1> g;
                g_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), g);
                g *= coef * input[constraintI];

                output_incremental.segment<dim>(MMCVIDI[0] * dim) += g.template segment<dim>(0);
                output_incremental.segment<dim>(MMCVIDI[1] * dim) += g.template segment<dim>(dim);
                output_incremental.segment<dim>(MMCVIDI[2] * dim) += g.template segment<dim>(dim * 2);
                output_incremental.segment<dim>(MMCVIDI[3] * dim) += g.template segment<dim>(dim * 3);
            }
            else {
                // point-triangle and degenerate edge-edge
                assert(MMCVIDI[1] >= 0);

                int v0I = -MMCVIDI[0] - 1;
                if (MMCVIDI[2] < 0) {
                    // PP
                    Eigen::Matrix<double, dim * 2, 1> g;
                    g_PP(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), g);
                    g *= coef * -MMCVIDI[3] * input[constraintI];

                    output_incremental.segment<dim>(v0I * dim) += g.template segment<dim>(0);
                    output_incremental.segment<dim>(MMCVIDI[1] * dim) += g.template segment<dim>(dim);
                }
                else if (MMCVIDI[3] < 0) {
                    // PE
                    Eigen::Matrix<double, dim * 3, 1> g;
                    g_PE(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), g);
                    g *= coef * -MMCVIDI[3] * input[constraintI];

                    output_incremental.segment<dim>(v0I * dim) += g.template segment<dim>(0);
                    output_incremental.segment<dim>(MMCVIDI[1] * dim) += g.template segment<dim>(dim);
                    output_incremental.segment<dim>(MMCVIDI[2] * dim) += g.template segment<dim>(dim * 2);
                }
                else {
                    // PT
                    Eigen::Matrix<double, dim * 4, 1> g;
                    g_PT(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), g);
                    g *= coef * input[constraintI];

                    output_incremental.segment<dim>(v0I * dim) += g.template segment<dim>(0);
                    output_incremental.segment<dim>(MMCVIDI[1] * dim) += g.template segment<dim>(dim);
                    output_incremental.segment<dim>(MMCVIDI[2] * dim) += g.template segment<dim>(dim * 2);
                    output_incremental.segment<dim>(MMCVIDI[3] * dim) += g.template segment<dim>(dim * 3);
                }
            }

            ++constraintI;
        }
    }
    else {
        // TODO: 2D collision
    }
}

// Evaluate the constraint for the SQP method.
template <int dim>
void SelfCollisionHandler<dim>::evaluateConstraintQP(
    const Mesh<dim>& mesh, const MMCVID& MMCVIDI,
    const CollisionConstraintType constraintType, double toi,
    double& val, double coef)
{
    if constexpr (dim == 3) {
        if (MMCVIDI[0] >= 0) {
            // edge-edge ++++
            compute_collision_constraint(
                // first edge at t = 0
                mesh.V_prev.row(MMCVIDI[0]), mesh.V_prev.row(MMCVIDI[1]),
                // second edge at t = 0
                mesh.V_prev.row(MMCVIDI[2]), mesh.V_prev.row(MMCVIDI[3]),
                // first edge at t = 1
                mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]),
                // second edge at t = 1
                mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]),
                constraintType, /*is_edge_edge=*/true, toi,
                val);
            val *= coef;
        }
        else {
            // point-triangle -+++
            // no triangle-point self collisions
            assert(MMCVIDI[1] >= 0 && MMCVIDI[2] >= 0 && MMCVIDI[3] >= 0);
            compute_collision_constraint(
                // point at t = 1
                mesh.V_prev.row(-MMCVIDI[0] - 1),
                // triangle at t = 1
                mesh.V_prev.row(MMCVIDI[1]),
                mesh.V_prev.row(MMCVIDI[2]),
                mesh.V_prev.row(MMCVIDI[3]),
                // point at t = 1
                mesh.V.row(-MMCVIDI[0] - 1),
                // triangle at t = 1
                mesh.V.row(MMCVIDI[1]),
                mesh.V.row(MMCVIDI[2]),
                mesh.V.row(MMCVIDI[3]),
                constraintType, /*is_edge_edge=*/false, toi,
                val);
            val *= coef;
        }
        if (!std::isfinite(val)) {
            spdlog::error("Constraint value for MMCVID {:s} is non-finite ({:g})", MMCVIDI.str(), val);
            spdlog::debug("toi={:g}", toi);
            throw "constraint value is non-finite";
        }
    }
    else {
        // TODO: 2D collision
    }
}

// Evaluate the constraints for the SQP method.
template <int dim>
void SelfCollisionHandler<dim>::evaluateConstraintsQP(
    const Mesh<dim>& mesh, const std::vector<MMCVID>& activeSet,
    const CollisionConstraintType constraintType,
    const std::unordered_map<MMCVID, double, MMCVIDHash>& mmcvid_to_toi,
    Eigen::VectorXd& val, double coef)
{
    int constraintStartInd = val.size();
    val.conservativeResize(constraintStartInd + activeSet.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)activeSet.size(), 1, [&](int cI)
#else
    for (int cI = 0; cI < activeSet.size(); ++cI)
#endif
        {
            double toi;
            try {
                toi = mmcvid_to_toi.at(activeSet[cI]);
            }
            catch (std::out_of_range err) {
                spdlog::error("where=\"SelfCollisionHandler:evaluateConstraintsQP\" msg=\"Cannot compute the constraint value of MMCVID without toi!\"");
                spdlog::debug("MMCVIDI={:s}", activeSet[cI].str());
                spdlog::debug("mmcvid_to_toi:");
                for (auto const& [mmcvid, toi] : mmcvid_to_toi) {
                    spdlog::debug("{:s}: {:g}", mmcvid.str(), toi);
                }
                throw "missing toi for constraint";
            }

            SelfCollisionHandler<dim>::evaluateConstraintQP(
                mesh, activeSet[cI], constraintType, toi,
                val[constraintStartInd + cI], coef);
        }
#ifdef USE_TBB
    );
#endif
}

// Premultiply the input vector with the Jacobian of the constraints, and
// store the results in output_incremental.
template <int dim>
void SelfCollisionHandler<dim>::leftMultiplyConstraintJacobianTQP(
    const Mesh<dim>& mesh,
    const std::vector<MMCVID>& activeSet,
    const Eigen::VectorXd& input,
    const CollisionConstraintType constraintType,
    const std::unordered_map<MMCVID, double, MMCVIDHash>& mmcvid_to_toi,
    Eigen::VectorXd& output_incremental,
    double coef)
{
    // TODO: parallelize
    if constexpr (dim == 3) {
        int constraintI = 0;
        for (const auto& MMCVIDI : activeSet) {
            Eigen::Matrix<double, dim * 4, 1> grad;

            double toi;
            try {
                toi = mmcvid_to_toi.at(MMCVIDI);
            }
            catch (std::out_of_range err) {
                spdlog::error("where=\"SelfCollisionHandler:evaluateConstraintsQP\" msg=\"Cannot compute the constraint value of MMCVID without toi!\"");
                spdlog::debug("MMCVIDI={:s}", MMCVIDI.str());
                spdlog::debug("mmcvid_to_toi:");
                for (auto const& [mmcvid, toi] : mmcvid_to_toi) {
                    spdlog::debug("{:s}: {:g}", mmcvid.str(), toi);
                }
                throw "missing toi for constraint";
            }

            if (MMCVIDI[0] >= 0) {
                // edge-edge ++++
                compute_collision_constraint_gradient(
                    // first edge at t = 0
                    mesh.V_prev.row(MMCVIDI[0]), mesh.V_prev.row(MMCVIDI[1]),
                    // second edge at t = 0
                    mesh.V_prev.row(MMCVIDI[2]), mesh.V_prev.row(MMCVIDI[3]),
                    // first edge at t = 1
                    mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]),
                    // second edge at t = 1
                    mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]),
                    constraintType, /*is_edge_edge=*/true, toi,
                    grad);
                grad *= coef * input[constraintI];
                for (int i = 0; i < 4; i++) {
                    output_incremental.segment<dim>(MMCVIDI[i] * dim) += grad.template segment<dim>(i * dim);
                }
            }
            else {
                // point-triangle -+++
                // no triangle-point self collisions
                assert(MMCVIDI[1] >= 0 && MMCVIDI[2] >= 0 && MMCVIDI[3] >= 0);

                compute_collision_constraint_gradient(
                    // point at t = 0
                    mesh.V_prev.row(-MMCVIDI[0] - 1),
                    // triangle at t = 0
                    mesh.V_prev.row(MMCVIDI[1]),
                    mesh.V_prev.row(MMCVIDI[2]),
                    mesh.V_prev.row(MMCVIDI[3]),
                    // point at t = 1
                    mesh.V.row(-MMCVIDI[0] - 1),
                    // triangle at t = 1
                    mesh.V.row(MMCVIDI[1]),
                    mesh.V.row(MMCVIDI[2]),
                    mesh.V.row(MMCVIDI[3]),
                    constraintType, /*is_edge_edge=*/false, toi,
                    grad);
                grad *= coef * input[constraintI];
                output_incremental.segment<dim>((-MMCVIDI[0] - 1) * dim) += grad.template segment<dim>(0);
                output_incremental.segment<dim>(MMCVIDI[1] * dim) += grad.template segment<dim>(dim);
                output_incremental.segment<dim>(MMCVIDI[2] * dim) += grad.template segment<dim>(2 * dim);
                output_incremental.segment<dim>(MMCVIDI[3] * dim) += grad.template segment<dim>(3 * dim);
            }

            ++constraintI;
        }
    }
    else {
        // TODO: 2D collision
    }
}

template <int dim>
void SelfCollisionHandler<dim>::augmentConnectivity(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& activeSet,
    std::vector<std::set<int>>& vNeighbor)
{
    // TODO: parallelize?

    for (const auto& MMCVIDI : activeSet) {
        if (MMCVIDI[0] >= 0) {
            // edge-edge
            vNeighbor[MMCVIDI[0]].insert(MMCVIDI[2]);
            vNeighbor[MMCVIDI[2]].insert(MMCVIDI[0]);
            vNeighbor[MMCVIDI[0]].insert(MMCVIDI[3]);
            vNeighbor[MMCVIDI[3]].insert(MMCVIDI[0]);
            vNeighbor[MMCVIDI[1]].insert(MMCVIDI[2]);
            vNeighbor[MMCVIDI[2]].insert(MMCVIDI[1]);
            vNeighbor[MMCVIDI[1]].insert(MMCVIDI[3]);
            vNeighbor[MMCVIDI[3]].insert(MMCVIDI[1]);
        }
        else {
            // point-triangle and degenerate edge-edge
            assert(MMCVIDI[1] >= 0);

            int v0I = -MMCVIDI[0] - 1;
            if (MMCVIDI[2] < 0) {
                // PP
                vNeighbor[v0I].insert(MMCVIDI[1]);
                vNeighbor[MMCVIDI[1]].insert(v0I);
            }
            else if (MMCVIDI[3] < 0) {
                // PE
                vNeighbor[v0I].insert(MMCVIDI[1]);
                vNeighbor[MMCVIDI[1]].insert(v0I);
                vNeighbor[v0I].insert(MMCVIDI[2]);
                vNeighbor[MMCVIDI[2]].insert(v0I);
            }
            else {
                // PT
                vNeighbor[v0I].insert(MMCVIDI[1]);
                vNeighbor[MMCVIDI[1]].insert(v0I);
                vNeighbor[v0I].insert(MMCVIDI[2]);
                vNeighbor[MMCVIDI[2]].insert(v0I);
                vNeighbor[v0I].insert(MMCVIDI[3]);
                vNeighbor[MMCVIDI[3]].insert(v0I);
            }
        }
    }
}

template <int dim>
void SelfCollisionHandler<dim>::augmentConnectivity(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& paraEEMMCVIDSet,
    const std::vector<std::pair<int, int>>& paraEEeIeJSet,
    std::vector<std::set<int>>& vNeighbor)
{
    // TODO: parallelize?
    for (const auto& eIeJ : paraEEeIeJSet) {
        if (eIeJ.first < 0 || eIeJ.second < 0) {
            continue;
        }
        // PE or PP
        const std::pair<int, int>& eI = mesh.SFEdges[eIeJ.first];
        const std::pair<int, int>& eJ = mesh.SFEdges[eIeJ.second];
        vNeighbor[eI.first].insert(eJ.first);
        vNeighbor[eJ.first].insert(eI.first);
        vNeighbor[eI.first].insert(eJ.second);
        vNeighbor[eJ.second].insert(eI.first);
        vNeighbor[eI.second].insert(eJ.first);
        vNeighbor[eJ.first].insert(eI.second);
        vNeighbor[eI.second].insert(eJ.second);
        vNeighbor[eJ.second].insert(eI.second);
    }
    for (const auto& MMCVIDI : paraEEMMCVIDSet) {
        if (MMCVIDI[0] < 0) {
            continue;
        }
        // EE
        vNeighbor[MMCVIDI[0]].insert(MMCVIDI[2]);
        vNeighbor[MMCVIDI[2]].insert(MMCVIDI[0]);
        vNeighbor[MMCVIDI[0]].insert(MMCVIDI[3]);
        vNeighbor[MMCVIDI[3]].insert(MMCVIDI[0]);
        vNeighbor[MMCVIDI[1]].insert(MMCVIDI[2]);
        vNeighbor[MMCVIDI[2]].insert(MMCVIDI[1]);
        vNeighbor[MMCVIDI[1]].insert(MMCVIDI[3]);
        vNeighbor[MMCVIDI[3]].insert(MMCVIDI[1]);
    }
}

template <int dim>
void SelfCollisionHandler<dim>::augmentIPHessian(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& activeSet,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* mtr_incremental,
    double dHat, double coef, bool projectDBC)
{
    if constexpr (dim == 3) {
        std::vector<Eigen::Matrix<double, 4 * dim, 4 * dim>> IPHessian(activeSet.size());
        std::vector<Eigen::Matrix<int, 4, 1>> rowIStart(activeSet.size());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)activeSet.size(), 1, [&](int cI)
#else
        for (int cI = 0; cI < activeSet.size(); ++cI)
#endif
            {
                const auto& MMCVIDI = activeSet[cI];
                if (MMCVIDI[0] >= 0) {
                    // edge-edge
                    double d;
                    d_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), d);
                    Eigen::Matrix<double, dim * 4, 1> g;
                    g_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), g);
                    Eigen::Matrix<double, dim * 4, dim * 4> H;
                    H_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), H);

                    double g_b, H_b;
                    compute_g_b(d, dHat, g_b);
                    compute_H_b(d, dHat, H_b);

                    IPHessian[cI] = ((coef * H_b) * g) * g.transpose() + (coef * g_b) * H;
                    IglUtils::makePD(IPHessian[cI]);

                    rowIStart[cI][0] = mesh.isProjectDBCVertex(MMCVIDI[0], projectDBC) ? -1 : (MMCVIDI[0] * dim);
                    rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                    rowIStart[cI][2] = mesh.isProjectDBCVertex(MMCVIDI[2], projectDBC) ? -1 : (MMCVIDI[2] * dim);
                    rowIStart[cI][3] = mesh.isProjectDBCVertex(MMCVIDI[3], projectDBC) ? -1 : (MMCVIDI[3] * dim);
                }
                else {
                    // point-triangle and degenerate edge-edge
                    assert(MMCVIDI[1] >= 0);

                    int v0I = -MMCVIDI[0] - 1;
                    if (MMCVIDI[2] < 0) {
                        // PP
                        double d;
                        d_PP(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), d);
                        Eigen::Matrix<double, dim * 2, 1> g;
                        g_PP(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), g);
                        Eigen::Matrix<double, dim * 2, dim * 2> H;
                        H_PP(H);

                        double g_b, H_b;
                        compute_g_b(d, dHat, g_b);
                        compute_H_b(d, dHat, H_b);

                        double coef_dup = coef * -MMCVIDI[3];
                        Eigen::Matrix<double, dim * 2, dim* 2> HessianBlock = ((coef_dup * H_b) * g) * g.transpose() + (coef_dup * g_b) * H;
                        IglUtils::makePD(HessianBlock);
                        IPHessian[cI].template block<dim * 2, dim * 2>(0, 0) = HessianBlock;

                        rowIStart[cI][0] = mesh.isProjectDBCVertex(v0I, projectDBC) ? -1 : (v0I * dim);
                        rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                        rowIStart[cI][2] = -1;
                        rowIStart[cI][3] = -1;
                    }
                    else if (MMCVIDI[3] < 0) {
                        // PE
                        double d;
                        d_PE(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), d);
                        Eigen::Matrix<double, dim * 3, 1> g;
                        g_PE(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), g);
                        Eigen::Matrix<double, dim * 3, dim * 3> H;
                        H_PE(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), H);

                        double g_b, H_b;
                        compute_g_b(d, dHat, g_b);
                        compute_H_b(d, dHat, H_b);

                        double coef_dup = coef * -MMCVIDI[3];
                        Eigen::Matrix<double, dim * 3, dim* 3> HessianBlock = ((coef_dup * H_b) * g) * g.transpose() + (coef_dup * g_b) * H;
                        IglUtils::makePD(HessianBlock);
                        IPHessian[cI].block(0, 0, dim * 3, dim * 3) = HessianBlock;

                        rowIStart[cI][0] = mesh.isProjectDBCVertex(v0I, projectDBC) ? -1 : (v0I * dim);
                        rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                        rowIStart[cI][2] = mesh.isProjectDBCVertex(MMCVIDI[2], projectDBC) ? -1 : (MMCVIDI[2] * dim);
                        rowIStart[cI][3] = -1;
                    }
                    else {
                        // PT
                        double d;
                        d_PT(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), d);
                        Eigen::Matrix<double, dim * 4, 1> g;
                        g_PT(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), g);
                        Eigen::Matrix<double, dim * 4, dim * 4> H;
                        H_PT(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), H);

                        double g_b, H_b;
                        compute_g_b(d, dHat, g_b);
                        compute_H_b(d, dHat, H_b);

                        IPHessian[cI] = ((coef * H_b) * g) * g.transpose() + (coef * g_b) * H;
                        IglUtils::makePD(IPHessian[cI]);

                        rowIStart[cI][0] = mesh.isProjectDBCVertex(v0I, projectDBC) ? -1 : (v0I * dim);
                        rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                        rowIStart[cI][2] = mesh.isProjectDBCVertex(MMCVIDI[2], projectDBC) ? -1 : (MMCVIDI[2] * dim);
                        rowIStart[cI][3] = mesh.isProjectDBCVertex(MMCVIDI[3], projectDBC) ? -1 : (MMCVIDI[3] * dim);
                    }
                }
            }
#ifdef USE_TBB
        );
#endif

        // TODO: parallelize
        for (int cI = 0; cI < activeSet.size(); ++cI) {
            for (int i = 0; i < rowIStart[cI].size(); ++i) {
                int rowIStartI = rowIStart[cI][i];
                if (rowIStartI >= 0) {
                    for (int j = 0; j < rowIStart[cI].size(); ++j) {
                        int colIStartI = rowIStart[cI][j];
                        if (colIStartI >= 0) {
                            mtr_incremental->addCoeff(rowIStartI, colIStartI, IPHessian[cI](i * dim, j * dim));
                            mtr_incremental->addCoeff(rowIStartI, colIStartI + 1, IPHessian[cI](i * dim, j * dim + 1));
                            mtr_incremental->addCoeff(rowIStartI + 1, colIStartI, IPHessian[cI](i * dim + 1, j * dim));
                            mtr_incremental->addCoeff(rowIStartI + 1, colIStartI + 1, IPHessian[cI](i * dim + 1, j * dim + 1));
                            if constexpr (dim == 3) {
                                mtr_incremental->addCoeff(rowIStartI, colIStartI + 2, IPHessian[cI](i * dim, j * dim + 2));
                                mtr_incremental->addCoeff(rowIStartI + 1, colIStartI + 2, IPHessian[cI](i * dim + 1, j * dim + 2));

                                mtr_incremental->addCoeff(rowIStartI + 2, colIStartI, IPHessian[cI](i * dim + 2, j * dim));
                                mtr_incremental->addCoeff(rowIStartI + 2, colIStartI + 1, IPHessian[cI](i * dim + 2, j * dim + 1));
                                mtr_incremental->addCoeff(rowIStartI + 2, colIStartI + 2, IPHessian[cI](i * dim + 2, j * dim + 2));
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        // TODO: 2D collision
    }
}

template <int dim>
void SelfCollisionHandler<dim>::largestFeasibleStepSize(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    const Eigen::VectorXd& searchDir,
    double slackness,
    const std::vector<std::pair<int, int>>& constraintSet,
    std::vector<std::pair<int, int>>& candidates,
    double& stepSize)
{
    const double CCDDistRatio = 1.0 - slackness;

#if (CFL_FOR_CCD != 0)
    if (constraintSet.size()) {
        Eigen::VectorXd largestAlphasAS(constraintSet.size());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)constraintSet.size(), 1, [&](int cI)
#else
        for (int cI = 0; cI < constraintSet.size(); ++cI)
#endif
            {
                MMCVID MMCVIDI;
                if (constraintSet[cI].first < 0) {
                    // PT
                    MMCVIDI[0] = -mesh.SVI[-constraintSet[cI].first - 1] - 1;
                    MMCVIDI[1] = mesh.SF(constraintSet[cI].second, 0);
                    MMCVIDI[2] = mesh.SF(constraintSet[cI].second, 1);
                    MMCVIDI[3] = mesh.SF(constraintSet[cI].second, 2);
                }
                else {
                    // EE
                    MMCVIDI[0] = mesh.SFEdges[constraintSet[cI].first].first;
                    MMCVIDI[1] = mesh.SFEdges[constraintSet[cI].first].second;
                    MMCVIDI[2] = mesh.SFEdges[constraintSet[cI].second].first;
                    MMCVIDI[3] = mesh.SFEdges[constraintSet[cI].second].second;
                }

                if (MMCVIDI[0] >= 0) {
                    // edge-edge
                    double d_sqrt;
                    computeEdgeEdgeD(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]),
                        mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), d_sqrt);
                    d_sqrt = std::sqrt(d_sqrt);

                    largestAlphasAS[cI] = 1.0;
                    if (CTCD::edgeEdgeCTCD(mesh.V.row(MMCVIDI[0]).transpose(),
                            mesh.V.row(MMCVIDI[1]).transpose(),
                            mesh.V.row(MMCVIDI[2]).transpose(),
                            mesh.V.row(MMCVIDI[3]).transpose(),
                            mesh.V.row(MMCVIDI[0]).transpose() + searchDir.segment<dim>(MMCVIDI[0] * dim),
                            mesh.V.row(MMCVIDI[1]).transpose() + searchDir.segment<dim>(MMCVIDI[1] * dim),
                            mesh.V.row(MMCVIDI[2]).transpose() + searchDir.segment<dim>(MMCVIDI[2] * dim),
                            mesh.V.row(MMCVIDI[3]).transpose() + searchDir.segment<dim>(MMCVIDI[3] * dim),
                            CCDDistRatio * d_sqrt,
                            largestAlphasAS[cI])) {
                        if (largestAlphasAS[cI] < 1.0e-6) {
                            if (CTCD::edgeEdgeCTCD(mesh.V.row(MMCVIDI[0]).transpose(),
                                    mesh.V.row(MMCVIDI[1]).transpose(),
                                    mesh.V.row(MMCVIDI[2]).transpose(),
                                    mesh.V.row(MMCVIDI[3]).transpose(),
                                    mesh.V.row(MMCVIDI[0]).transpose() + searchDir.segment<dim>(MMCVIDI[0] * dim),
                                    mesh.V.row(MMCVIDI[1]).transpose() + searchDir.segment<dim>(MMCVIDI[1] * dim),
                                    mesh.V.row(MMCVIDI[2]).transpose() + searchDir.segment<dim>(MMCVIDI[2] * dim),
                                    mesh.V.row(MMCVIDI[3]).transpose() + searchDir.segment<dim>(MMCVIDI[3] * dim),
                                    0.0, largestAlphasAS[cI])) {
                                if (largestAlphasAS[cI] == 0.0) {
                                    // numerically parallel edge-edge causes CCD code to fail
                                    largestAlphasAS[cI] = 1.0;
                                }
                                largestAlphasAS[cI] *= slackness;
                            }
                            else {
                                largestAlphasAS[cI] = 1.0;
                            }
                        }
                    }
                }
                else {
                    // point-triangle
                    int vI = -MMCVIDI[0] - 1;
                    assert(MMCVIDI[1] >= 0);

                    double d_sqrt;
                    computePointTriD(mesh.V.row(vI), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), d_sqrt);
                    d_sqrt = std::sqrt(d_sqrt);

                    largestAlphasAS[cI] = 1.0;
                    if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                            mesh.V.row(MMCVIDI[1]).transpose(),
                            mesh.V.row(MMCVIDI[2]).transpose(),
                            mesh.V.row(MMCVIDI[3]).transpose(),
                            mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                            mesh.V.row(MMCVIDI[1]).transpose() + searchDir.segment<dim>(MMCVIDI[1] * dim),
                            mesh.V.row(MMCVIDI[2]).transpose() + searchDir.segment<dim>(MMCVIDI[2] * dim),
                            mesh.V.row(MMCVIDI[3]).transpose() + searchDir.segment<dim>(MMCVIDI[3] * dim),
                            CCDDistRatio * d_sqrt,
                            largestAlphasAS[cI])) {
                        if (largestAlphasAS[cI] < 1.0e-6) {
                            if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                                    mesh.V.row(MMCVIDI[1]).transpose(),
                                    mesh.V.row(MMCVIDI[2]).transpose(),
                                    mesh.V.row(MMCVIDI[3]).transpose(),
                                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                                    mesh.V.row(MMCVIDI[1]).transpose() + searchDir.segment<dim>(MMCVIDI[1] * dim),
                                    mesh.V.row(MMCVIDI[2]).transpose() + searchDir.segment<dim>(MMCVIDI[2] * dim),
                                    mesh.V.row(MMCVIDI[3]).transpose() + searchDir.segment<dim>(MMCVIDI[3] * dim),
                                    0.0, largestAlphasAS[cI])) {
                                largestAlphasAS[cI] *= slackness;
                            }
                            else {
                                largestAlphasAS[cI] = 1.0;
                            }
                        }
                    }
                }
            }
#ifdef USE_TBB
        );
#endif
        stepSize = std::min(stepSize, largestAlphasAS.minCoeff());
    }
#else
    largestFeasibleStepSize_CCD(mesh, sh, searchDir, slackness, candidates, stepSize);
#endif
}

#ifdef IPC_WITH_TIGHT_INCLUSION
template <int dim>
void SelfCollisionHandler<dim>::largestFeasibleStepSize_TightInclusion(
    const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    const Eigen::VectorXd& searchDir,
    double tolerance,
    const std::vector<std::pair<int, int>>& constraintSet,
    std::vector<std::pair<int, int>>& candidates,
    double& stepSize)
{
#if (CFL_FOR_CCD != 0)
    if (!constraintSet.size()) { return; }

#ifdef USE_TBB
    std::mutex stepSizeLock;
    tbb::parallel_for(0, (int)constraintSet.size(), 1, [&](int cI) {
#else
    for (int cI = 0; cI < constraintSet.size(); ++cI) {
#endif
        // #endif
        MMCVID MMCVIDI;
        if (constraintSet[cI].first < 0) {
            // PT
            MMCVIDI[0] = -mesh.SVI[-constraintSet[cI].first - 1] - 1;
            MMCVIDI[1] = mesh.SF(constraintSet[cI].second, 0);
            MMCVIDI[2] = mesh.SF(constraintSet[cI].second, 1);
            MMCVIDI[3] = mesh.SF(constraintSet[cI].second, 2);
        }
        else {
            // EE
            MMCVIDI[0] = mesh.SFEdges[constraintSet[cI].first].first;
            MMCVIDI[1] = mesh.SFEdges[constraintSet[cI].first].second;
            MMCVIDI[2] = mesh.SFEdges[constraintSet[cI].second].first;
            MMCVIDI[3] = mesh.SFEdges[constraintSet[cI].second].second;
        }

        if (MMCVIDI[0] >= 0) { // edge-edge
            double d_sqrt;
            computeEdgeEdgeD(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]),
                mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), d_sqrt);
            d_sqrt = std::sqrt(d_sqrt);
            if (d_sqrt == 0) {
                spdlog::error("Initial CCD distance is zero! Returning 0 stepSize.");
#ifdef USE_TBB
                std::scoped_lock lock(stepSizeLock);
#endif
                stepSize = 0;
                return;
            }

            double toi, output_tolerance;
            bool has_collision = inclusion_ccd::edgeEdgeCCD_double(
                mesh.V.row(MMCVIDI[0]).transpose(),
                mesh.V.row(MMCVIDI[1]).transpose(),
                mesh.V.row(MMCVIDI[2]).transpose(),
                mesh.V.row(MMCVIDI[3]).transpose(),
                mesh.V.row(MMCVIDI[0]).transpose() + searchDir.segment<dim>(MMCVIDI[0] * dim),
                mesh.V.row(MMCVIDI[1]).transpose() + searchDir.segment<dim>(MMCVIDI[1] * dim),
                mesh.V.row(MMCVIDI[2]).transpose() + searchDir.segment<dim>(MMCVIDI[2] * dim),
                mesh.V.row(MMCVIDI[3]).transpose() + searchDir.segment<dim>(MMCVIDI[3] * dim),
                /*err=*/tight_inclusion_ee_err,
                /*ms=*/std::min(TIGHT_INCLUSION_DIST_P * d_sqrt, TIGHT_INCLUSION_MIN_DIST),
                toi,
                tolerance,
                /*max_t=*/stepSize,
                /* max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                output_tolerance,
                /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);

            if (has_collision && toi < 1e-6) {
                has_collision = inclusion_ccd::edgeEdgeCCD_double(
                    mesh.V.row(MMCVIDI[0]).transpose(),
                    mesh.V.row(MMCVIDI[1]).transpose(),
                    mesh.V.row(MMCVIDI[2]).transpose(),
                    mesh.V.row(MMCVIDI[3]).transpose(),
                    mesh.V.row(MMCVIDI[0]).transpose() + searchDir.segment<dim>(MMCVIDI[0] * dim),
                    mesh.V.row(MMCVIDI[1]).transpose() + searchDir.segment<dim>(MMCVIDI[1] * dim),
                    mesh.V.row(MMCVIDI[2]).transpose() + searchDir.segment<dim>(MMCVIDI[2] * dim),
                    mesh.V.row(MMCVIDI[3]).transpose() + searchDir.segment<dim>(MMCVIDI[3] * dim),
                    /*err=*/tight_inclusion_ee_err,
                    /*ms=*/0,
                    toi,
                    tolerance,
                    /*max_t=*/stepSize,
                    /* max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                    output_tolerance,
                    /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                    /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);
                if (has_collision) {
                    toi *= 0.8;
                }
            }

            if (has_collision) {
#ifdef USE_TBB
                std::scoped_lock lock(stepSizeLock);
#endif
                if (toi < stepSize) {
                    stepSize = toi;
                }
            }
        }
        else { // point-triangle
            int vI = -MMCVIDI[0] - 1;
            assert(MMCVIDI[1] >= 0);

            double d_sqrt;
            computePointTriD(mesh.V.row(vI), mesh.V.row(MMCVIDI[1]),
                mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), d_sqrt);
            d_sqrt = std::sqrt(d_sqrt);
            if (d_sqrt == 0) {
                spdlog::error("Initial CCD distance is zero! Returning 0 stepSize.");
                stepSize = 0;
                return;
            }

            double toi, output_tolerance;
            bool has_collision = inclusion_ccd::vertexFaceCCD_double(
                mesh.V.row(vI).transpose(),
                mesh.V.row(MMCVIDI[1]).transpose(),
                mesh.V.row(MMCVIDI[2]).transpose(),
                mesh.V.row(MMCVIDI[3]).transpose(),
                mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                mesh.V.row(MMCVIDI[1]).transpose() + searchDir.segment<dim>(MMCVIDI[1] * dim),
                mesh.V.row(MMCVIDI[2]).transpose() + searchDir.segment<dim>(MMCVIDI[2] * dim),
                mesh.V.row(MMCVIDI[3]).transpose() + searchDir.segment<dim>(MMCVIDI[3] * dim),
                /*err=*/tight_inclusion_vf_err,
                /*ms=*/std::min(TIGHT_INCLUSION_DIST_P * d_sqrt, TIGHT_INCLUSION_MIN_DIST),
                toi,
                tolerance,
                /*max_t=*/stepSize,
                /* max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                output_tolerance,
                /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);

            if (has_collision && toi < 1e-6) {
                has_collision = inclusion_ccd::vertexFaceCCD_double(
                    mesh.V.row(vI).transpose(),
                    mesh.V.row(MMCVIDI[1]).transpose(),
                    mesh.V.row(MMCVIDI[2]).transpose(),
                    mesh.V.row(MMCVIDI[3]).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                    mesh.V.row(MMCVIDI[1]).transpose() + searchDir.segment<dim>(MMCVIDI[1] * dim),
                    mesh.V.row(MMCVIDI[2]).transpose() + searchDir.segment<dim>(MMCVIDI[2] * dim),
                    mesh.V.row(MMCVIDI[3]).transpose() + searchDir.segment<dim>(MMCVIDI[3] * dim),
                    /*err=*/tight_inclusion_vf_err,
                    /*ms=*/0,
                    toi,
                    tolerance,
                    /*max_t=*/stepSize,
                    /* max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                    output_tolerance,
                    /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                    /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);
                if (has_collision) {
                    toi *= 0.8;
                }
            }

            if (has_collision) {
#ifdef USE_TBB
                std::scoped_lock lock(stepSizeLock);
#endif
                if (toi < stepSize) {
                    stepSize = toi;
                }
            }
        }
    }
#ifdef USE_TBB
    );
#endif
#else
    largestFeasibleStepSize_CCD(mesh, sh, searchDir, tolerance, candidates, stepSize);
#endif
}
#endif // IPC_WITH_TIGHT_INCLUSION

template <int dim>
void SelfCollisionHandler<dim>::largestFeasibleStepSize_exact(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    const Eigen::VectorXd& searchDir,
    const ccd::CCDMethod method,
    const std::vector<std::pair<int, int>>& constraintSet,
    double& stepSize)
{
#if (CFL_FOR_CCD != 0)
    if (constraintSet.size()) {
        Eigen::VectorXd largestAlphasAS(constraintSet.size());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)constraintSet.size(), 1, [&](int cI)
#else
        for (int cI = 0; cI < constraintSet.size(); ++cI)
#endif
            {
                MMCVID MMCVIDI;
                if (constraintSet[cI].first < 0) {
                    // PT
                    MMCVIDI[0] = -mesh.SVI[-constraintSet[cI].first - 1] - 1;
                    MMCVIDI[1] = mesh.SF(constraintSet[cI].second, 0);
                    MMCVIDI[2] = mesh.SF(constraintSet[cI].second, 1);
                    MMCVIDI[3] = mesh.SF(constraintSet[cI].second, 2);
                }
                else {
                    // EE
                    MMCVIDI[0] = mesh.SFEdges[constraintSet[cI].first].first;
                    MMCVIDI[1] = mesh.SFEdges[constraintSet[cI].first].second;
                    MMCVIDI[2] = mesh.SFEdges[constraintSet[cI].second].first;
                    MMCVIDI[3] = mesh.SFEdges[constraintSet[cI].second].second;
                }

                if (MMCVIDI[0] >= 0) {
                    // edge-edge
                    largestAlphasAS[cI] = stepSize;
                    while (ccd::edgeEdgeCCD(
                        mesh.V.row(MMCVIDI[0]).transpose(),
                        mesh.V.row(MMCVIDI[1]).transpose(),
                        mesh.V.row(MMCVIDI[2]).transpose(),
                        mesh.V.row(MMCVIDI[3]).transpose(),
                        mesh.V.row(MMCVIDI[0]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[0] * dim),
                        mesh.V.row(MMCVIDI[1]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[1] * dim),
                        mesh.V.row(MMCVIDI[2]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[2] * dim),
                        mesh.V.row(MMCVIDI[3]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[3] * dim),
                        method)) {
                        largestAlphasAS[cI] /= 2.0;
                    }

#ifdef CHECK_RATIONAL_CCD
                    // largestAlphasAS[cI] is now a step size exact doesn't detect interpenetration
                    double inexactAlpha = 1.0;
                    if (CTCD::edgeEdgeCTCD(mesh.V.row(MMCVIDI[0]).transpose(),
                            mesh.V.row(MMCVIDI[1]).transpose(),
                            mesh.V.row(MMCVIDI[2]).transpose(),
                            mesh.V.row(MMCVIDI[3]).transpose(),
                            mesh.V.row(MMCVIDI[0]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[0] * dim),
                            mesh.V.row(MMCVIDI[1]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[1] * dim),
                            mesh.V.row(MMCVIDI[2]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[2] * dim),
                            mesh.V.row(MMCVIDI[3]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[3] * dim),
                            0.0, inexactAlpha)) {
                        std::cout << "exact EE false positive" << std::endl;
                    }
#endif // CHECK_RATIONAL_CCD
                }
                else {
                    // point-triangle
                    int vI = -MMCVIDI[0] - 1;
                    assert(MMCVIDI[1] >= 0);

                    largestAlphasAS[cI] = stepSize;
                    while (ccd::vertexFaceCCD(
                        mesh.V.row(vI).transpose(),
                        mesh.V.row(MMCVIDI[1]).transpose(),
                        mesh.V.row(MMCVIDI[2]).transpose(),
                        mesh.V.row(MMCVIDI[3]).transpose(),
                        mesh.V.row(vI).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(vI * dim),
                        mesh.V.row(MMCVIDI[1]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[1] * dim),
                        mesh.V.row(MMCVIDI[2]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[2] * dim),
                        mesh.V.row(MMCVIDI[3]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[3] * dim),
                        method)) {
                        largestAlphasAS[cI] /= 2.0;
                    }

#ifdef CHECK_RATIONAL_CCD
                    // largestAlphasAS[cI] is now a step size exact doesn't detect interpenetration
                    double inexactAlpha = 1.0;
                    if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                            mesh.V.row(MMCVIDI[1]).transpose(),
                            mesh.V.row(MMCVIDI[2]).transpose(),
                            mesh.V.row(MMCVIDI[3]).transpose(),
                            mesh.V.row(vI).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(vI * dim),
                            mesh.V.row(MMCVIDI[1]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[1] * dim),
                            mesh.V.row(MMCVIDI[2]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[2] * dim),
                            mesh.V.row(MMCVIDI[3]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(MMCVIDI[3] * dim),
                            0.0, inexactAlpha)) {
                        std::cout << "exact PT false positive" << std::endl;
                    }
#endif // CHECK_RATIONAL_CCD
                }
            }
#ifdef USE_TBB
        );
#endif
        stepSize = std::min(stepSize, largestAlphasAS.minCoeff());
    }
    return;
#endif

    largestFeasibleStepSize_CCD_exact(mesh, sh, searchDir, method, stepSize);
}

template <int dim>
void SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh, const Eigen::VectorXd& searchDir,
    double slackness, std::vector<std::pair<int, int>>& candidates,
    double& stepSize)
{
    const double CCDDistRatio = 1.0 - slackness;

    // point-point,edge,triangle
    timer_temp3.start(12);
    Eigen::VectorXd largestAlphasPPET(mesh.SVI.size());
#ifdef CCD_FILTERED_CS
    std::vector<std::vector<int>> PTCandidates(mesh.SVI.size());
#endif
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SVI.size(), 1, [&](int svI) {
#else
    for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
#endif
        int vI = mesh.SVI[svI];
        largestAlphasPPET[svI] = 1.0;
        int vICoDim = mesh.vICoDim(vI);

#ifdef USE_SH_LFSS
        std::unordered_set<int> sVInds, sEdgeInds, sTriInds;
        sh.queryPointForPrimitives(svI, sVInds, sEdgeInds, sTriInds);
        // NOTE: results may differ when computing step size with large eta as long-distance pairs are dropped
#endif

        // point-point
#ifdef USE_SH_LFSS
        for (const auto& svJ : sVInds) {
            if (svJ <= svI) { continue; }
#else
        for (int svJ = svI + 1; svJ < mesh.SVI.size(); ++svJ) {
#endif
            int vJ = mesh.SVI[svJ];
            if ((vICoDim < 3 && mesh.vICoDim(vJ) < 3) || (mesh.isDBCVertex(vI) && mesh.isDBCVertex(vJ))) {
                continue;
            }

            double largestAlpha = 1.0;
            double d_sqrt = (mesh.V.row(vI) - mesh.V.row(vJ)).norm();
            if (CTCD::vertexVertexCTCD(mesh.V.row(vI).transpose(),
                    mesh.V.row(vJ).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                    mesh.V.row(vJ).transpose() + searchDir.segment<dim>(vJ * dim),
                    CCDDistRatio * d_sqrt,
                    largestAlpha)) {
                if (largestAlpha < 1.0e-6) {
                    std::cout << "PP CCD tiny: " << vI << " " << vJ << std::endl;
                    if (!CTCD::vertexVertexCTCD(mesh.V.row(vI).transpose(),
                            mesh.V.row(vJ).transpose(),
                            mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                            mesh.V.row(vJ).transpose() + searchDir.segment<dim>(vJ * dim),
                            0.0, largestAlpha)) {
                        continue;
                    }
                    largestAlpha *= slackness;
                }
                if (largestAlpha < largestAlphasPPET[svI]) {
                    largestAlphasPPET[svI] = largestAlpha;
                }
            }
        }

        // point-edge
#ifdef USE_SH_LFSS
        for (const auto& seI : sEdgeInds) {
            const auto& meshEI = mesh.SFEdges[seI];
#else
        for (const auto& meshEI : mesh.SFEdges) {
#endif
            if (!(meshEI.first == vI || meshEI.second == vI)) {
                if ((vICoDim < 3 && mesh.vICoDim(meshEI.first) < 3) || (mesh.isDBCVertex(vI) && mesh.isDBCVertex(meshEI.first) && mesh.isDBCVertex(meshEI.second))) {
                    continue;
                }

                double d_sqrt;
                computePointEdgeD(mesh.V.row(vI), mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlpha = 1.0;
                if (CTCD::vertexEdgeCTCD(mesh.V.row(vI).transpose(),
                        mesh.V.row(meshEI.first).transpose(),
                        mesh.V.row(meshEI.second).transpose(),
                        mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                        mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                        mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                        CCDDistRatio * d_sqrt,
                        largestAlpha)) {
                    if (largestAlpha < 1.0e-6) {
                        std::cout << "PE CCD tiny: " << vI << " " << meshEI.first << " " << meshEI.second << std::endl;
                        if (!CTCD::vertexEdgeCTCD(mesh.V.row(vI).transpose(),
                                mesh.V.row(meshEI.first).transpose(),
                                mesh.V.row(meshEI.second).transpose(),
                                mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                                mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                                mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                                0.0, largestAlpha)) {
                            continue;
                        }
                        largestAlpha *= slackness;
                    }
                    if (largestAlpha < largestAlphasPPET[svI]) {
                        largestAlphasPPET[svI] = largestAlpha;
                    }
                }
            }
        }

        // point-triangle
#ifdef USE_SH_LFSS
        for (const auto& sfI : sTriInds) {
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI) {
#endif
            const RowVector3i& sfVInd = mesh.SF.row(sfI);
            if (!(vI == sfVInd[0] || vI == sfVInd[1] || vI == sfVInd[2])) {
                if ((vICoDim < 3 && mesh.sfICoDim(sfI) < 3) || (mesh.isDBCVertex(vI) && mesh.isDBCVertex(sfVInd[0]) && mesh.isDBCVertex(sfVInd[1]) && mesh.isDBCVertex(sfVInd[2]))) {
                    continue;
                }

                double d_sqrt;
                computePointTriD(mesh.V.row(vI), mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlpha = 1.0;
                if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose(),
                        mesh.V.row(sfVInd[1]).transpose(),
                        mesh.V.row(sfVInd[2]).transpose(),
                        mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                        mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                        mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                        mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                        CCDDistRatio * d_sqrt,
                        largestAlpha)) {
                    if (largestAlpha < 1.0e-6) {
                        std::cout << "PT CCD tiny: " << vI << " " << sfVInd[0] << " " << sfVInd[1] << " " << sfVInd[2] << std::endl;
                        if (!CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                                mesh.V.row(sfVInd[0]).transpose(),
                                mesh.V.row(sfVInd[1]).transpose(),
                                mesh.V.row(sfVInd[2]).transpose(),
                                mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                                mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                                mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                                mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                                0.0, largestAlpha)) {
                            continue;
                        }
                        largestAlpha *= slackness;
                    }
                    // activeSet_next.emplace_back(-vI - 1, sfVInd[0], sfVInd[2], sfVInd[1]);

                    if (largestAlpha < largestAlphasPPET[svI]) {
                        largestAlphasPPET[svI] = largestAlpha;
                    }
                }

#ifdef CHECK_RATIONAL_CCD
                if (ccd::vertexFaceCCD(mesh.V.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose(),
                        mesh.V.row(sfVInd[1]).transpose(),
                        mesh.V.row(sfVInd[2]).transpose(),
                        mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                        mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                        mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                        mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                        ccd::CCDMethod::RATIONAL_ROOT_PARITY)
                    && largestAlpha == 1.0) {
                    std::cout << "PT false negative type1" << std::endl;
                }

                if (ccd::vertexFaceCCD(mesh.V.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose(),
                        mesh.V.row(sfVInd[1]).transpose(),
                        mesh.V.row(sfVInd[2]).transpose(),
                        mesh.V.row(vI).transpose() + largestAlpha * searchDir.segment<dim>(vI * dim),
                        mesh.V.row(sfVInd[0]).transpose() + largestAlpha * searchDir.segment<dim>(sfVInd[0] * dim),
                        mesh.V.row(sfVInd[1]).transpose() + largestAlpha * searchDir.segment<dim>(sfVInd[1] * dim),
                        mesh.V.row(sfVInd[2]).transpose() + largestAlpha * searchDir.segment<dim>(sfVInd[2] * dim),
                        ccd::CCDMethod::RATIONAL_ROOT_PARITY)) {
                    std::cout << "PT false negative type2" << std::endl;
                }
#endif // CHECK_RATIONAL_CCD

#ifdef CCD_FILTERED_CS
                if (largestAlpha < 1.0) {
                    PTCandidates[svI].emplace_back(sfI);
                }
#endif

                // std::cout << "PT CCD before and after:" << std::endl;
                // std::cout << mesh.V.row(vI) << std::endl;
                // std::cout << mesh.V.row(sfVInd[0]) << std::endl;
                // std::cout << mesh.V.row(sfVInd[1]) << std::endl;
                // std::cout << mesh.V.row(sfVInd[2]) << std::endl;
                // std::cout << mesh.V.row(vI) + largestAlpha * searchDir.segment<dim>(vI * dim).transpose() << std::endl;
                // std::cout << mesh.V.row(sfVInd[0]) + largestAlpha * searchDir.segment<dim>(sfVInd[0] * dim).transpose() << std::endl;
                // std::cout << mesh.V.row(sfVInd[1]) + largestAlpha * searchDir.segment<dim>(sfVInd[1] * dim).transpose() << std::endl;
                // std::cout << mesh.V.row(sfVInd[2]) + largestAlpha * searchDir.segment<dim>(sfVInd[2] * dim).transpose() << std::endl;
            }
        }
    }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphasPPET.minCoeff());
    timer_temp3.stop();

    // #ifdef USE_SH_LFSS
    //     timer_temp3.start(11);
    //     sh.build(mesh, searchDir, stepSize, mesh.avgEdgeLen / 3.0); // rebuild with possibly smaller stepSize for less redundant candidates
    //     timer_temp3.stop();
    // #endif

    // edge-edge
    timer_temp3.start(13);
    Eigen::VectorXd largestAlphasEE(mesh.SFEdges.size());
#ifdef CCD_FILTERED_CS
    std::vector<std::vector<int>> EECandidates(mesh.SFEdges.size());
#endif
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SFEdges.size(), 1, [&](int eI)
#else
    for (int eI = 0; eI < mesh.SFEdges.size(); ++eI)
#endif
        {
            timer_mt.start(7);
            const auto& meshEI = mesh.SFEdges[eI];
            timer_mt.stop();
            int eICoDim = mesh.vICoDim(meshEI.first);

            largestAlphasEE[eI] = 1.0;
#ifdef USE_SH_LFSS
            std::unordered_set<int> sEdgeInds;
            timer_mt.start(3);
            sh.queryEdgeForEdgesWithBBoxCheck(mesh, searchDir, stepSize, eI, sEdgeInds);
            timer_mt.stop();
            // NOTE: results may differ when computing step size with large eta as long-distance pairs are dropped
            for (const auto& eJ : sEdgeInds) {
                timer_mt.start(7);
                const auto& meshEJ = mesh.SFEdges[eJ];
                timer_mt.stop();
#else
        for (int eJ = eI + 1; eJ < mesh.SFEdges.size(); ++eJ) {
            const auto& meshEJ = mesh.SFEdges[eJ];
#endif
                if (!(meshEI.first == meshEJ.first || meshEI.first == meshEJ.second || meshEI.second == meshEJ.first || meshEI.second == meshEJ.second || eI > eJ)) {
                    if ((eICoDim < 3 && mesh.vICoDim(meshEJ.first) < 3) || (mesh.isDBCVertex(meshEI.first) && mesh.isDBCVertex(meshEI.second) && mesh.isDBCVertex(meshEJ.first) && mesh.isDBCVertex(meshEJ.second))) {
                        continue;
                    }

                    double d_sqrt;
                    timer_mt.start(4);
                    computeEdgeEdgeD(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                        mesh.V.row(meshEJ.first), mesh.V.row(meshEJ.second), d_sqrt);
                    d_sqrt = std::sqrt(d_sqrt);
                    timer_mt.stop();

                    double largestAlpha = 1.0;
                    timer_mt.start(5);
                    if (CTCD::edgeEdgeCTCD(mesh.V.row(meshEI.first).transpose(),
                            mesh.V.row(meshEI.second).transpose(),
                            mesh.V.row(meshEJ.first).transpose(),
                            mesh.V.row(meshEJ.second).transpose(),
                            mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                            mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                            mesh.V.row(meshEJ.first).transpose() + searchDir.segment<dim>(meshEJ.first * dim),
                            mesh.V.row(meshEJ.second).transpose() + searchDir.segment<dim>(meshEJ.second * dim),
                            CCDDistRatio * d_sqrt,
                            largestAlpha)) {
                        if (largestAlpha < 1.0e-6) {
                            logFile << "EE CCD tiny: " << meshEI.first << " " << meshEI.second << " " << meshEJ.first << " " << meshEJ.second << std::endl;
                            // std::cout << "eta = " << CCDDistRatio * d_sqrt << std::endl;
                            // std::cout << mesh.V.row(meshEI.first) << std::endl;
                            // std::cout << mesh.V.row(meshEI.second) << std::endl;
                            // std::cout << mesh.V.row(meshEJ.first) << std::endl;
                            // std::cout << mesh.V.row(meshEJ.second) << std::endl;
                            // std::cout << mesh.V.row(meshEI.first) + searchDir.segment<dim>(meshEI.first * dim).transpose() << std::endl;
                            // std::cout << mesh.V.row(meshEI.second) + searchDir.segment<dim>(meshEI.second * dim).transpose() << std::endl;
                            // std::cout << mesh.V.row(meshEJ.first) + searchDir.segment<dim>(meshEJ.first * dim).transpose() << std::endl;
                            // std::cout << mesh.V.row(meshEJ.second) + searchDir.segment<dim>(meshEJ.second * dim).transpose() << std::endl;
                            if (!CTCD::edgeEdgeCTCD(mesh.V.row(meshEI.first).transpose(),
                                    mesh.V.row(meshEI.second).transpose(),
                                    mesh.V.row(meshEJ.first).transpose(),
                                    mesh.V.row(meshEJ.second).transpose(),
                                    mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                                    mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                                    mesh.V.row(meshEJ.first).transpose() + searchDir.segment<dim>(meshEJ.first * dim),
                                    mesh.V.row(meshEJ.second).transpose() + searchDir.segment<dim>(meshEJ.second * dim),
                                    0.0, largestAlpha)) {
                                continue;
                            }
                            largestAlpha *= slackness;
                        }
                        // MMCVID cVID(meshEI.first, meshEI.second, meshEJ.first, meshEJ.second);
                        // double c;
                        // compute_c<dim>(mesh.V.row(cVID[0]),
                        //     mesh.V.row(cVID[1]),
                        //     mesh.V.row(cVID[2]),
                        //     mesh.V.row(cVID[3]),
                        //     c, 1.0);
                        // if (c != 0.0) {
                        //     if (c < 0.0) {
                        //         int tmp = cVID[0];
                        //         cVID[0] = cVID[1];
                        //         cVID[1] = tmp;
                        //     }
                        //     activeSet_next.emplace_back(cVID);
                        // }

                        if (largestAlpha < largestAlphasEE[eI]) {
                            largestAlphasEE[eI] = largestAlpha;
                        }
                    }
                    timer_mt.stop();

#ifdef CHECK_RATIONAL_CCD
                    if (ccd::edgeEdgeCCD(mesh.V.row(meshEI.first).transpose(),
                                         mesh.V.row(meshEI.second).transpose(),
                                         mesh.V.row(meshEJ.first).transpose(),
                                         mesh.V.row(meshEJ.second).transpose(),
                                         mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                                         mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                                         mesh.V.row(meshEJ.first).transpose() + searchDir.segment<dim>(meshEJ.first * dim),
                                         mesh.V.row(meshEJ.second).transpose() + searchDir.segment<dim>(meshEJ.second * dim),
                                         ccd
                                         : CCDMethod::RATIONAL_ROOT_PARITY)
                        && largestAlpha == 1.0) {
                        std::cout << "EE false negative type1" << std::endl;
                    }

                    if (ccd::edgeEdgeCCD(mesh.V.row(meshEI.first).transpose(),
                                         mesh.V.row(meshEI.second).transpose(),
                                         mesh.V.row(meshEJ.first).transpose(),
                                         mesh.V.row(meshEJ.second).transpose(),
                                         mesh.V.row(meshEI.first).transpose() + largestAlpha * searchDir.segment<dim>(meshEI.first * dim),
                                         mesh.V.row(meshEI.second).transpose() + largestAlpha * searchDir.segment<dim>(meshEI.second * dim),
                                         mesh.V.row(meshEJ.first).transpose() + largestAlpha * searchDir.segment<dim>(meshEJ.first * dim),
                                         mesh.V.row(meshEJ.second).transpose() + largestAlpha * searchDir.segment<dim>(meshEJ.second * dim),
                                         ccd
                                         : CCDMethod::RATIONAL_ROOT_PARITY)) {
                        std::cout << "EE false negative type2" << std::endl;
                    }
#endif

#ifdef CCD_FILTERED_CS
                    if (largestAlpha < 1.0) {
                        EECandidates[eI].emplace_back(eJ);
                    }
#endif

                    // std::cout << "EE CCD before and after:" << std::endl;
                    // std::cout << mesh.V.row(meshEI.first) << std::endl;
                    // std::cout << mesh.V.row(meshEI.second) << std::endl;
                    // std::cout << mesh.V.row(meshEJ.first) << std::endl;
                    // std::cout << mesh.V.row(meshEJ.second) << std::endl;
                    // std::cout << mesh.V.row(meshEI.first) + largestAlpha * searchDir.segment<dim>(meshEI.first * dim).transpose() << std::endl;
                    // std::cout << mesh.V.row(meshEI.second) + largestAlpha * searchDir.segment<dim>(meshEI.second * dim).transpose() << std::endl;
                    // std::cout << mesh.V.row(meshEJ.first) + largestAlpha * searchDir.segment<dim>(meshEJ.first * dim).transpose() << std::endl;
                    // std::cout << mesh.V.row(meshEJ.second) + largestAlpha * searchDir.segment<dim>(meshEJ.second * dim).transpose() << std::endl;
                }
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphasEE.minCoeff());

#ifdef CCD_FILTERED_CS
    for (int svI = 0; svI < PTCandidates.size(); ++svI) {
        for (const auto& sfI : PTCandidates[svI]) {
            candidates.emplace_back(-svI - 1, sfI);
        }
    }
    for (int eI = 0; eI < EECandidates.size(); ++eI) {
        for (const auto& eJ : EECandidates[eI]) {
            candidates.emplace_back(eI, eJ);
        }
    }
    spdlog::info("# of new candidates {:d}", candidates.size());
#endif
    timer_temp3.stop();
}

#ifdef IPC_WITH_TIGHT_INCLUSION
template <int dim>
void SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD_TightInclusion(
    const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    const Eigen::VectorXd& searchDir,
    double tolerance,
    std::vector<std::pair<int, int>>& candidates,
    double& stepSize)
{
    timer_temp3.start(12);

#ifdef CCD_FILTERED_CS
    std::vector<std::vector<int>> PTCandidates(mesh.SVI.size());
#endif
#ifdef USE_TBB
    std::mutex stepSizeLock;
    tbb::parallel_for(0, (int)mesh.SVI.size(), 1, [&](int svI) {
#else
    for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
#endif
        int vI = mesh.SVI[svI];
        int vICoDim = mesh.vICoDim(vI);

#ifdef USE_SH_LFSS
        std::unordered_set<int> sVInds, sEdgeInds, sTriInds;
        sh.queryPointForPrimitives(svI, sVInds, sEdgeInds, sTriInds);
        // NOTE: results may differ when computing step size with large eta as long-distance pairs are dropped
#endif

        // point-triangle
#ifdef USE_SH_LFSS
        for (const auto& sfI : sTriInds) {
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI) {
#endif
            const RowVector3i& sfVInd = mesh.SF.row(sfI);
            if (!(vI == sfVInd[0] || vI == sfVInd[1] || vI == sfVInd[2])) {
                if ((vICoDim < 3 && mesh.sfICoDim(sfI) < 3)
                    || (mesh.isDBCVertex(vI) && mesh.isDBCVertex(sfVInd[0])
                        && mesh.isDBCVertex(sfVInd[1])
                        && mesh.isDBCVertex(sfVInd[2]))) {
                    continue;
                }

                double d_sqrt;
                computePointTriD(mesh.V.row(vI), mesh.V.row(sfVInd[0]),
                    mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);
                if (d_sqrt == 0) {
                    spdlog::error("Initial CCD distance is zero! Returning 0 stepSize.");
#ifdef USE_TBB
                    std::scoped_lock lock(stepSizeLock);
#endif
                    stepSize = 0;
                    return;
                }

                double toi, output_tolerance;
                bool has_collision = inclusion_ccd::vertexFaceCCD_double(
                    mesh.V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose(),
                    mesh.V.row(sfVInd[1]).transpose(),
                    mesh.V.row(sfVInd[2]).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                    mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                    mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                    mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                    /*err=*/tight_inclusion_vf_err,
                    /*ms=*/std::min(TIGHT_INCLUSION_DIST_P * d_sqrt, TIGHT_INCLUSION_MIN_DIST),
                    toi,
                    tolerance,
                    /*max_t=*/stepSize,
                    /* max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                    output_tolerance,
                    /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                    /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);

                if (has_collision && toi < 1e-6) {
                    has_collision = inclusion_ccd::vertexFaceCCD_double(
                        mesh.V.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose(),
                        mesh.V.row(sfVInd[1]).transpose(),
                        mesh.V.row(sfVInd[2]).transpose(),
                        mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                        mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                        mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                        mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                        /*err=*/tight_inclusion_vf_err,
                        /*ms=*/0,
                        toi,
                        tolerance,
                        /*max_t=*/stepSize,
                        /* max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                        output_tolerance,
                        /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                        /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);
                    if (has_collision) {
                        toi *= 0.8;
                    }
                }

                if (has_collision) {
#ifdef USE_TBB
                    std::scoped_lock lock(stepSizeLock);
#endif
                    if (toi < stepSize) {
                        stepSize = toi;
                    }
                }

#ifdef CCD_FILTERED_CS
                if (toi < 1.0) {
                    PTCandidates[svI].emplace_back(sfI);
                }
#endif
            }
        }
    }
#ifdef USE_TBB
    );
#endif
    timer_temp3.stop();

    // edge-edge
    timer_temp3.start(13);
#ifdef CCD_FILTERED_CS
    std::vector<std::vector<int>> EECandidates(mesh.SFEdges.size());
#endif
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SFEdges.size(), 1, [&](int eI) {
#else
    for (int eI = 0; eI < mesh.SFEdges.size(); ++eI) {
#endif
        timer_mt.start(7);
        const auto& meshEI = mesh.SFEdges[eI];
        timer_mt.stop();
        int eICoDim = mesh.vICoDim(meshEI.first);

#ifdef USE_SH_LFSS
        std::unordered_set<int> sEdgeInds;
        timer_mt.start(3);
        sh.queryEdgeForEdgesWithBBoxCheck(mesh, searchDir, stepSize, eI, sEdgeInds);
        timer_mt.stop();
        // NOTE: results may differ when computing step size with large eta as long-distance pairs are dropped
        for (const auto& eJ : sEdgeInds) {
            timer_mt.start(7);
            const auto& meshEJ = mesh.SFEdges[eJ];
            timer_mt.stop();
#else
        for (int eJ = eI + 1; eJ < mesh.SFEdges.size(); ++eJ) {
            const auto& meshEJ = mesh.SFEdges[eJ];
#endif
            if (!(meshEI.first == meshEJ.first
                    || meshEI.first == meshEJ.second
                    || meshEI.second == meshEJ.first
                    || meshEI.second == meshEJ.second
                    || eI > eJ)) {
                if ((eICoDim < 3 && mesh.vICoDim(meshEJ.first) < 3)
                    || (mesh.isDBCVertex(meshEI.first)
                        && mesh.isDBCVertex(meshEI.second)
                        && mesh.isDBCVertex(meshEJ.first)
                        && mesh.isDBCVertex(meshEJ.second))) {
                    continue;
                }

                double d_sqrt;
                timer_mt.start(4);
                computeEdgeEdgeD(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                    mesh.V.row(meshEJ.first), mesh.V.row(meshEJ.second), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);
                timer_mt.stop();
                if (d_sqrt == 0) {
                    spdlog::error("Initial CCD distance is zero! Returning 0 stepSize.");
#ifdef USE_TBB
                    std::scoped_lock lock(stepSizeLock);
#endif
                    stepSize = 0;
                    return;
                }

                double toi, output_tolerance;
                bool has_collision = inclusion_ccd::edgeEdgeCCD_double(
                    mesh.V.row(meshEI.first).transpose(),
                    mesh.V.row(meshEI.second).transpose(),
                    mesh.V.row(meshEJ.first).transpose(),
                    mesh.V.row(meshEJ.second).transpose(),
                    mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                    mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                    mesh.V.row(meshEJ.first).transpose() + searchDir.segment<dim>(meshEJ.first * dim),
                    mesh.V.row(meshEJ.second).transpose() + searchDir.segment<dim>(meshEJ.second * dim),
                    /*err=*/tight_inclusion_ee_err,
                    /*ms=*/std::min(TIGHT_INCLUSION_DIST_P * d_sqrt, TIGHT_INCLUSION_MIN_DIST),
                    toi,
                    tolerance,
                    /*max_t=*/stepSize,
                    /* max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                    output_tolerance,
                    /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                    /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);

                if (has_collision && toi < 1e-6) {
                    has_collision = inclusion_ccd::edgeEdgeCCD_double(
                        mesh.V.row(meshEI.first).transpose(),
                        mesh.V.row(meshEI.second).transpose(),
                        mesh.V.row(meshEJ.first).transpose(),
                        mesh.V.row(meshEJ.second).transpose(),
                        mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                        mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                        mesh.V.row(meshEJ.first).transpose() + searchDir.segment<dim>(meshEJ.first * dim),
                        mesh.V.row(meshEJ.second).transpose() + searchDir.segment<dim>(meshEJ.second * dim),
                        /*err=*/tight_inclusion_ee_err,
                        /*ms=*/0,
                        toi,
                        tolerance,
                        /*max_t=*/stepSize,
                        /* max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                        output_tolerance,
                        /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                        /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);
                    if (has_collision) {
                        toi *= 0.8;
                    }
                }

                if (has_collision) {
#ifdef USE_TBB
                    std::scoped_lock lock(stepSizeLock);
#endif
                    if (toi < stepSize) {
                        stepSize = toi;
                    }
                }

                timer_mt.stop();

#ifdef CCD_FILTERED_CS
                if (toi < 1.0) {
                    EECandidates[eI].emplace_back(eJ);
                }
#endif
            }
        }
    }
#ifdef USE_TBB
    );
#endif

#ifdef CCD_FILTERED_CS
    for (int svI = 0; svI < PTCandidates.size(); ++svI) {
        for (const auto& sfI : PTCandidates[svI]) {
            candidates.emplace_back(-svI - 1, sfI);
        }
    }
    for (int eI = 0; eI < EECandidates.size(); ++eI) {
        for (const auto& eJ : EECandidates[eI]) {
            candidates.emplace_back(eI, eJ);
        }
    }
    spdlog::info("# of new candidates {:d}", candidates.size());
#endif
    timer_temp3.stop();
}
#endif // IPC_WITH_TIGHT_INCLUSION

template <int dim>
void SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD_exact(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh, const Eigen::VectorXd& searchDir,
    const ccd::CCDMethod method, double& stepSize)
{
    // point-point,edge,triangle
    timer_temp3.start(12);
    Eigen::VectorXd largestAlphasPT(mesh.SVI.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SVI.size(), 1, [&](int svI)
#else
    for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
        {
            int vI = mesh.SVI[svI];
            largestAlphasPT[svI] = stepSize;
            int vICoDim = mesh.vICoDim(vI);

#ifdef USE_SH_LFSS
            std::unordered_set<int> sTriInds;
            sh.queryPointForTriangles(svI, sTriInds);
        // NOTE: results may differ when computing step size with large eta as long-distance pairs are dropped
#endif

        // point-triangle
#ifdef USE_SH_LFSS
            for (const auto& sfI : sTriInds)
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
            {
                const RowVector3i& sfVInd = mesh.SF.row(sfI);
                if (!(vI == sfVInd[0] || vI == sfVInd[1] || vI == sfVInd[2])) {
                    if ((vICoDim < 3 && mesh.sfICoDim(sfI) < 3) || (mesh.isDBCVertex(vI) && mesh.isDBCVertex(sfVInd[0]) && mesh.isDBCVertex(sfVInd[1]) && mesh.isDBCVertex(sfVInd[2]))) {
                        continue;
                    }

                    while (ccd::vertexFaceCCD(
                        mesh.V.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose(),
                        mesh.V.row(sfVInd[1]).transpose(),
                        mesh.V.row(sfVInd[2]).transpose(),
                        mesh.V.row(vI).transpose() + largestAlphasPT[svI] * searchDir.segment<dim>(vI * dim),
                        mesh.V.row(sfVInd[0]).transpose() + largestAlphasPT[svI] * searchDir.segment<dim>(sfVInd[0] * dim),
                        mesh.V.row(sfVInd[1]).transpose() + largestAlphasPT[svI] * searchDir.segment<dim>(sfVInd[1] * dim),
                        mesh.V.row(sfVInd[2]).transpose() + largestAlphasPT[svI] * searchDir.segment<dim>(sfVInd[2] * dim),
                        method)) {
                        largestAlphasPT[svI] /= 2.0;
                    }

#ifdef CHECK_RATIONAL_CCD
                    // largestAlphasPT[svI] is now a step size exact doesn't detect interpenetration
                    double inexactAlpha = 1.0;
                    if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                            mesh.V.row(sfVInd[0]).transpose(),
                            mesh.V.row(sfVInd[1]).transpose(),
                            mesh.V.row(sfVInd[2]).transpose(),
                            mesh.V.row(vI).transpose() + largestAlphasPT[svI] * searchDir.segment<dim>(vI * dim),
                            mesh.V.row(sfVInd[0]).transpose() + largestAlphasPT[svI] * searchDir.segment<dim>(sfVInd[0] * dim),
                            mesh.V.row(sfVInd[1]).transpose() + largestAlphasPT[svI] * searchDir.segment<dim>(sfVInd[1] * dim),
                            mesh.V.row(sfVInd[2]).transpose() + largestAlphasPT[svI] * searchDir.segment<dim>(sfVInd[2] * dim),
                            0.0, inexactAlpha)) {
                        std::cout << "exact PT false positive" << std::endl;
                    }
#endif // CHECK_RATIONAL_CCD
                }
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphasPT.minCoeff());
    timer_temp3.stop();

    // #ifdef USE_SH_LFSS
    //     timer_temp3.start(11);
    //     sh.build(mesh, searchDir, stepSize, mesh.avgEdgeLen / 3.0); // rebuild with possibly smaller stepSize for less redundant candidates
    //     timer_temp3.stop();
    // #endif

    // edge-edge
    timer_temp3.start(13);
    Eigen::VectorXd largestAlphasEE(mesh.SFEdges.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SFEdges.size(), 1, [&](int eI)
#else
    for (int eI = 0; eI < mesh.SFEdges.size(); ++eI)
#endif
        {
            timer_mt.start(7);
            const auto& meshEI = mesh.SFEdges[eI];
            timer_mt.stop();
            int eICoDim = mesh.vICoDim(meshEI.first);

            largestAlphasEE[eI] = stepSize;
#ifdef USE_SH_LFSS
            std::unordered_set<int> sEdgeInds;
            timer_mt.start(3);
            sh.queryEdgeForEdgesWithBBoxCheck(mesh, searchDir, stepSize, eI, sEdgeInds);
            timer_mt.stop();
            // NOTE: results may differ when computing step size with large eta as long-distance pairs are dropped
            for (const auto& eJ : sEdgeInds) {
                timer_mt.start(7);
                const auto& meshEJ = mesh.SFEdges[eJ];
                timer_mt.stop();
#else
        for (int eJ = eI + 1; eJ < mesh.SFEdges.size(); ++eJ) {
            const auto& meshEJ = mesh.SFEdges[eJ];
#endif
                if (!(meshEI.first == meshEJ.first || meshEI.first == meshEJ.second || meshEI.second == meshEJ.first || meshEI.second == meshEJ.second) || (eI > eJ)) {
                    if ((eICoDim < 3 && mesh.vICoDim(meshEJ.first) < 3) || (mesh.isDBCVertex(meshEI.first) && mesh.isDBCVertex(meshEI.second) && mesh.isDBCVertex(meshEJ.first) && mesh.isDBCVertex(meshEJ.second))) {
                        continue;
                    }

                    timer_mt.start(5);
                    while (ccd::edgeEdgeCCD(
                        mesh.V.row(meshEI.first).transpose(),
                        mesh.V.row(meshEI.second).transpose(),
                        mesh.V.row(meshEJ.first).transpose(),
                        mesh.V.row(meshEJ.second).transpose(),
                        mesh.V.row(meshEI.first).transpose() + largestAlphasEE[eI] * searchDir.segment<dim>(meshEI.first * dim),
                        mesh.V.row(meshEI.second).transpose() + largestAlphasEE[eI] * searchDir.segment<dim>(meshEI.second * dim),
                        mesh.V.row(meshEJ.first).transpose() + largestAlphasEE[eI] * searchDir.segment<dim>(meshEJ.first * dim),
                        mesh.V.row(meshEJ.second).transpose() + largestAlphasEE[eI] * searchDir.segment<dim>(meshEJ.second * dim),
                        method)) {
                        largestAlphasEE[eI] /= 2.0;
                    }
                    timer_mt.stop();

#ifdef CHECK_RATIONAL_CCD
                    // largestAlphasEE[eI] is now a step size exact doesn't detect interpenetration
                    double inexactAlpha = 1.0;
                    if (CTCD::edgeEdgeCTCD(mesh.V.row(meshEI.first).transpose(),
                            mesh.V.row(meshEI.second).transpose(),
                            mesh.V.row(meshEJ.first).transpose(),
                            mesh.V.row(meshEJ.second).transpose(),
                            mesh.V.row(meshEI.first).transpose() + largestAlphasEE[eI] * searchDir.segment<dim>(meshEI.first * dim),
                            mesh.V.row(meshEI.second).transpose() + largestAlphasEE[eI] * searchDir.segment<dim>(meshEI.second * dim),
                            mesh.V.row(meshEJ.first).transpose() + largestAlphasEE[eI] * searchDir.segment<dim>(meshEJ.first * dim),
                            mesh.V.row(meshEJ.second).transpose() + largestAlphasEE[eI] * searchDir.segment<dim>(meshEJ.second * dim),
                            0.0, inexactAlpha)) {
                        std::cout << "exact EE false positive" << std::endl;
                    }
#endif // CHECK_RATIONAL_CCD
                }
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphasEE.minCoeff());
    timer_temp3.stop();
}

template <int dim>
void SelfCollisionHandler<dim>::updateConstraints_QP(
    const Mesh<dim>& mesh,
    const std::vector<MMCVID>& activeSet,
    const CollisionConstraintType constraintType,
    const std::unordered_map<MMCVID, double, MMCVIDHash>& mmcvid_to_toi,
    std::vector<Eigen::Triplet<double>>& A_triplet, Eigen::VectorXd& l)
{
    if constexpr (dim == 3) {
        A_triplet.reserve(A_triplet.size() + activeSet.size() * 4 * dim);
        int constraintI = l.size();
        for (const auto& MMCVIDI : activeSet) {
            Eigen::Matrix<double, dim * 4, 1> gradc;

            double toi;
            try {
                toi = mmcvid_to_toi.at(MMCVIDI);
            }
            catch (std::out_of_range err) {
                spdlog::error("where=\"SelfCollisionHandler:evaluateConstraintsQP\" msg=\"Cannot compute the constraint value of MMCVID without toi!\"");
                spdlog::debug("MMCVIDI={:s}", MMCVIDI.str());
                spdlog::debug("mmcvid_to_toi:");
                for (auto const& [mmcvid, toi] : mmcvid_to_toi) {
                    spdlog::debug("{:s}: {:g}", mmcvid.str(), toi);
                }
                throw "missing toi for constraint";
            }

            if (MMCVIDI[0] >= 0) {
                // edge-edge ++++
                compute_collision_constraint_gradient(
                    // first edge at t = 0
                    mesh.V_prev.row(MMCVIDI[0]), mesh.V_prev.row(MMCVIDI[1]),
                    // second edge at t = 0
                    mesh.V_prev.row(MMCVIDI[2]), mesh.V_prev.row(MMCVIDI[3]),
                    // first edge at t = 1
                    mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]),
                    // second edge at t = 1
                    mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]),
                    constraintType, /*is_edge_edge=*/true, toi,
                    gradc);

                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < dim; j++) {
                        A_triplet.emplace_back(
                            constraintI,
                            MMCVIDI[i] * dim + j,
                            gradc[i * dim + j]);
                    }
                }
            }
            else {
                // point-triangle -+++
                assert(MMCVIDI[1] >= 0); // no triangle-point self collisions

                compute_collision_constraint_gradient(
                    // point at t = 0
                    mesh.V_prev.row(-MMCVIDI[0] - 1),
                    // triangle at t = 0
                    mesh.V_prev.row(MMCVIDI[1]),
                    mesh.V_prev.row(MMCVIDI[2]),
                    mesh.V_prev.row(MMCVIDI[3]),
                    // point at t = 1
                    mesh.V.row(-MMCVIDI[0] - 1),
                    // triangle at t = 1
                    mesh.V.row(MMCVIDI[1]),
                    mesh.V.row(MMCVIDI[2]),
                    mesh.V.row(MMCVIDI[3]),
                    constraintType, /*is_edge_edge=*/false, toi,
                    gradc);

                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < dim; j++) {
                        A_triplet.emplace_back(
                            constraintI,
                            (i == 0 ? (-MMCVIDI[i] - 1) : (MMCVIDI[i])) * dim + j,
                            gradc[i * dim + j]);
                    }
                }
            }

            ++constraintI;
        }

        // Use coef=-1 to move the linear component to the right side of the
        // inequality:  ∇g(x₀) * x + g(x₀) ≥ ϵ ≡ ∇g(x₀) * x ≥ -g(x₀) + ϵ
        SelfCollisionHandler<dim>::evaluateConstraintsQP(
            mesh, activeSet, constraintType, mmcvid_to_toi, l, /*coef=*/-1.0);
    }
    else {
        // TODO: 2D collision
    }
}

template <int dim>
void SelfCollisionHandler<dim>::filterSearchDir_QP(
    const Mesh<dim>& mesh,
    Eigen::VectorXd& searchDir,
    std::vector<MMCVID>& activeSet_next)
{
}

// Update the active set of constraints.
template <int dim>
bool SelfCollisionHandler<dim>::updateActiveSet_QP(
    const Mesh<dim>& mesh,
    const Eigen::VectorXd& searchDir,
    const CollisionConstraintType constraintType,
    std::vector<MMCVID>& activeSet,
    std::unordered_map<MMCVID, double, MMCVIDHash>& mmcvid_to_toi,
    const ccd::CCDMethod ccdMethod,
    const double eta,
    const double ccd_tol)
{
    bool newConstraintsAdded = false;
    mmcvid_to_toi.clear();
    std::unordered_set<MMCVID, MMCVIDHash> prevActiveSet;
    for (const auto& mmcvid : activeSet) {
        prevActiveSet.insert(mmcvid);
        mmcvid_to_toi.insert_or_assign(mmcvid, std::numeric_limits<double>::infinity());
    }
    // Verschoor does not clear the previous active set every iteration.
    // At the start of the time-step, make sure to clear the active set.
    // bool wasActiveSetCleared = constraintType != CollisionConstraintType::VERSCHOOR && constraintType != CollisionConstraintType::CMR;
    bool wasActiveSetCleared = false;
    if (wasActiveSetCleared) {
        activeSet.clear();
    }

#ifdef USE_SH_CCS
    Eigen::MatrixXd searchDir0 = mesh.V - mesh.V_prev;
    searchDir0 = Eigen::Map<Eigen::VectorXd>(searchDir0.data(), searchDir0.size());
    searchDir0 += searchDir;
    SpatialHash<dim> sh(mesh, searchDir0, 1, mesh.avgEdgeLen / 3.0, true);
#endif

    // point-triangle
    // Loop over mesh surface vertices
    for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
        int vI = mesh.SVI[svI];
        // Loop over mesh surface faces
#ifdef USE_SH_CCS
        std::unordered_set<int> sTriInds;
        sh.queryPointForTriangles(
            mesh.V_prev.row(vI),
            mesh.V.row(vI) + searchDir.segment<dim>(vI * dim).transpose() - mesh.V_prev.row(vI),
            sTriInds,
            /*radius=*/eta);
        for (const int& sfI : sTriInds) {
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI) {
#endif
            const RowVector3i& sfVInd = mesh.SF.row(sfI);
            if (vI == sfVInd[0] || vI == sfVInd[1] || vI == sfVInd[2]) {
                continue; // Skip triangles that contain the point
            }

            double toi;
            bool intersects;
            switch (ccdMethod) {
            case ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER: {
                intersects = CTCD::vertexFaceCTCD(
                    mesh.V_prev.row(vI).transpose(),
                    mesh.V_prev.row(sfVInd[0]).transpose(),
                    mesh.V_prev.row(sfVInd[1]).transpose(),
                    mesh.V_prev.row(sfVInd[2]).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                    mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                    mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                    mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                    eta,
                    toi);
                break;
            }
            case ccd::CCDMethod::TIGHT_INCLUSION: {
#ifdef IPC_WITH_TIGHT_INCLUSION
                double output_tolerance;
                intersects = inclusion_ccd::vertexFaceCCD_double(
                    mesh.V_prev.row(vI).transpose(),
                    mesh.V_prev.row(sfVInd[0]).transpose(),
                    mesh.V_prev.row(sfVInd[1]).transpose(),
                    mesh.V_prev.row(sfVInd[2]).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                    mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                    mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                    mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                    tight_inclusion_vf_err,
                    eta,
                    toi,
                    /*tolerance=*/ccd_tol,
                    /*max_t=*/1,
                    /* max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                    output_tolerance,
                    /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                    /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);
                break;
#else
                spdlog::error("Tight Inclusion CCD is disabled in CMake (CCD_WRAPPER_WITH_TIGHT_INCLUSION=OFF)!");
                throw std::runtime_error("Tight Inclusion CCD is disabled in CMake (CCD_WRAPPER_WITH_TIGHT_INCLUSION=OFF)!");
#endif
            }
            default:
                intersects = vertexFaceToIBisection(
                    mesh.V_prev.row(vI).transpose(),
                    mesh.V_prev.row(sfVInd[0]).transpose(),
                    mesh.V_prev.row(sfVInd[1]).transpose(),
                    mesh.V_prev.row(sfVInd[2]).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                    mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                    mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                    mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                    ccdMethod,
                    toi);
            }

            MMCVID mmcvid(
                -vI - 1, // mesh point
                sfVInd[0], sfVInd[2], sfVInd[1]); // mesh triangle

            mmcvid_to_toi.insert_or_assign(mmcvid,
                intersects ? toi : std::numeric_limits<double>::infinity());
            bool isNewConstraint = prevActiveSet.find(mmcvid) == prevActiveSet.end();
            if (intersects && (wasActiveSetCleared || isNewConstraint)) {
                newConstraintsAdded |= isNewConstraint;
                newConstraintsAdded = true;
                activeSet.push_back(mmcvid);
            }
        }
    }

    // edge-edge
    for (int i = 0; i < mesh.SFEdges.size(); i++) {
        const std::pair<int, int>& edge1 = mesh.SFEdges[i];
        Eigen::Matrix<double, dim, 1> edge1_dir = (mesh.V.row(edge1.second) - mesh.V.row(edge1.first)).normalized();
#ifdef USE_SH_CCS
        std::vector<int> sEdgeInds;
        sh.queryEdgeForEdges(
            mesh.V_prev.row(edge1.first).transpose(),
            mesh.V_prev.row(edge1.second).transpose(),
            mesh.V.row(edge1.first).transpose() + searchDir.segment<dim>(edge1.first * dim) - mesh.V_prev.row(edge1.first).transpose(),
            mesh.V.row(edge1.second).transpose() + searchDir.segment<dim>(edge1.second * dim) - mesh.V_prev.row(edge1.second).transpose(),
            sEdgeInds, eta, i);
        for (const int& j : sEdgeInds) {
#else
        for (int j = i + 1; j < mesh.SFEdges.size(); j++) {
#endif
            const std::pair<int, int>& edge2 = mesh.SFEdges[j];
            if (edge1.first == edge2.first || edge1.first == edge2.second
                || edge1.second == edge2.first || edge1.second == edge2.second) {
                continue;
            }

            Eigen::Matrix<double, dim, 1> edge2_dir = (mesh.V.row(edge2.second) - mesh.V.row(edge2.first)).normalized();
            // Skip parallel edges
            const double tol = 1e-12;
            if (abs(edge1_dir.dot(edge2_dir)) > 1 - tol) {
                continue;
            }

            double toi;
            bool intersects;
            switch (ccdMethod) {
            case ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER:
                intersects = CTCD::edgeEdgeCTCD(
                    mesh.V_prev.row(edge1.first).transpose(),
                    mesh.V_prev.row(edge1.second).transpose(),
                    mesh.V_prev.row(edge2.first).transpose(),
                    mesh.V_prev.row(edge2.second).transpose(),
                    mesh.V.row(edge1.first).transpose() + searchDir.segment<dim>(edge1.first * dim),
                    mesh.V.row(edge1.second).transpose() + searchDir.segment<dim>(edge1.second * dim),
                    mesh.V.row(edge2.first).transpose() + searchDir.segment<dim>(edge2.first * dim),
                    mesh.V.row(edge2.second).transpose() + searchDir.segment<dim>(edge2.first * dim),
                    eta,
                    toi);
                break;
            case ccd::CCDMethod::TIGHT_INCLUSION: {
#ifdef IPC_WITH_TIGHT_INCLUSION
                double output_tolerance;
                intersects = inclusion_ccd::edgeEdgeCCD_double(
                    mesh.V_prev.row(edge1.first).transpose(),
                    mesh.V_prev.row(edge1.second).transpose(),
                    mesh.V_prev.row(edge2.first).transpose(),
                    mesh.V_prev.row(edge2.second).transpose(),
                    mesh.V.row(edge1.first).transpose() + searchDir.segment<dim>(edge1.first * dim),
                    mesh.V.row(edge1.second).transpose() + searchDir.segment<dim>(edge1.second * dim),
                    mesh.V.row(edge2.first).transpose() + searchDir.segment<dim>(edge2.first * dim),
                    mesh.V.row(edge2.second).transpose() + searchDir.segment<dim>(edge2.first * dim),
                    tight_inclusion_ee_err,
                    eta,
                    toi,
                    /*tolerance=*/ccd_tol,
                    /*max_t=*/1,
                    /* max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                    output_tolerance,
                    /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                    /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);
                break;
#else
                spdlog::error("Tight Inclusion CCD is disabled in CMake (CCD_WRAPPER_WITH_TIGHT_INCLUSION=OFF)!");
                throw std::runtime_error("Tight Inclusion CCD is disabled in CMake (CCD_WRAPPER_WITH_TIGHT_INCLUSION=OFF)!");
#endif
            }
            default:
                intersects = edgeEdgeToIBisection(
                    mesh.V_prev.row(edge1.first).transpose(),
                    mesh.V_prev.row(edge1.second).transpose(),
                    mesh.V_prev.row(edge2.first).transpose(),
                    mesh.V_prev.row(edge2.second).transpose(),
                    mesh.V.row(edge1.first).transpose() + searchDir.segment<dim>(edge1.first * dim),
                    mesh.V.row(edge1.second).transpose() + searchDir.segment<dim>(edge1.second * dim),
                    mesh.V.row(edge2.first).transpose() + searchDir.segment<dim>(edge2.first * dim),
                    mesh.V.row(edge2.second).transpose() + searchDir.segment<dim>(edge2.first * dim),
                    ccdMethod,
                    toi);
            }

            if (intersects && constraintType == CollisionConstraintType::VERSCHOOR) {
                Eigen::Vector3d v0_toi = (mesh.V.row(edge1.first) - mesh.V_prev.row(edge1.first)) * toi + mesh.V_prev.row(edge1.first);
                Eigen::Vector3d v1_toi = (mesh.V.row(edge1.second) - mesh.V_prev.row(edge1.second)) * toi + mesh.V_prev.row(edge1.second);
                Eigen::Vector3d v2_toi = (mesh.V.row(edge2.first) - mesh.V_prev.row(edge2.first)) * toi + mesh.V_prev.row(edge2.first);
                Eigen::Vector3d v3_toi = (mesh.V.row(edge2.second) - mesh.V_prev.row(edge2.second)) * toi + mesh.V_prev.row(edge2.second);
                Eigen::Vector3d edge1_dir_toi = (v1_toi - v0_toi).normalized();
                Eigen::Vector3d edge2_dir_toi = (v3_toi - v2_toi).normalized();
                if (abs(edge1_dir_toi.dot(edge2_dir_toi)) > 1 - tol) {
                    continue;
                }
            }

            MMCVID mmcvid;
            // Check the orientation of the edge.
            double test_volume;
            compute_c<dim>(
                mesh.V_prev.row(edge1.first),
                mesh.V_prev.row(edge1.second),
                mesh.V_prev.row(edge2.first),
                mesh.V_prev.row(edge2.second),
                test_volume, 1.0);
            if (test_volume >= 0) {
                mmcvid = MMCVID(
                    edge1.first, edge1.second,
                    edge2.first, edge2.second);
            }
            else {
                mmcvid = MMCVID(
                    edge1.second, edge1.first,
                    edge2.first, edge2.second);
            }

            mmcvid_to_toi.insert_or_assign(mmcvid,
                intersects ? toi : std::numeric_limits<double>::infinity());
            bool isNewConstraint = prevActiveSet.find(mmcvid) == prevActiveSet.end();
            if (intersects && (wasActiveSetCleared || isNewConstraint)) {
                newConstraintsAdded |= isNewConstraint;
                activeSet.push_back(mmcvid);
            }
        }
    }

    return newConstraintsAdded;
}

template <int dim>
void SelfCollisionHandler<dim>::computeConstraintSet(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    double dHat, std::vector<MMCVID>& constraintSet,
    std::vector<MMCVID>& paraEEMMCVIDSet,
    std::vector<std::pair<int, int>>& paraEEeIeJSet,
    bool getPTEE, std::vector<std::pair<int, int>>& cs_PTEE)
{
#ifdef USE_SH_CCS
    double sqrtDHat = std::sqrt(dHat);
#endif

    // point-triangle
    timer_temp3.start(8);
    std::vector<std::vector<MMCVID>> constraintSetPT(mesh.SVI.size());
    std::vector<std::vector<int>> cs_PT;
    if (getPTEE) {
        cs_PT.resize(mesh.SVI.size());
    }
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SVI.size(), 1, [&](int svI)
#else
    for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
        {
            int vI = mesh.SVI[svI];
            int vICoDim = mesh.vICoDim(vI);
#ifdef USE_SH_CCS
            std::unordered_set<int> triInds; // NOTE: different constraint order will result in numerically different results
            sh.queryPointForTriangles(mesh.V.row(vI), sqrtDHat, triInds);
            for (const auto& sfI : triInds)
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
            {
                const RowVector3i& sfVInd = mesh.SF.row(sfI);
                if (!(vI == sfVInd[0] || vI == sfVInd[1] || vI == sfVInd[2])) {
                    if ((vICoDim < 3 && mesh.sfICoDim(sfI) < 3) || (mesh.isDBCVertex(vI) && mesh.isDBCVertex(sfVInd[0]) && mesh.isDBCVertex(sfVInd[1]) && mesh.isDBCVertex(sfVInd[2]))) {
                        continue;
                    }

                    int dtype = dType_PT(mesh.V.row(vI), mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]));
                    double d;
                    switch (dtype) {
                    case 0: {
                        d_PP(mesh.V.row(vI), mesh.V.row(sfVInd[0]), d);
                        if (d < dHat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[0], -1, -1);
                        }
                        break;
                    }

                    case 1: {
                        d_PP(mesh.V.row(vI), mesh.V.row(sfVInd[1]), d);
                        if (d < dHat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[1], -1, -1);
                        }
                        break;
                    }

                    case 2: {
                        d_PP(mesh.V.row(vI), mesh.V.row(sfVInd[2]), d);
                        if (d < dHat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[2], -1, -1);
                        }
                        break;
                    }

                    case 3: {
                        d_PE(mesh.V.row(vI), mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), d);
                        if (d < dHat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[0], sfVInd[1], -1);
                        }
                        break;
                    }

                    case 4: {
                        d_PE(mesh.V.row(vI), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), d);
                        if (d < dHat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[1], sfVInd[2], -1);
                        }
                        break;
                    }

                    case 5: {
                        d_PE(mesh.V.row(vI), mesh.V.row(sfVInd[2]), mesh.V.row(sfVInd[0]), d);
                        if (d < dHat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[2], sfVInd[0], -1);
                        }
                        break;
                    }

                    case 6: {
                        d_PT(mesh.V.row(vI), mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), d);
                        if (d < dHat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[0], sfVInd[1], sfVInd[2]);
                        }
                        break;
                    }

                    default:
                        break;
                    }

                    if (getPTEE && d < dHat) {
                        cs_PT[svI].emplace_back(sfI);
                    }
                }
            }
        }
#ifdef USE_TBB
    );
#endif
    timer_temp3.stop();

    // edge-edge
    timer_temp3.start(9);
    std::vector<std::vector<MMCVID>> constraintSetEE(mesh.SFEdges.size());
    std::vector<std::vector<int>> cs_EE;
    if (getPTEE) {
        cs_EE.resize(mesh.SFEdges.size());
    }
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SFEdges.size(), 1, [&](int eI)
#else
    for (int eI = 0; eI < mesh.SFEdges.size(); ++eI)
#endif
        {
            timer_mt.start(6);
            const auto& meshEI = mesh.SFEdges[eI];
            timer_mt.stop();
            int eICoDim = mesh.vICoDim(meshEI.first);

#ifdef USE_SH_CCS
            std::vector<int> edgeInds; // NOTE: different constraint order will result in numerically different results
            // timer_mt.start(0);
            sh.queryEdgeForEdgesWithBBoxCheck(mesh, mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), sqrtDHat, edgeInds, eI);
            // timer_mt.stop();
            timer_mt.start(23);
            for (const auto& eJ : edgeInds) {
#else
        for (int eJ = 0; eJ < mesh.SFEdges.size(); eJ++) {
#endif
                timer_mt.start(6);
                const auto& meshEJ = mesh.SFEdges[eJ];
                timer_mt.stop();
                if (!(meshEI.first == meshEJ.first || meshEI.first == meshEJ.second || meshEI.second == meshEJ.first || meshEI.second == meshEJ.second || eI > eJ)) {
                    if ((eICoDim < 3 && mesh.vICoDim(meshEJ.first) < 3) || (mesh.isDBCVertex(meshEI.first) && mesh.isDBCVertex(meshEI.second) && mesh.isDBCVertex(meshEJ.first) && mesh.isDBCVertex(meshEJ.second))) {
                        continue;
                    }

                    timer_mt.start(1);
                    int dtype = dType_EE(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), mesh.V.row(meshEJ.first), mesh.V.row(meshEJ.second));
                    timer_mt.stop();

                    timer_mt.start(22);
                    double EECrossSqNorm, eps_x;
                    computeEECrossSqNorm(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), mesh.V.row(meshEJ.first), mesh.V.row(meshEJ.second), EECrossSqNorm);
                    compute_eps_x(mesh, meshEI.first, meshEI.second, meshEJ.first, meshEJ.second, eps_x);
                    // int add_e = -1;
                    int add_e = (EECrossSqNorm < eps_x) ? -eJ - 2 : -1;
                    // == -1: regular,
                    // <= -2 && >= -mesh.SFEdges.size()-1: nearly parallel PP or PE,
                    // <= -mesh.SFEdges.size()-2: nearly parallel EE

                    double d;
                    timer_mt.start(2);
                    switch (dtype) {
                    case 0: {
                        d_PP(mesh.V.row(meshEI.first), mesh.V.row(meshEJ.first), d);
                        if (d < dHat) {
                            constraintSetEE[eI].emplace_back(-meshEI.first - 1, meshEJ.first, -1, add_e);
                        }
                        break;
                    }

                    case 1: {
                        d_PP(mesh.V.row(meshEI.first), mesh.V.row(meshEJ.second), d);
                        if (d < dHat) {
                            constraintSetEE[eI].emplace_back(-meshEI.first - 1, meshEJ.second, -1, add_e);
                        }
                        break;
                    }

                    case 2: {
                        d_PE(mesh.V.row(meshEI.first), mesh.V.row(meshEJ.first), mesh.V.row(meshEJ.second), d);
                        if (d < dHat) {
                            constraintSetEE[eI].emplace_back(-meshEI.first - 1, meshEJ.first, meshEJ.second, add_e);
                        }
                        break;
                    }

                    case 3: {
                        d_PP(mesh.V.row(meshEI.second), mesh.V.row(meshEJ.first), d);
                        if (d < dHat) {
                            constraintSetEE[eI].emplace_back(-meshEI.second - 1, meshEJ.first, -1, add_e);
                        }
                        break;
                    }

                    case 4: {
                        d_PP(mesh.V.row(meshEI.second), mesh.V.row(meshEJ.second), d);
                        if (d < dHat) {
                            constraintSetEE[eI].emplace_back(-meshEI.second - 1, meshEJ.second, -1, add_e);
                        }
                        break;
                    }

                    case 5: {
                        d_PE(mesh.V.row(meshEI.second), mesh.V.row(meshEJ.first), mesh.V.row(meshEJ.second), d);
                        if (d < dHat) {
                            constraintSetEE[eI].emplace_back(-meshEI.second - 1, meshEJ.first, meshEJ.second, add_e);
                        }
                        break;
                    }

                    case 6: {
                        d_PE(mesh.V.row(meshEJ.first), mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), d);
                        if (d < dHat) {
                            constraintSetEE[eI].emplace_back(-meshEJ.first - 1, meshEI.first, meshEI.second, add_e);
                        }
                        break;
                    }

                    case 7: {
                        d_PE(mesh.V.row(meshEJ.second), mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), d);
                        if (d < dHat) {
                            constraintSetEE[eI].emplace_back(-meshEJ.second - 1, meshEI.first, meshEI.second, add_e);
                        }
                        break;
                    }

                    case 8: {
                        d_EE(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), mesh.V.row(meshEJ.first), mesh.V.row(meshEJ.second), d);
                        if (d < dHat) {
                            if (add_e <= -2) {
                                constraintSetEE[eI].emplace_back(meshEI.first, meshEI.second, meshEJ.first, -meshEJ.second - mesh.SFEdges.size() - 2);
                            }
                            else {
                                constraintSetEE[eI].emplace_back(meshEI.first, meshEI.second, meshEJ.first, meshEJ.second);
                            }
                        }
                        break;
                    }

                    default:
                        break;
                    }

                    if (getPTEE && d < dHat) {
                        cs_EE[eI].emplace_back(eJ);
                    }
                    timer_mt.stop();
                }
                timer_mt.start(23);
            }
        }
#ifdef USE_TBB
    );
#endif
    timer_temp3.stop();

    timer_temp3.start(10);
    if (getPTEE) {
        cs_PTEE.resize(0);
        cs_PTEE.reserve(cs_PT.size() + cs_EE.size());
        for (int svI = 0; svI < cs_PT.size(); ++svI) {
            for (const auto& sfI : cs_PT[svI]) {
                cs_PTEE.emplace_back(-svI - 1, sfI);
            }
        }
        for (int eI = 0; eI < cs_EE.size(); ++eI) {
            for (const auto& eJ : cs_EE[eI]) {
                cs_PTEE.emplace_back(eI, eJ);
            }
        }
    }

    constraintSet.resize(0);
    constraintSet.reserve(constraintSetPT.size() + constraintSetEE.size());
    // for (const auto& csI : constraintSetPT) {
    //     constraintSet.insert(constraintSet.end(), csI.begin(), csI.end());
    // }
    // for (const auto& csI : constraintSetEE) {
    //     constraintSet.insert(constraintSet.end(), csI.begin(), csI.end());
    // }
    std::map<MMCVID, int> constraintCounter;
    for (const auto& csI : constraintSetPT) {
        for (const auto& cI : csI) {
            if (cI[3] < 0) {
                // PP or PE
                ++constraintCounter[cI];
            }
            else {
                constraintSet.emplace_back(cI);
            }
        }
    }
    paraEEMMCVIDSet.resize(0);
    paraEEeIeJSet.resize(0);
    int eI = 0;
    for (const auto& csI : constraintSetEE) {
        for (const auto& cI : csI) {
            if (cI[3] >= 0) {
                // regular EE
                constraintSet.emplace_back(cI);
            }
            else if (cI[3] == -1) {
                // regular PP or PE
                ++constraintCounter[cI];
            }
            else if (cI[3] >= -mesh.SFEdges.size() - 1) {
                // nearly parallel PP or PE
                paraEEMMCVIDSet.emplace_back(cI[0], cI[1], cI[2], -1);
                paraEEeIeJSet.emplace_back(eI, -cI[3] - 2);
            }
            else {
                // nearly parallel EE
                paraEEMMCVIDSet.emplace_back(cI[0], cI[1], cI[2], -cI[3] - mesh.SFEdges.size() - 2);
                paraEEeIeJSet.emplace_back(-1, -1);
            }
        }
        ++eI;
    }

    constraintSet.reserve(constraintSet.size() + constraintCounter.size());
    for (const auto& ccI : constraintCounter) {
        constraintSet.emplace_back(MMCVID(ccI.first[0], ccI.first[1], ccI.first[2], -ccI.second));
    }
    timer_temp3.stop();
}

template <int dim>
void SelfCollisionHandler<dim>::computeDistCoordAndTanBasis(
    const Mesh<dim>& mesh,
    const std::vector<MMCVID>& constraintSet,
    std::vector<Eigen::Vector2d>& MMDistCoord,
    std::vector<Eigen::Matrix<double, 3, 2>>& MMTanBasis)
{
    // TODO: parallelize
    MMDistCoord.resize(constraintSet.size());
    MMTanBasis.resize(constraintSet.size());
    for (int cI = 0; cI < constraintSet.size(); ++cI) {
        const auto& MMCVIDI = constraintSet[cI];
        if (MMCVIDI[0] >= 0) {
            // edge-edge
            computeClosestPoint_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]),
                mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), MMDistCoord[cI]);
            computeTangentBasis_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]),
                mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), MMTanBasis[cI]);
        }
        else {
            // point-triangle and degenerate edge-edge
            assert(MMCVIDI[1] >= 0);
            if (MMCVIDI[2] < 0) {
                // PP
                MMDistCoord[cI].setZero(); // Store something instead of random memory
                computeTangentBasis_PP(mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(MMCVIDI[1]),
                    MMTanBasis[cI]);
            }
            else if (MMCVIDI[3] < 0) {
                // PE
                computeClosestPoint_PE(mesh.V.row(-MMCVIDI[0] - 1),
                    mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), MMDistCoord[cI][0]);
                MMDistCoord[cI][1] = 0; // Store something instead of random memory
                computeTangentBasis_PE(mesh.V.row(-MMCVIDI[0] - 1),
                    mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), MMTanBasis[cI]);
            }
            else {
                // PT
                computeClosestPoint_PT(mesh.V.row(-MMCVIDI[0] - 1),
                    mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]),
                    MMDistCoord[cI]);
                computeTangentBasis_PT(mesh.V.row(-MMCVIDI[0] - 1),
                    mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]),
                    MMTanBasis[cI]);
            }
        }
    }
}

template <int dim>
void SelfCollisionHandler<dim>::computeFrictionEnergy(const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& Vt, const std::vector<MMCVID>& constraintSet,
    const Eigen::VectorXd& multipliers,
    const std::vector<Eigen::Vector2d>& MMDistCoord,
    const std::vector<Eigen::Matrix<double, 3, 2>>& MMTanBasis,
    double& Ef, double eps2, double coef)
{
    double eps = std::sqrt(eps2);

    Eigen::VectorXd EI(constraintSet.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)constraintSet.size(), 1, [&](int cI)
#else
    for (int cI = 0; cI < constraintSet.size(); ++cI)
#endif
        {
            Eigen::RowVector3d relDX3D;

            const auto& MMCVIDI = constraintSet[cI];
            if (MMCVIDI[0] >= 0) {
                // edge-edge
                computeRelDX_EE(V.row(MMCVIDI[0]) - Vt.row(MMCVIDI[0]),
                    V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]),
                    V.row(MMCVIDI[2]) - Vt.row(MMCVIDI[2]),
                    V.row(MMCVIDI[3]) - Vt.row(MMCVIDI[3]),
                    MMDistCoord[cI][0], MMDistCoord[cI][1], relDX3D);
            }
            else {
                // point-triangle and degenerate edge-edge
                assert(MMCVIDI[1] >= 0);
                if (MMCVIDI[2] < 0) {
                    // PP
                    computeRelDX_PP(V.row(-MMCVIDI[0] - 1) - Vt.row(-MMCVIDI[0] - 1),
                        V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]), relDX3D);
                }
                else if (MMCVIDI[3] < 0) {
                    // PE
                    computeRelDX_PE(V.row(-MMCVIDI[0] - 1) - Vt.row(-MMCVIDI[0] - 1),
                        V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]),
                        V.row(MMCVIDI[2]) - Vt.row(MMCVIDI[2]),
                        MMDistCoord[cI][0], relDX3D);
                }
                else {
                    // PT
                    computeRelDX_PT(V.row(-MMCVIDI[0] - 1) - Vt.row(-MMCVIDI[0] - 1),
                        V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]),
                        V.row(MMCVIDI[2]) - Vt.row(MMCVIDI[2]),
                        V.row(MMCVIDI[3]) - Vt.row(MMCVIDI[3]),
                        MMDistCoord[cI][0], MMDistCoord[cI][1], relDX3D);
                }
            }

            EI[cI] = 0.0;
            double relDXSqNorm = (relDX3D * MMTanBasis[cI]).squaredNorm();
            if (relDXSqNorm > eps2) {
                EI[cI] += multipliers[cI] * std::sqrt(relDXSqNorm);
            }
            else {
                double f0;
                f0_SF(relDXSqNorm, eps, f0);
                EI[cI] += multipliers[cI] * f0;
            }
        }
#ifdef USE_TBB
    );
#endif
    Ef = EI.sum() * coef;
}

template <int dim>
void SelfCollisionHandler<dim>::augmentFrictionGradient(const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& Vt, const std::vector<MMCVID>& constraintSet,
    const Eigen::VectorXd& multipliers,
    const std::vector<Eigen::Vector2d>& MMDistCoord,
    const std::vector<Eigen::Matrix<double, 3, 2>>& MMTanBasis,
    Eigen::VectorXd& grad_inc, double eps2, double coef)
{
    double eps = std::sqrt(eps2);

    // TODO: parallelize
    for (int cI = 0; cI < constraintSet.size(); ++cI) {
        Eigen::RowVector3d relDX3D;

        const auto& MMCVIDI = constraintSet[cI];
        if (MMCVIDI[0] >= 0) {
            // edge-edge
            computeRelDX_EE(V.row(MMCVIDI[0]) - Vt.row(MMCVIDI[0]),
                V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]),
                V.row(MMCVIDI[2]) - Vt.row(MMCVIDI[2]),
                V.row(MMCVIDI[3]) - Vt.row(MMCVIDI[3]),
                MMDistCoord[cI][0], MMDistCoord[cI][1], relDX3D);

            Eigen::Vector2d relDX = (relDX3D * MMTanBasis[cI]).transpose();
            double relDXSqNorm = relDX.squaredNorm();
            if (relDXSqNorm > eps2) {
                relDX /= std::sqrt(relDXSqNorm);
            }
            else {
                double f1_div_relDXNorm;
                f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);
                relDX *= f1_div_relDXNorm;
            }

            Eigen::Matrix<double, 12, 1> TTTDX;
            liftRelDXTanToMesh_EE(relDX, MMTanBasis[cI],
                MMDistCoord[cI][0], MMDistCoord[cI][1], TTTDX);
            TTTDX *= coef * multipliers[cI];

            grad_inc.template segment<dim>(MMCVIDI[0] * dim) += TTTDX.template segment<dim>(0);
            grad_inc.template segment<dim>(MMCVIDI[1] * dim) += TTTDX.template segment<dim>(3);
            grad_inc.template segment<dim>(MMCVIDI[2] * dim) += TTTDX.template segment<dim>(6);
            grad_inc.template segment<dim>(MMCVIDI[3] * dim) += TTTDX.template segment<dim>(9);
        }
        else {
            // point-triangle and degenerate edge-edge
            assert(MMCVIDI[1] >= 0);
            if (MMCVIDI[2] < 0) {
                // PP
                computeRelDX_PP(V.row(-MMCVIDI[0] - 1) - Vt.row(-MMCVIDI[0] - 1),
                    V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]), relDX3D);

                Eigen::Vector2d relDX = (relDX3D * MMTanBasis[cI]).transpose();
                double relDXSqNorm = relDX.squaredNorm();
                if (relDXSqNorm > eps2) {
                    relDX /= std::sqrt(relDXSqNorm);
                }
                else {
                    double f1_div_relDXNorm;
                    f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);
                    relDX *= f1_div_relDXNorm;
                }

                Eigen::Matrix<double, 6, 1> TTTDX;
                liftRelDXTanToMesh_PP(relDX, MMTanBasis[cI], TTTDX);
                TTTDX *= coef * multipliers[cI];

                grad_inc.template segment<dim>((-MMCVIDI[0] - 1) * dim) += TTTDX.template segment<dim>(0);
                grad_inc.template segment<dim>(MMCVIDI[1] * dim) += TTTDX.template segment<dim>(3);
            }
            else if (MMCVIDI[3] < 0) {
                // PE
                computeRelDX_PE(V.row(-MMCVIDI[0] - 1) - Vt.row(-MMCVIDI[0] - 1),
                    V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]),
                    V.row(MMCVIDI[2]) - Vt.row(MMCVIDI[2]),
                    MMDistCoord[cI][0], relDX3D);

                Eigen::Vector2d relDX = (relDX3D * MMTanBasis[cI]).transpose();
                double relDXSqNorm = relDX.squaredNorm();
                if (relDXSqNorm > eps2) {
                    relDX /= std::sqrt(relDXSqNorm);
                }
                else {
                    double f1_div_relDXNorm;
                    f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);
                    relDX *= f1_div_relDXNorm;
                }

                Eigen::Matrix<double, 9, 1> TTTDX;
                liftRelDXTanToMesh_PE(relDX, MMTanBasis[cI], MMDistCoord[cI][0], TTTDX);
                TTTDX *= coef * multipliers[cI];

                grad_inc.template segment<dim>((-MMCVIDI[0] - 1) * dim) += TTTDX.template segment<dim>(0);
                grad_inc.template segment<dim>(MMCVIDI[1] * dim) += TTTDX.template segment<dim>(3);
                grad_inc.template segment<dim>(MMCVIDI[2] * dim) += TTTDX.template segment<dim>(6);
            }
            else {
                // PT
                computeRelDX_PT(V.row(-MMCVIDI[0] - 1) - Vt.row(-MMCVIDI[0] - 1),
                    V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]),
                    V.row(MMCVIDI[2]) - Vt.row(MMCVIDI[2]),
                    V.row(MMCVIDI[3]) - Vt.row(MMCVIDI[3]),
                    MMDistCoord[cI][0], MMDistCoord[cI][1], relDX3D);

                Eigen::Vector2d relDX = (relDX3D * MMTanBasis[cI]).transpose();
                double relDXSqNorm = relDX.squaredNorm();
                if (relDXSqNorm > eps2) {
                    relDX /= std::sqrt(relDXSqNorm);
                }
                else {
                    double f1_div_relDXNorm;
                    f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);
                    relDX *= f1_div_relDXNorm;
                }

                Eigen::Matrix<double, 12, 1> TTTDX;
                liftRelDXTanToMesh_PT(relDX, MMTanBasis[cI],
                    MMDistCoord[cI][0], MMDistCoord[cI][1], TTTDX);
                TTTDX *= coef * multipliers[cI];

                grad_inc.template segment<dim>((-MMCVIDI[0] - 1) * dim) += TTTDX.template segment<dim>(0);
                grad_inc.template segment<dim>(MMCVIDI[1] * dim) += TTTDX.template segment<dim>(3);
                grad_inc.template segment<dim>(MMCVIDI[2] * dim) += TTTDX.template segment<dim>(6);
                grad_inc.template segment<dim>(MMCVIDI[3] * dim) += TTTDX.template segment<dim>(9);
            }
        }
    }
}

template <int dim>
void SelfCollisionHandler<dim>::augmentFrictionHessian(const Mesh<dim>& mesh,
    const Eigen::MatrixXd& Vt, const std::vector<MMCVID>& constraintSet,
    const Eigen::VectorXd& multipliers,
    const std::vector<Eigen::Vector2d>& MMDistCoord,
    const std::vector<Eigen::Matrix<double, 3, 2>>& MMTanBasis,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* mtr_incremental,
    double eps2, double coef, bool projectDBC)
{
    augmentFrictionHessian(
        mesh, Vt, constraintSet, multipliers, MMDistCoord, MMTanBasis,
        [&mtr_incremental](size_t row, size_t col, double val) {
            mtr_incremental->addCoeff(row, col, val);
        },
        eps2, coef, projectDBC);
}

template <int dim>
void SelfCollisionHandler<dim>::augmentFrictionHessian(const Mesh<dim>& mesh,
    const Eigen::MatrixXd& Vt, const std::vector<MMCVID>& constraintSet,
    const Eigen::VectorXd& multipliers,
    const std::vector<Eigen::Vector2d>& MMDistCoord,
    const std::vector<Eigen::Matrix<double, 3, 2>>& MMTanBasis,
    const std::function<void(size_t, size_t, double)>& addCoeff,
    double eps2, double coef, bool projectDBC)
{
    double eps = std::sqrt(eps2);

    std::vector<Eigen::Matrix<double, 12, 12>> IPHessian(constraintSet.size());
    std::vector<Eigen::Matrix<int, 4, 1>> rowIStart(constraintSet.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)constraintSet.size(), 1, [&](int cI)
#else
    for (int cI = 0; cI < constraintSet.size(); ++cI)
#endif
        {
            const auto& MMCVIDI = constraintSet[cI];
            if (MMCVIDI[0] >= 0) {
                // edge-edge
                Eigen::RowVector3d relDX3D;
                computeRelDX_EE(mesh.V.row(MMCVIDI[0]) - Vt.row(MMCVIDI[0]),
                    mesh.V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]),
                    mesh.V.row(MMCVIDI[2]) - Vt.row(MMCVIDI[2]),
                    mesh.V.row(MMCVIDI[3]) - Vt.row(MMCVIDI[3]),
                    MMDistCoord[cI][0], MMDistCoord[cI][1], relDX3D);

                Eigen::Vector2d relDX = (relDX3D * MMTanBasis[cI]).transpose();
                double relDXSqNorm = relDX.squaredNorm();
                double relDXNorm = std::sqrt(relDXSqNorm);

                computeTTT_EE(MMTanBasis[cI], MMDistCoord[cI][0], MMDistCoord[cI][1], IPHessian[cI]);
                if (relDXSqNorm > eps2) {
                    IPHessian[cI] *= coef * multipliers[cI] / relDXNorm;

                    Eigen::Matrix<double, 12, 1> TTTDX;
                    liftRelDXTanToMesh_EE(relDX, MMTanBasis[cI],
                        MMDistCoord[cI][0], MMDistCoord[cI][1], TTTDX);
                    IPHessian[cI] -= (TTTDX * (coef * multipliers[cI] / (relDXSqNorm * relDXNorm))) * TTTDX.transpose();

                    IglUtils::makePD(IPHessian[cI]);
                }
                else {
                    double f1_div_relDXNorm;
                    f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);

                    IPHessian[cI] *= coef * multipliers[cI] * f1_div_relDXNorm;

                    double f2;
                    f2_SF(relDXSqNorm, eps, f2);
                    if (f2 != f1_div_relDXNorm && relDXSqNorm) {
                        // higher order clamping and not exactly static
                        Eigen::Matrix<double, 12, 1> TTTDX;
                        liftRelDXTanToMesh_EE(relDX, MMTanBasis[cI],
                            MMDistCoord[cI][0], MMDistCoord[cI][1], TTTDX);
                        IPHessian[cI] += (TTTDX * (coef * multipliers[cI] * (f2 - f1_div_relDXNorm) / relDXSqNorm)) * TTTDX.transpose();

                        IglUtils::makePD(IPHessian[cI]);
                    }
                }

                rowIStart[cI][0] = mesh.isProjectDBCVertex(MMCVIDI[0], projectDBC) ? -1 : (MMCVIDI[0] * dim);
                rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                rowIStart[cI][2] = mesh.isProjectDBCVertex(MMCVIDI[2], projectDBC) ? -1 : (MMCVIDI[2] * dim);
                rowIStart[cI][3] = mesh.isProjectDBCVertex(MMCVIDI[3], projectDBC) ? -1 : (MMCVIDI[3] * dim);
            }
            else {
                // point-triangle and degenerate edge-edge
                assert(MMCVIDI[1] >= 0);
                int v0I = -MMCVIDI[0] - 1;
                if (MMCVIDI[2] < 0) {
                    // PP
                    Eigen::RowVector3d relDX3D;
                    computeRelDX_PP(mesh.V.row(v0I) - Vt.row(v0I),
                        mesh.V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]), relDX3D);

                    Eigen::Vector2d relDX = (relDX3D * MMTanBasis[cI]).transpose();
                    double relDXSqNorm = relDX.squaredNorm();
                    double relDXNorm = std::sqrt(relDXSqNorm);

                    Eigen::Matrix<double, 6, 6> HessianBlock;
                    computeTTT_PP(MMTanBasis[cI], HessianBlock);
                    if (relDXSqNorm > eps2) {
                        HessianBlock *= coef * multipliers[cI] / relDXNorm;

                        Eigen::Matrix<double, 6, 1> TTTDX;
                        liftRelDXTanToMesh_PP(relDX, MMTanBasis[cI], TTTDX);
                        HessianBlock -= (TTTDX * (coef * multipliers[cI] / (relDXSqNorm * relDXNorm))) * TTTDX.transpose();

                        IglUtils::makePD(HessianBlock);
                    }
                    else {
                        double f1_div_relDXNorm;
                        f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);

                        HessianBlock *= coef * multipliers[cI] * f1_div_relDXNorm;

                        double f2;
                        f2_SF(relDXSqNorm, eps, f2);
                        if (f2 != f1_div_relDXNorm && relDXSqNorm) {
                            // higher order clamping and not exactly static
                            Eigen::Matrix<double, 6, 1> TTTDX;
                            liftRelDXTanToMesh_PP(relDX, MMTanBasis[cI], TTTDX);
                            HessianBlock += (TTTDX * (coef * multipliers[cI] * (f2 - f1_div_relDXNorm) / relDXSqNorm)) * TTTDX.transpose();

                            IglUtils::makePD(HessianBlock);
                        }
                    }
                    IPHessian[cI].template block<6, 6>(0, 0) = HessianBlock;

                    rowIStart[cI][0] = mesh.isProjectDBCVertex(v0I, projectDBC) ? -1 : (v0I * dim);
                    rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                    rowIStart[cI][2] = -1;
                    rowIStart[cI][3] = -1;
                }
                else if (MMCVIDI[3] < 0) {
                    // PE
                    Eigen::RowVector3d relDX3D;
                    computeRelDX_PE(mesh.V.row(v0I) - Vt.row(v0I),
                        mesh.V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]),
                        mesh.V.row(MMCVIDI[2]) - Vt.row(MMCVIDI[2]),
                        MMDistCoord[cI][0], relDX3D);

                    Eigen::Vector2d relDX = (relDX3D * MMTanBasis[cI]).transpose();
                    double relDXSqNorm = relDX.squaredNorm();
                    double relDXNorm = std::sqrt(relDXSqNorm);

                    Eigen::Matrix<double, 9, 9> HessianBlock;
                    computeTTT_PE(MMTanBasis[cI], MMDistCoord[cI][0], HessianBlock);
                    if (relDXSqNorm > eps2) {
                        HessianBlock *= coef * multipliers[cI] / relDXNorm;

                        Eigen::Matrix<double, 9, 1> TTTDX;
                        liftRelDXTanToMesh_PE(relDX, MMTanBasis[cI], MMDistCoord[cI][0], TTTDX);
                        HessianBlock -= (TTTDX * (coef * multipliers[cI] / (relDXSqNorm * relDXNorm))) * TTTDX.transpose();

                        IglUtils::makePD(HessianBlock);
                    }
                    else {
                        double f1_div_relDXNorm;
                        f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);

                        HessianBlock *= coef * multipliers[cI] * f1_div_relDXNorm;

                        double f2;
                        f2_SF(relDXSqNorm, eps, f2);
                        if (f2 != f1_div_relDXNorm && relDXSqNorm) {
                            // higher order clamping and not exactly static
                            Eigen::Matrix<double, 9, 1> TTTDX;
                            liftRelDXTanToMesh_PE(relDX, MMTanBasis[cI], MMDistCoord[cI][0], TTTDX);
                            HessianBlock += (TTTDX * (coef * multipliers[cI] * (f2 - f1_div_relDXNorm) / relDXSqNorm)) * TTTDX.transpose();

                            IglUtils::makePD(HessianBlock);
                        }
                    }
                    IPHessian[cI].template block<9, 9>(0, 0) = HessianBlock;

                    rowIStart[cI][0] = mesh.isProjectDBCVertex(v0I, projectDBC) ? -1 : (v0I * dim);
                    rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                    rowIStart[cI][2] = mesh.isProjectDBCVertex(MMCVIDI[2], projectDBC) ? -1 : (MMCVIDI[2] * dim);
                    rowIStart[cI][3] = -1;
                }
                else {
                    // PT
                    Eigen::RowVector3d relDX3D;
                    computeRelDX_PT(mesh.V.row(v0I) - Vt.row(v0I),
                        mesh.V.row(MMCVIDI[1]) - Vt.row(MMCVIDI[1]),
                        mesh.V.row(MMCVIDI[2]) - Vt.row(MMCVIDI[2]),
                        mesh.V.row(MMCVIDI[3]) - Vt.row(MMCVIDI[3]),
                        MMDistCoord[cI][0], MMDistCoord[cI][1], relDX3D);

                    Eigen::Vector2d relDX = (relDX3D * MMTanBasis[cI]).transpose();
                    double relDXSqNorm = relDX.squaredNorm();
                    double relDXNorm = std::sqrt(relDXSqNorm);

                    computeTTT_PT(MMTanBasis[cI], MMDistCoord[cI][0], MMDistCoord[cI][1], IPHessian[cI]);
                    if (relDXSqNorm > eps2) {
                        IPHessian[cI] *= coef * multipliers[cI] / relDXNorm;
                        Eigen::Matrix<double, 12, 1> TTTDX;
                        liftRelDXTanToMesh_PT(relDX, MMTanBasis[cI],
                            MMDistCoord[cI][0], MMDistCoord[cI][1], TTTDX);
                        IPHessian[cI] -= (TTTDX * (coef * multipliers[cI] / (relDXSqNorm * relDXNorm))) * TTTDX.transpose();

                        IglUtils::makePD(IPHessian[cI]);
                    }
                    else {
                        double f1_div_relDXNorm;
                        f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);

                        IPHessian[cI] *= coef * multipliers[cI] * f1_div_relDXNorm;

                        double f2;
                        f2_SF(relDXSqNorm, eps, f2);
                        if (f2 != f1_div_relDXNorm && relDXSqNorm) {
                            // higher order clamping and not exactly static
                            Eigen::Matrix<double, 12, 1> TTTDX;
                            liftRelDXTanToMesh_PT(relDX, MMTanBasis[cI],
                                MMDistCoord[cI][0], MMDistCoord[cI][1], TTTDX);
                            IPHessian[cI] += (TTTDX * (coef * multipliers[cI] * (f2 - f1_div_relDXNorm) / relDXSqNorm)) * TTTDX.transpose();

                            IglUtils::makePD(IPHessian[cI]);
                        }
                    }

                    rowIStart[cI][0] = mesh.isProjectDBCVertex(v0I, projectDBC) ? -1 : (v0I * dim);
                    rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                    rowIStart[cI][2] = mesh.isProjectDBCVertex(MMCVIDI[2], projectDBC) ? -1 : (MMCVIDI[2] * dim);
                    rowIStart[cI][3] = mesh.isProjectDBCVertex(MMCVIDI[3], projectDBC) ? -1 : (MMCVIDI[3] * dim);
                }
            }
        }
#ifdef USE_TBB
    );
#endif

    // TODO: parallelize
    for (int cI = 0; cI < constraintSet.size(); ++cI) {
        for (int i = 0; i < rowIStart[cI].size(); ++i) {
            int rowIStartI = rowIStart[cI][i];
            if (rowIStartI >= 0) {
                for (int j = 0; j < rowIStart[cI].size(); ++j) {
                    int colIStartI = rowIStart[cI][j];
                    if (colIStartI >= 0) {
                        addCoeff(rowIStartI, colIStartI, IPHessian[cI](i * dim, j * dim));
                        addCoeff(rowIStartI, colIStartI + 1, IPHessian[cI](i * dim, j * dim + 1));
                        addCoeff(rowIStartI + 1, colIStartI, IPHessian[cI](i * dim + 1, j * dim));
                        addCoeff(rowIStartI + 1, colIStartI + 1, IPHessian[cI](i * dim + 1, j * dim + 1));
                        if constexpr (dim == 3) {
                            addCoeff(rowIStartI, colIStartI + 2, IPHessian[cI](i * dim, j * dim + 2));
                            addCoeff(rowIStartI + 1, colIStartI + 2, IPHessian[cI](i * dim + 1, j * dim + 2));

                            addCoeff(rowIStartI + 2, colIStartI, IPHessian[cI](i * dim + 2, j * dim));
                            addCoeff(rowIStartI + 2, colIStartI + 1, IPHessian[cI](i * dim + 2, j * dim + 1));
                            addCoeff(rowIStartI + 2, colIStartI + 2, IPHessian[cI](i * dim + 2, j * dim + 2));
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void SelfCollisionHandler<dim>::augmentParaEEGradient(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& paraEEMMCVIDSet,
    const std::vector<std::pair<int, int>>& paraEEeIeJSet,
    Eigen::VectorXd& grad_inc, double dHat, double coef)
{
    Eigen::VectorXd e_db_div_dd;
    SelfCollisionHandler<dim>::evaluateConstraints(mesh, paraEEMMCVIDSet, e_db_div_dd);
    for (int cI = 0; cI < e_db_div_dd.size(); ++cI) {
        double b;
        compute_b(e_db_div_dd[cI], dHat, b);
        compute_g_b(e_db_div_dd[cI], dHat, e_db_div_dd[cI]);

        const MMCVID& MMCVIDI = paraEEMMCVIDSet[cI];
        double eps_x, e;
        double coef_b = coef * b;
        if (MMCVIDI[3] >= 0) {
            // EE
            compute_eps_x(mesh, MMCVIDI[0], MMCVIDI[1], MMCVIDI[2], MMCVIDI[3], eps_x);
            compute_e(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), eps_x, e);

            Eigen::Matrix<double, 12, 1> e_g;
            compute_e_g(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), eps_x, e_g);
            grad_inc.template segment<dim>(MMCVIDI[0] * dim) += coef_b * e_g.template segment<dim>(0);
            grad_inc.template segment<dim>(MMCVIDI[1] * dim) += coef_b * e_g.template segment<dim>(dim);
            grad_inc.template segment<dim>(MMCVIDI[2] * dim) += coef_b * e_g.template segment<dim>(dim * 2);
            grad_inc.template segment<dim>(MMCVIDI[3] * dim) += coef_b * e_g.template segment<dim>(dim * 3);
        }
        else {
            // PP or PE
            const std::pair<int, int>& eIeJ = paraEEeIeJSet[cI];
            const std::pair<int, int>& eI = mesh.SFEdges[eIeJ.first];
            const std::pair<int, int>& eJ = mesh.SFEdges[eIeJ.second];
            compute_eps_x(mesh, eI.first, eI.second, eJ.first, eJ.second, eps_x);
            compute_e(mesh.V.row(eI.first), mesh.V.row(eI.second), mesh.V.row(eJ.first), mesh.V.row(eJ.second), eps_x, e);

            Eigen::Matrix<double, 12, 1> e_g;
            compute_e_g(mesh.V.row(eI.first), mesh.V.row(eI.second), mesh.V.row(eJ.first), mesh.V.row(eJ.second), eps_x, e_g);
            grad_inc.template segment<dim>(eI.first * dim) += coef_b * e_g.template segment<dim>(0);
            grad_inc.template segment<dim>(eI.second * dim) += coef_b * e_g.template segment<dim>(dim);
            grad_inc.template segment<dim>(eJ.first * dim) += coef_b * e_g.template segment<dim>(dim * 2);
            grad_inc.template segment<dim>(eJ.second * dim) += coef_b * e_g.template segment<dim>(dim * 3);
        }
        e_db_div_dd[cI] *= e;
    }
    SelfCollisionHandler<dim>::leftMultiplyConstraintJacobianT(mesh, paraEEMMCVIDSet,
        e_db_div_dd, grad_inc, coef);
}

template <int dim>
void SelfCollisionHandler<dim>::augmentParaEEHessian(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& paraEEMMCVIDSet,
    const std::vector<std::pair<int, int>>& paraEEeIeJSet,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* H_inc,
    double dHat, double coef, bool projectDBC)
{
    if constexpr (dim == 3) {
        std::vector<Eigen::Matrix<double, 4 * dim, 4 * dim>> PEEHessian(paraEEMMCVIDSet.size());
        std::vector<Eigen::Matrix<int, 4, 1>> rowIStart(paraEEMMCVIDSet.size());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)paraEEMMCVIDSet.size(), 1, [&](int cI)
#else
        for (int cI = 0; cI < paraEEMMCVIDSet.size(); ++cI)
#endif
            {
                double b, g_b, H_b;

                double eps_x, e;
                Eigen::Matrix<double, 12, 1> e_g;
                Eigen::Matrix<double, 12, 12> e_H;

                const MMCVID& MMCVIDI = paraEEMMCVIDSet[cI];
                double d;
                evaluateConstraint(mesh, MMCVIDI, d);
                Eigen::Matrix<double, 12, 1> grad_d;
                Eigen::Matrix<double, 12, 12> H_d;

                if (MMCVIDI[3] >= 0) {
                    // EE
                    compute_b(d, dHat, b);
                    compute_g_b(d, dHat, g_b);
                    compute_H_b(d, dHat, H_b);

                    compute_eps_x(mesh, MMCVIDI[0], MMCVIDI[1], MMCVIDI[2], MMCVIDI[3], eps_x);
                    compute_e(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), eps_x, e);
                    compute_e_g(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), eps_x, e_g);
                    compute_e_H(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), eps_x, e_H);

                    g_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), grad_d);
                    H_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), mesh.V.row(MMCVIDI[3]), H_d);

                    rowIStart[cI][0] = mesh.isProjectDBCVertex(MMCVIDI[0], projectDBC) ? -1 : (MMCVIDI[0] * dim);
                    rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                    rowIStart[cI][2] = mesh.isProjectDBCVertex(MMCVIDI[2], projectDBC) ? -1 : (MMCVIDI[2] * dim);
                    rowIStart[cI][3] = mesh.isProjectDBCVertex(MMCVIDI[3], projectDBC) ? -1 : (MMCVIDI[3] * dim);
                }
                else {
                    // PP or PE
                    compute_b(d, dHat, b);
                    compute_g_b(d, dHat, g_b);
                    compute_H_b(d, dHat, H_b);

                    const std::pair<int, int>& eIeJ = paraEEeIeJSet[cI];
                    const std::pair<int, int>& eI = mesh.SFEdges[eIeJ.first];
                    const std::pair<int, int>& eJ = mesh.SFEdges[eIeJ.second];
                    compute_eps_x(mesh, eI.first, eI.second, eJ.first, eJ.second, eps_x);
                    compute_e(mesh.V.row(eI.first), mesh.V.row(eI.second), mesh.V.row(eJ.first), mesh.V.row(eJ.second), eps_x, e);
                    compute_e_g(mesh.V.row(eI.first), mesh.V.row(eI.second), mesh.V.row(eJ.first), mesh.V.row(eJ.second), eps_x, e_g);
                    compute_e_H(mesh.V.row(eI.first), mesh.V.row(eI.second), mesh.V.row(eJ.first), mesh.V.row(eJ.second), eps_x, e_H);

                    rowIStart[cI][0] = mesh.isProjectDBCVertex(eI.first, projectDBC) ? -1 : (eI.first * dim);
                    rowIStart[cI][1] = mesh.isProjectDBCVertex(eI.second, projectDBC) ? -1 : (eI.second * dim);
                    rowIStart[cI][2] = mesh.isProjectDBCVertex(eJ.first, projectDBC) ? -1 : (eJ.first * dim);
                    rowIStart[cI][3] = mesh.isProjectDBCVertex(eJ.second, projectDBC) ? -1 : (eJ.second * dim);

                    int v0I = -MMCVIDI[0] - 1;
                    if (MMCVIDI[2] >= 0) {
                        // PE
                        Eigen::Matrix<double, 9, 1> grad_d_PE;
                        g_PE(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), grad_d_PE);
                        Eigen::Matrix<double, 9, 9> H_d_PE;
                        H_PE(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), mesh.V.row(MMCVIDI[2]), H_d_PE);

                        // fill in
                        int ind[4] = { eI.first, eI.second, eJ.first, eJ.second };
                        int indMap[3];
                        for (int i = 0; i < 4; ++i) {
                            if (v0I == ind[i]) {
                                indMap[0] = i;
                            }
                            else if (MMCVIDI[1] == ind[i]) {
                                indMap[1] = i;
                            }
                            else if (MMCVIDI[2] == ind[i]) {
                                indMap[2] = i;
                            }
                        }

                        grad_d.setZero();
                        H_d.setZero();
                        for (int i = 0; i < 3; ++i) {
                            grad_d.template segment<dim>(indMap[i] * dim) = grad_d_PE.template segment<dim>(i * dim);
                            for (int j = 0; j < 3; ++j) {
                                H_d.template block<dim, dim>(indMap[i] * dim, indMap[j] * dim) = H_d_PE.template block<dim, dim>(i * dim, j * dim);
                            }
                        }
                    }
                    else {
                        // PP
                        Eigen::Matrix<double, 6, 1> grad_d_PP;
                        g_PP(mesh.V.row(v0I), mesh.V.row(MMCVIDI[1]), grad_d_PP);
                        Eigen::Matrix<double, 6, 6> H_d_PP;
                        H_PP(H_d_PP);

                        int ind[4] = { eI.first, eI.second, eJ.first, eJ.second };
                        int indMap[2];
                        for (int i = 0; i < 4; ++i) {
                            if (v0I == ind[i]) {
                                indMap[0] = i;
                            }
                            else if (MMCVIDI[1] == ind[i]) {
                                indMap[1] = i;
                            }
                        }

                        grad_d.setZero();
                        H_d.setZero();
                        for (int i = 0; i < 2; ++i) {
                            grad_d.template segment<dim>(indMap[i] * dim) = grad_d_PP.template segment<dim>(i * dim);
                            for (int j = 0; j < 2; ++j) {
                                H_d.template block<dim, dim>(indMap[i] * dim, indMap[j] * dim) = H_d_PP.template block<dim, dim>(i * dim, j * dim);
                            }
                        }
                    }
                }

                Eigen::Matrix<double, 12, 12> kappa_gradb_gradeT;
                kappa_gradb_gradeT = ((coef * g_b) * grad_d) * e_g.transpose();

                PEEHessian[cI] = kappa_gradb_gradeT + kappa_gradb_gradeT.transpose() + (coef * b) * e_H + ((coef * e * H_b) * grad_d) * grad_d.transpose() + (coef * e * g_b) * H_d;
                IglUtils::makePD(PEEHessian[cI]);
            }
#ifdef USE_TBB
        );
#endif

        // TODO: parallelize
        for (int cI = 0; cI < paraEEMMCVIDSet.size(); ++cI) {
            for (int i = 0; i < rowIStart[cI].size(); ++i) {
                int rowIStartI = rowIStart[cI][i];
                if (rowIStartI >= 0) {
                    for (int j = 0; j < rowIStart[cI].size(); ++j) {
                        int colIStartI = rowIStart[cI][j];
                        if (colIStartI >= 0) {
                            H_inc->addCoeff(rowIStartI, colIStartI, PEEHessian[cI](i * dim, j * dim));
                            H_inc->addCoeff(rowIStartI, colIStartI + 1, PEEHessian[cI](i * dim, j * dim + 1));
                            H_inc->addCoeff(rowIStartI + 1, colIStartI, PEEHessian[cI](i * dim + 1, j * dim));
                            H_inc->addCoeff(rowIStartI + 1, colIStartI + 1, PEEHessian[cI](i * dim + 1, j * dim + 1));
                            if constexpr (dim == 3) {
                                H_inc->addCoeff(rowIStartI, colIStartI + 2, PEEHessian[cI](i * dim, j * dim + 2));
                                H_inc->addCoeff(rowIStartI + 1, colIStartI + 2, PEEHessian[cI](i * dim + 1, j * dim + 2));

                                H_inc->addCoeff(rowIStartI + 2, colIStartI, PEEHessian[cI](i * dim + 2, j * dim));
                                H_inc->addCoeff(rowIStartI + 2, colIStartI + 1, PEEHessian[cI](i * dim + 2, j * dim + 1));
                                H_inc->addCoeff(rowIStartI + 2, colIStartI + 2, PEEHessian[cI](i * dim + 2, j * dim + 2));
                            }
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
bool SelfCollisionHandler<dim>::checkEdgeTriIntersection(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh)
{
    std::vector<std::vector<int>> intersectedEdges(mesh.SF.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SF.rows(), 1, [&](int sfI)
#else
    for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
        {
            const Eigen::RowVector3i& sfVInd = mesh.SF.row(sfI);
#ifdef USE_SH_CCS
            std::unordered_set<int> sEdgeInds;
            sh.queryTriangleForEdges(mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), 0.0, sEdgeInds);
            for (const auto& eI : sEdgeInds)
#else
        for (int eI = 0; eI < mesh.SFEdges.size(); ++eI)
#endif
            {
                const auto& meshEI = mesh.SFEdges[eI];
                if (meshEI.first == sfVInd[0] || meshEI.first == sfVInd[1] || meshEI.first == sfVInd[2] || meshEI.second == sfVInd[0] || meshEI.second == sfVInd[1] || meshEI.second == sfVInd[2]) {
                    continue;
                }
                if (IglUtils::segTriIntersect(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                        mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]))) {
                    intersectedEdges[sfI].emplace_back(eI);
                }
            }
        }
#ifdef USE_TBB
    );
#endif

    bool intersected = false;
    for (int sfI = 0; sfI < intersectedEdges.size(); ++sfI) {
        for (const auto& eI : intersectedEdges[sfI]) {
            intersected = true;
            std::cout << "self edge - triangle intersection detected" << std::endl;
            std::cout << mesh.SFEdges[eI].first << " " << mesh.SFEdges[eI].second << std::endl;
            std::cout << mesh.SF.row(sfI) << std::endl;
            std::cout << mesh.V.row(mesh.SFEdges[eI].first) << std::endl;
            std::cout << mesh.V.row(mesh.SFEdges[eI].second) << std::endl;
            std::cout << mesh.V.row(mesh.SF(sfI, 0)) << std::endl;
            std::cout << mesh.V.row(mesh.SF(sfI, 1)) << std::endl;
            std::cout << mesh.V.row(mesh.SF(sfI, 2)) << std::endl;
        }
    }
    return !intersected;
}

template <int dim>
bool SelfCollisionHandler<dim>::checkEdgeTriIntersectionIfAny(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh)
{
    Eigen::ArrayXi intersected(mesh.SF.rows());
    intersected.setZero();
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SF.rows(), 1, [&](int sfI)
#else
    for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
        {
            const Eigen::RowVector3i& sfVInd = mesh.SF.row(sfI);
            int coDim_sfI = mesh.sfICoDim(sfI);
#ifdef USE_SH_CCS
            std::unordered_set<int> sEdgeInds;
            sh.queryTriangleForEdges(mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), 0.0, sEdgeInds);
            for (const auto& eI : sEdgeInds)
#else
        for (int eI = 0; eI < mesh.SFEdges.size(); ++eI)
#endif
            {
                const auto& meshEI = mesh.SFEdges[eI];
                if (meshEI.first == sfVInd[0] || meshEI.first == sfVInd[1] || meshEI.first == sfVInd[2] || meshEI.second == sfVInd[0] || meshEI.second == sfVInd[1] || meshEI.second == sfVInd[2]) {
                    continue;
                }

                int coDim_eI = mesh.vICoDim(meshEI.first);
                if ((coDim_sfI < 3 && coDim_eI < 3) || (mesh.isDBCVertex(meshEI.first) && mesh.isDBCVertex(meshEI.second) && mesh.isDBCVertex(sfVInd[0]) && mesh.isDBCVertex(sfVInd[1]) && mesh.isDBCVertex(sfVInd[2]))) {
                    continue;
                }

                if (IglUtils::segTriIntersect(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                        mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]))) {
                    intersected[sfI] = 1;
                    break;
                }
            }
        }
#ifdef USE_TBB
    );
#endif

    if ((intersected != 0).any()) {
        return false;
    }

    for (int compI = 0; compI < mesh.componentCoDim.size(); ++compI) {
        if (mesh.componentCoDim[compI] == 0) {
            Eigen::ArrayXi intersected_pTet(mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI]);
            intersected_pTet.setZero();
#ifdef USE_TBB
            tbb::parallel_for(mesh.componentNodeRange[compI], mesh.componentNodeRange[compI + 1], 1, [&](int vI)
#else
            for (int vI = mesh.componentNodeRange[compI]; vI < mesh.componentNodeRange[compI + 1]; ++vI)
#endif
                {
                    for (int fI = 0; fI < mesh.F.rows(); ++fI) {
                        const Eigen::Matrix<int, 1, dim + 1>& fVInd = mesh.F.row(fI);
                        Eigen::Array<double, 1, dim> mins = mesh.V.row(fVInd[0]).array().min(
                                                                                            mesh.V.row(fVInd[1]).array())
                                                                .min(mesh.V.row(fVInd[2]).array())
                                                                .min(mesh.V.row(fVInd[3]).array());
                        Eigen::Array<double, 1, dim> maxs = mesh.V.row(fVInd[0]).array().max(
                                                                                            mesh.V.row(fVInd[1]).array())
                                                                .max(mesh.V.row(fVInd[2]).array())
                                                                .max(mesh.V.row(fVInd[3]).array());
                        if ((mins <= mesh.V.row(vI).array()).all() && (maxs >= mesh.V.row(vI).array()).all()) {
                            if (IglUtils::pointInsideTetrahedron(mesh.V.row(vI),
                                    mesh.V.row(fVInd[0]), mesh.V.row(fVInd[1]),
                                    mesh.V.row(fVInd[2]), mesh.V.row(fVInd[3]))) {
                                intersected_pTet[vI - mesh.componentNodeRange[compI]] = 1;
                            }
                        }
                    }
                }
#ifdef USE_TBB
            );
#endif
            if ((intersected_pTet != 0).any()) {
                return false;
            }
        }
    }

    return true;
}

// Determine if the mesh is in an intersected state.
template <int dim>
bool SelfCollisionHandler<dim>::isIntersected(
    const Mesh<dim>& mesh, const Eigen::MatrixXd& V0,
    const ccd::CCDMethod method)
{
#ifdef USE_EXACT_CCD
    if (method == ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER) {
        return false;
    }

#ifdef USE_SH_INTERSECTED
    MatrixXd searchDir = mesh.V - V0;
    searchDir = Map<VectorXd>(searchDir.data(), searchDir.size());
    SpatialHash<dim> sh(mesh, -searchDir, 1, mesh.avgEdgeLen / 3.0);
#endif

    // point-triangle
    // Loop over mesh surface vertices
    for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
        int vI = mesh.SVI[svI];
        // Loop over mesh faces
#ifdef USE_SH_INTERSECTED
        // sVInds and sEdgeInds are not used as we only check point-triangles
        std::unordered_set<int> sVInds, sEdgeInds, sTriInds;
        sh.queryPointForPrimitives(svI, sVInds, sEdgeInds, sTriInds);
        for (const auto& sfI : sTriInds) {
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI) {
#endif
            const RowVector3i& sfVInd = mesh.SF.row(sfI);
            if (vI == sfVInd[0] || vI == sfVInd[1] || vI == sfVInd[2]) {
                continue; // Skip triangles that contain the point
            }
            if (ccd::vertexFaceCCD(
                    V0.row(vI).transpose(),
                    V0.row(sfVInd[0]).transpose(),
                    V0.row(sfVInd[1]).transpose(),
                    V0.row(sfVInd[2]).transpose(),
                    mesh.V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose(),
                    mesh.V.row(sfVInd[1]).transpose(),
                    mesh.V.row(sfVInd[2]).transpose(),
                    method)) {
                return true;
            }
        }
    }

    // edge-edge
    // Loop over FEM mesh edges
    // Get the FEM mesh edge as a pair of indiced into mesh.V
    for (int i = 0; i < mesh.SFEdges.size(); i++) {
        const std::pair<int, int>& edge1 = mesh.SFEdges[i];
        // Loop over meshCO edges
        // Get the mesh CO edge as a pair of indiced into Base::V
#ifdef USE_SH_INTERSECTED
        std::unordered_set<int> sEdgeInds;
        sh.queryEdgeForEdges(i, sEdgeInds);
        for (const int& j : sEdgeInds) {
#else
        for (int j = i + 1; j < mesh.SFEdges.size(); j++) {
#endif
            const std::pair<int, int>& edge2 = mesh.SFEdges[j];
            if (edge1.first == edge2.first || edge1.first == edge2.second
                || edge1.second == edge2.first || edge1.second == edge2.second) {
                continue;
            }
            if (ccd::edgeEdgeCCD(
                    V0.row(edge1.first).transpose(),
                    V0.row(edge1.second).transpose(),
                    V0.row(edge2.first).transpose(),
                    V0.row(edge2.second).transpose(),
                    mesh.V.row(edge1.first).transpose(),
                    mesh.V.row(edge1.second).transpose(),
                    mesh.V.row(edge2.first).transpose(),
                    mesh.V.row(edge2.second).transpose(),
                    method)) {
                return true;
            }
        }
    }
#endif
    return false;
}

template class SelfCollisionHandler<DIM>;

} // namespace IPC
