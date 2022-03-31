//
//  MeshCO.cpp
//  IPC
//
//  Created by Minchen Li on 4/25/19.
//

#include "MeshCO.hpp"

#include <stdexcept>

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

#include <igl/readOBJ.h>
#include <spdlog/spdlog.h>

#ifdef USE_TBB
#include <tbb/parallel_for.h>
#include <mutex>
#endif

extern Timer timer_temp3, timer_mt;

namespace IPC {

template <int dim>
MeshCO<dim>::MeshCO(const char* objFilePath,
    const Eigen::Matrix<double, dim, 1>& p_position,
    const Eigen::Matrix3d& p_rotation,
    double scale, double p_friction)
{
    this->objFilePath = objFilePath;

    Base::origin = p_position;
    this->scale = scale;

    Base::friction = p_friction;

    igl::readOBJ(objFilePath, Base::V, Base::F);

    Eigen::Matrix<double, dim, 1> curPos = (Base::V.colwise().sum() / Base::V.rows()).transpose();
    double curScale = (Base::V.colwise().maxCoeff() - Base::V.colwise().minCoeff()).maxCoeff();
    Base::V.rowwise() -= curPos.transpose();
    for (int i = 0; i < Base::V.rows(); ++i) {
        Base::V.row(i) *= p_rotation.transpose();
    }
    Base::V *= scale / curScale;
    Base::V.rowwise() += Base::origin.transpose();

    if constexpr (dim == 3) {
        std::set<std::pair<int, int>> edgesSet;
        for (int sfI = 0; sfI < Base::F.rows(); ++sfI) {
            auto finder = edgesSet.find(std::pair<int, int>(Base::F(sfI, 1), Base::F(sfI, 0)));
            if (finder == edgesSet.end()) {
                edgesSet.insert(std::pair<int, int>(Base::F(sfI, 0), Base::F(sfI, 1)));
            }

            finder = edgesSet.find(std::pair<int, int>(Base::F(sfI, 2), Base::F(sfI, 1)));
            if (finder == edgesSet.end()) {
                edgesSet.insert(std::pair<int, int>(Base::F(sfI, 1), Base::F(sfI, 2)));
            }

            finder = edgesSet.find(std::pair<int, int>(Base::F(sfI, 0), Base::F(sfI, 2)));
            if (finder == edgesSet.end()) {
                edgesSet.insert(std::pair<int, int>(Base::F(sfI, 2), Base::F(sfI, 0)));
            }
        }
        edges = std::vector<std::pair<int, int>>(edgesSet.begin(), edgesSet.end());
    }
}

template <int dim>
void MeshCO<dim>::evaluateConstraint(const Mesh<dim>& mesh,
    const MMCVID& MMCVIDI, double& val, double coef) const
{
    if (MMCVIDI[0] >= 0) {
        // edge-edge ++++
        d_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), val);
    }
    else {
        // point and triangle and degenerate edge-edge
        if (MMCVIDI[1] >= 0) {
            int v0I = -MMCVIDI[0] - 1;
            if (MMCVIDI[2] < 0) {
                // point-point -+--
                d_PP(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), val);
            }
            else if (MMCVIDI[3] < 0) {
                // point-edge -++-
                d_PE(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), val);
            }
            else {
                // point-triangle -+++
                d_PT(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), val);
            }
        }
        else {
            if (MMCVIDI[2] < 0) {
                // triangle-point ---+
                d_PT(Base::V.row(MMCVIDI[3]), mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(-MMCVIDI[1] - 1), mesh.V.row(-MMCVIDI[2] - 1), val);
            }
            else {
                // edge-point --+-
                d_PE(Base::V.row(MMCVIDI[2]), mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(-MMCVIDI[1] - 1), val);
            }
        }
    }
    val *= coef;
}

template <int dim>
void MeshCO<dim>::leftMultiplyConstraintJacobianT(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& activeSet,
    const Eigen::VectorXd& input,
    Eigen::VectorXd& output_incremental,
    double coef) const
{
    // TODO: parallelize

    if constexpr (dim == 3) {
        int constraintI = 0;
        for (const auto& MMCVIDI : activeSet) {
            if (MMCVIDI[0] >= 0) {
                // edge-edge
                Eigen::Matrix<double, dim * 4, 1> g;
                g_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), g);
                g *= coef * input[constraintI];

                output_incremental.segment<dim>(MMCVIDI[0] * dim) += g.template segment<dim>(0);
                output_incremental.segment<dim>(MMCVIDI[1] * dim) += g.template segment<dim>(dim);
            }
            else {
                // point-triangle and degenerate edge-edge
                if (MMCVIDI[1] >= 0) {
                    int v0I = -MMCVIDI[0] - 1;
                    if (MMCVIDI[2] < 0) {
                        // PP
                        Eigen::Matrix<double, dim * 2, 1> g;
                        g_PP(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), g);
                        g *= coef * -MMCVIDI[3] * input[constraintI];

                        output_incremental.segment<dim>(v0I * dim) += g.template segment<dim>(0);
                    }
                    else if (MMCVIDI[3] < 0) {
                        // PE
                        Eigen::Matrix<double, dim * 3, 1> g;
                        g_PE(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), g);
                        g *= coef * -MMCVIDI[3] * input[constraintI];

                        output_incremental.segment<dim>(v0I * dim) += g.template segment<dim>(0);
                    }
                    else {
                        // PT
                        Eigen::Matrix<double, dim * 4, 1> g;
                        g_PT(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), g);
                        g *= coef * input[constraintI];

                        output_incremental.segment<dim>(v0I * dim) += g.template segment<dim>(0);
                    }
                }
                else {
                    if (MMCVIDI[2] < 0) {
                        // triangle-point
                        Eigen::Matrix<double, dim * 4, 1> g;
                        g_PT(Base::V.row(MMCVIDI[3]), mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(-MMCVIDI[1] - 1), mesh.V.row(-MMCVIDI[2] - 1), g);
                        g *= coef * input[constraintI];

                        output_incremental.segment<dim>((-MMCVIDI[0] - 1) * dim) += g.template segment<dim>(dim);
                        output_incremental.segment<dim>((-MMCVIDI[1] - 1) * dim) += g.template segment<dim>(dim * 2);
                        output_incremental.segment<dim>((-MMCVIDI[2] - 1) * dim) += g.template segment<dim>(dim * 3);
                    }
                    else {
                        // edge-point
                        Eigen::Matrix<double, dim * 3, 1> g;
                        g_PE(Base::V.row(MMCVIDI[2]), mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(-MMCVIDI[1] - 1), g);
                        g *= coef * -MMCVIDI[3] * input[constraintI];

                        output_incremental.segment<dim>((-MMCVIDI[0] - 1) * dim) += g.template segment<dim>(dim);
                        output_incremental.segment<dim>((-MMCVIDI[1] - 1) * dim) += g.template segment<dim>(dim * 2);
                    }
                }
            }

            ++constraintI;
        }
    }
    else {
        // TODO: 2D collision
    }
}

template <int dim>
void MeshCO<dim>::evaluateConstraintQP(
    const Mesh<dim>& mesh, // FEM mesh
    const MMCVID& MMCVIDI, // Constraint pair
    const CollisionConstraintType constraintType, // Constraint type to use
    double& val, // Computed value.
    double coef) const // Coefficent to multiply the value by.
{
    double toi;
    try {
        toi = mmcvid_to_toi.at(MMCVIDI);
    }
    catch (std::out_of_range err) {
        spdlog::error("where=\"MeshCO:evaluateConstraintsQP\" msg=\"Cannot compute the constraint value of MMCVID without toi!\"");
        spdlog::debug("MMCVIDI={:s}", MMCVIDI.str());
        spdlog::debug("mmcvid_to_toi:");
        for (auto const& [mmcvid, toi] : mmcvid_to_toi) {
            spdlog::debug("{:s}: {:g}", mmcvid.str(), toi);
        }
        throw "missing toi for constraint";
    }
    if constexpr (dim == 3) {
        if (MMCVIDI[0] >= 0) {
            // edge-edge ++++
            compute_collision_constraint(
                mesh.V_prev.row(MMCVIDI[0]),
                mesh.V_prev.row(MMCVIDI[1]),
                Base::V.row(MMCVIDI[2]),
                Base::V.row(MMCVIDI[3]),
                mesh.V.row(MMCVIDI[0]),
                mesh.V.row(MMCVIDI[1]),
                Base::V.row(MMCVIDI[2]),
                Base::V.row(MMCVIDI[3]),
                constraintType, /*is_edge_edge=*/true, toi,
                val);
            val *= coef;
        }
        else {
            // point and triangle and degenerate edge-edge
            if (MMCVIDI[1] >= 0) {
                int v0I = -MMCVIDI[0] - 1;
                if (MMCVIDI[2] < 0) {
                    // point-point -+--
                    spdlog::warn(
                        "SQP does not handle point-point collisions explicitly.");
                }
                else if (MMCVIDI[3] < 0) {
                    // point-edge -++-
                    spdlog::warn(
                        "SQP does not handle point-edge collisions explicitly.");
                }
                else {
                    // point-triangle -+++
                    compute_collision_constraint(
                        mesh.V_prev.row(v0I),
                        Base::V.row(MMCVIDI[1]),
                        Base::V.row(MMCVIDI[2]),
                        Base::V.row(MMCVIDI[3]),
                        mesh.V.row(v0I),
                        Base::V.row(MMCVIDI[1]),
                        Base::V.row(MMCVIDI[2]),
                        Base::V.row(MMCVIDI[3]),
                        constraintType, /*is_edge_edge=*/false, toi,
                        val);
                    val *= coef;
                }
            }
            else {
                if (MMCVIDI[2] < 0) {
                    // triangle-point ---+
                    compute_collision_constraint(
                        Base::V.row(MMCVIDI[3]),
                        mesh.V_prev.row(-MMCVIDI[0] - 1),
                        mesh.V_prev.row(-MMCVIDI[1] - 1),
                        mesh.V_prev.row(-MMCVIDI[2] - 1),
                        Base::V.row(MMCVIDI[3]),
                        mesh.V.row(-MMCVIDI[0] - 1),
                        mesh.V.row(-MMCVIDI[1] - 1),
                        mesh.V.row(-MMCVIDI[2] - 1),
                        constraintType, /*is_edge_edge=*/false, toi,
                        val);
                    val *= coef;
                }
                else {
                    // edge-point --+-
                    spdlog::warn(
                        "SQP does not handle edge-point collisions explicitly.");
                }
            }
        }
    }
    else {
        // TODO: 2D collisions
        spdlog::error("IPC does not handle 2D collisions.");
    }
}

template <int dim>
void MeshCO<dim>::leftMultiplyConstraintJacobianTQP(
    const Mesh<dim>& mesh,
    const std::vector<MMCVID>& activeSet,
    const Eigen::VectorXd& input,
    const CollisionConstraintType constraintType,
    Eigen::VectorXd& output_incremental,
    double coef) const
{
    // TODO: parallelize
    if constexpr (dim == 3) {
        int constraintI = 0;
        for (const auto& MMCVIDI : activeSet) {
            double toi;
            try {
                toi = mmcvid_to_toi.at(MMCVIDI);
            }
            catch (std::out_of_range err) {
                spdlog::error("where=\"MeshCO:evaluateConstraintsQP\" msg=\"Cannot compute the constraint value of MMCVID without toi!\"");
                spdlog::debug("MMCVIDI={:s}", MMCVIDI.str());
                spdlog::debug("mmcvid_to_toi:");
                for (auto const& [mmcvid, toi] : mmcvid_to_toi) {
                    spdlog::debug("{:s}: {:g}", mmcvid.str(), toi);
                }
                throw "missing toi for constraint";
            }

            if (MMCVIDI[0] >= 0) {
                // edge-edge
                Eigen::Matrix<double, dim * 4, 1> grad;
                compute_collision_constraint_gradient(
                    mesh.V_prev.row(MMCVIDI[0]),
                    mesh.V_prev.row(MMCVIDI[1]),
                    Base::V.row(MMCVIDI[2]),
                    Base::V.row(MMCVIDI[3]),
                    mesh.V.row(MMCVIDI[0]),
                    mesh.V.row(MMCVIDI[1]),
                    Base::V.row(MMCVIDI[2]),
                    Base::V.row(MMCVIDI[3]),
                    constraintType, /*is_edge_edge=*/true, toi,
                    grad);
                grad *= coef * input[constraintI];
                output_incremental.segment<dim>(MMCVIDI[0] * dim) += grad.template segment<dim>(0);
                output_incremental.segment<dim>(MMCVIDI[1] * dim) += grad.template segment<dim>(dim);
            }
            else {
                // point-triangle and degenerate edge-edge
                if (MMCVIDI[1] >= 0) {
                    int v0I = -MMCVIDI[0] - 1;
                    if (MMCVIDI[2] < 0) {
                        // PP
                    }
                    else if (MMCVIDI[3] < 0) {
                        // PE
                    }
                    else {
                        // PT
                        Eigen::Matrix<double, dim * 4, 1> grad;
                        compute_collision_constraint_gradient(
                            mesh.V_prev.row(v0I),
                            Base::V.row(MMCVIDI[1]),
                            Base::V.row(MMCVIDI[2]),
                            Base::V.row(MMCVIDI[3]),
                            mesh.V.row(v0I),
                            Base::V.row(MMCVIDI[1]),
                            Base::V.row(MMCVIDI[2]),
                            Base::V.row(MMCVIDI[3]),
                            constraintType, /*is_edge_edge=*/false, toi,
                            grad);
                        grad *= coef * input[constraintI];
                        output_incremental.segment<dim>(v0I * dim) += grad.template segment<dim>(0);
                    }
                }
                else {
                    if (MMCVIDI[2] < 0) {
                        // triangle-point
                        Eigen::Matrix<double, dim * 4, 1> grad;
                        compute_collision_constraint_gradient(
                            Base::V.row(MMCVIDI[3]),
                            mesh.V_prev.row(-MMCVIDI[0] - 1),
                            mesh.V_prev.row(-MMCVIDI[1] - 1),
                            mesh.V_prev.row(-MMCVIDI[2] - 1),
                            Base::V.row(MMCVIDI[3]),
                            mesh.V.row(-MMCVIDI[0] - 1),
                            mesh.V.row(-MMCVIDI[1] - 1),
                            mesh.V.row(-MMCVIDI[2] - 1),
                            constraintType, /*is_edge_edge=*/false, toi,
                            grad);
                        grad *= coef * input[constraintI];
                        output_incremental.segment<dim>((-MMCVIDI[0] - 1) * dim) += grad.template segment<dim>(dim);
                        output_incremental.segment<dim>((-MMCVIDI[1] - 1) * dim) += grad.template segment<dim>(2 * dim);
                        output_incremental.segment<dim>((-MMCVIDI[2] - 1) * dim) += grad.template segment<dim>(3 * dim);
                    }
                    else {
                        // edge-point
                    }
                }
            }

            ++constraintI;
        }
    }
    else {
        // TODO: 2D collision
    }
}

template <int dim>
void MeshCO<dim>::augmentIPHessian(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& activeSet,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* mtr_incremental,
    double dHat, double coef, bool projectDBC) const
{
    // TODO: parallelize

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
                    d_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), d);
                    Eigen::Matrix<double, dim * 4, 1> g;
                    g_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), g);
                    Eigen::Matrix<double, dim * 4, dim * 4> H;
                    H_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), H);

                    double g_b, H_b;
                    compute_g_b(d, dHat, g_b);
                    compute_H_b(d, dHat, H_b);

                    IPHessian[cI] = ((coef * H_b) * g) * g.transpose() + (coef * g_b) * H;
                    IglUtils::makePD(IPHessian[cI]);

                    rowIStart[cI][0] = mesh.isProjectDBCVertex(MMCVIDI[0], projectDBC) ? -1 : (MMCVIDI[0] * dim);
                    rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                    rowIStart[cI][2] = -1;
                    rowIStart[cI][3] = -1;
                }
                else {
                    // point-triangle and degenerate edge-edge
                    if (MMCVIDI[1] >= 0) {
                        int v0I = -MMCVIDI[0] - 1;
                        if (MMCVIDI[2] < 0) {
                            // PP
                            double d;
                            d_PP(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), d);
                            Eigen::Matrix<double, dim * 2, 1> g;
                            g_PP(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), g);
                            Eigen::Matrix<double, dim * 2, dim * 2> H;
                            H_PP(H);

                            double g_b, H_b;
                            compute_g_b(d, dHat, g_b);
                            compute_H_b(d, dHat, H_b);

                            double coef_dup = coef * -MMCVIDI[3];
                            Eigen::Matrix<double, dim * 2, dim* 2> HessianBlock = ((coef_dup * H_b) * g) * g.transpose() + (coef_dup * g_b) * H;
                            IglUtils::makePD(HessianBlock);
                            IPHessian[cI].template block<dim, dim>(0, 0) = HessianBlock.template block<dim, dim>(0, 0);
                        }
                        else if (MMCVIDI[3] < 0) {
                            // PE
                            double d;
                            d_PE(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), d);
                            Eigen::Matrix<double, dim * 3, 1> g;
                            g_PE(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), g);
                            Eigen::Matrix<double, dim * 3, dim * 3> H;
                            H_PE(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), H);

                            double g_b, H_b;
                            compute_g_b(d, dHat, g_b);
                            compute_H_b(d, dHat, H_b);

                            double coef_dup = coef * -MMCVIDI[3];
                            Eigen::Matrix<double, dim * 3, dim* 3> HessianBlock = ((coef_dup * H_b) * g) * g.transpose() + (coef_dup * g_b) * H;
                            IglUtils::makePD(HessianBlock);
                            IPHessian[cI].template block<dim, dim>(0, 0) = HessianBlock.template block<dim, dim>(0, 0);
                        }
                        else {
                            // PT
                            double d;
                            d_PT(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), d);
                            Eigen::Matrix<double, dim * 4, 1> g;
                            g_PT(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), g);
                            Eigen::Matrix<double, dim * 4, dim * 4> H;
                            H_PT(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), H);

                            double g_b, H_b;
                            compute_g_b(d, dHat, g_b);
                            compute_H_b(d, dHat, H_b);

                            IPHessian[cI] = ((coef * H_b) * g) * g.transpose() + (coef * g_b) * H;
                            IglUtils::makePD(IPHessian[cI]);
                        }
                        rowIStart[cI][0] = mesh.isProjectDBCVertex(v0I, projectDBC) ? -1 : (v0I * dim);
                        rowIStart[cI][1] = -1;
                        rowIStart[cI][2] = -1;
                        rowIStart[cI][3] = -1;
                    }
                    else {
                        if (MMCVIDI[2] < 0) {
                            // triangle-point
                            double d;
                            d_PT(Base::V.row(MMCVIDI[3]), mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(-MMCVIDI[1] - 1), mesh.V.row(-MMCVIDI[2] - 1), d);
                            Eigen::Matrix<double, dim * 4, 1> g;
                            g_PT(Base::V.row(MMCVIDI[3]), mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(-MMCVIDI[1] - 1), mesh.V.row(-MMCVIDI[2] - 1), g);
                            Eigen::Matrix<double, dim * 4, dim * 4> H;
                            H_PT(Base::V.row(MMCVIDI[3]), mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(-MMCVIDI[1] - 1), mesh.V.row(-MMCVIDI[2] - 1), H);

                            double g_b, H_b;
                            compute_g_b(d, dHat, g_b);
                            compute_H_b(d, dHat, H_b);

                            IPHessian[cI] = ((coef * H_b) * g) * g.transpose() + (coef * g_b) * H;
                            IglUtils::makePD(IPHessian[cI]);

                            rowIStart[cI][0] = -1;
                            rowIStart[cI][1] = mesh.isProjectDBCVertex(-MMCVIDI[0] - 1, projectDBC) ? -1 : ((-MMCVIDI[0] - 1) * dim);
                            rowIStart[cI][2] = mesh.isProjectDBCVertex(-MMCVIDI[1] - 1, projectDBC) ? -1 : ((-MMCVIDI[1] - 1) * dim);
                            rowIStart[cI][3] = mesh.isProjectDBCVertex(-MMCVIDI[2] - 1, projectDBC) ? -1 : ((-MMCVIDI[2] - 1) * dim);
                        }
                        else {
                            // edge-point
                            double d;
                            d_PE(Base::V.row(MMCVIDI[2]), mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(-MMCVIDI[1] - 1), d);
                            Eigen::Matrix<double, dim * 3, 1> g;
                            g_PE(Base::V.row(MMCVIDI[2]), mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(-MMCVIDI[1] - 1), g);
                            Eigen::Matrix<double, dim * 3, dim * 3> H;
                            H_PE(Base::V.row(MMCVIDI[2]), mesh.V.row(-MMCVIDI[0] - 1), mesh.V.row(-MMCVIDI[1] - 1), H);

                            double g_b, H_b;
                            compute_g_b(d, dHat, g_b);
                            compute_H_b(d, dHat, H_b);

                            double coef_dup = coef * -MMCVIDI[3];
                            Eigen::Matrix<double, dim * 3, dim* 3> HessianBlock = ((coef_dup * H_b) * g) * g.transpose() + (coef_dup * g_b) * H;
                            IglUtils::makePD(HessianBlock);
                            IPHessian[cI].template block<dim * 2, dim * 2>(dim, dim) = HessianBlock.template block<dim * 2, dim * 2>(dim, dim);

                            rowIStart[cI][0] = -1;
                            rowIStart[cI][1] = mesh.isProjectDBCVertex(-MMCVIDI[0] - 1, projectDBC) ? -1 : ((-MMCVIDI[0] - 1) * dim);
                            rowIStart[cI][2] = mesh.isProjectDBCVertex(-MMCVIDI[1] - 1, projectDBC) ? -1 : ((-MMCVIDI[1] - 1) * dim);
                            rowIStart[cI][3] = -1;
                        }
                    }
                }
            }
#ifdef USE_TBB
        );
#endif

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
void MeshCO<dim>::largestFeasibleStepSize(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    const Eigen::VectorXd& searchDir,
    double slackness,
    const std::vector<std::pair<int, int>>& constraintSet,
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
                if (constraintSet[cI].first < 0) {
                    if (constraintSet[cI].second < 0) {
                        // TP
                        int vI = -constraintSet[cI].second - 1;
                        const RowVector3i& sfVInd = mesh.SF.row(-constraintSet[cI].first - 1);

                        double d_sqrt;
                        computePointTriD(Base::V.row(vI), mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), d_sqrt);
                        d_sqrt = std::sqrt(d_sqrt);

                        largestAlphasAS[cI] = 1.0;
                        if (CTCD::vertexFaceCTCD(Base::V.row(vI).transpose(),
                                mesh.V.row(sfVInd[0]).transpose(),
                                mesh.V.row(sfVInd[1]).transpose(),
                                mesh.V.row(sfVInd[2]).transpose(),
                                Base::V.row(vI).transpose(),
                                mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                                mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                                mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                                CCDDistRatio * d_sqrt,
                                largestAlphasAS[cI])) {
                            if (largestAlphasAS[cI] < 1.0e-6) {
                                if (CTCD::vertexFaceCTCD(Base::V.row(vI).transpose(),
                                        mesh.V.row(sfVInd[0]).transpose(),
                                        mesh.V.row(sfVInd[1]).transpose(),
                                        mesh.V.row(sfVInd[2]).transpose(),
                                        Base::V.row(vI).transpose(),
                                        mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                                        mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                                        mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                                        0.0, largestAlphasAS[cI])) {
                                    largestAlphasAS[cI] *= slackness;
                                }
                                else {
                                    largestAlphasAS[cI] = 1.0;
                                }
                            }
                        }
                    }
                    else {
                        // PT
                        int vI = mesh.SVI[-constraintSet[cI].first - 1];
                        const RowVector3i& sfVInd = Base::F.row(constraintSet[cI].second);

                        double d_sqrt;
                        computePointTriD(mesh.V.row(vI), Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]), d_sqrt);
                        d_sqrt = std::sqrt(d_sqrt);

                        largestAlphasAS[cI] = 1.0;
                        if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                                Base::V.row(sfVInd[0]).transpose(),
                                Base::V.row(sfVInd[1]).transpose(),
                                Base::V.row(sfVInd[2]).transpose(),
                                mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                                Base::V.row(sfVInd[0]).transpose(),
                                Base::V.row(sfVInd[1]).transpose(),
                                Base::V.row(sfVInd[2]).transpose(),
                                CCDDistRatio * d_sqrt,
                                largestAlphasAS[cI])) {
                            if (largestAlphasAS[cI] < 1.0e-6) {
                                if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                                        Base::V.row(sfVInd[0]).transpose(),
                                        Base::V.row(sfVInd[1]).transpose(),
                                        Base::V.row(sfVInd[2]).transpose(),
                                        mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                                        Base::V.row(sfVInd[0]).transpose(),
                                        Base::V.row(sfVInd[1]).transpose(),
                                        Base::V.row(sfVInd[2]).transpose(),
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
                else {
                    // EE
                    const auto& meshEI = mesh.SFEdges[constraintSet[cI].first];
                    const auto& meshEJ = edges[constraintSet[cI].second];

                    double d_sqrt;
                    computeEdgeEdgeD(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                        Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), d_sqrt);
                    d_sqrt = std::sqrt(d_sqrt);

                    largestAlphasAS[cI] = 1.0;
                    if (CTCD::edgeEdgeCTCD(mesh.V.row(meshEI.first).transpose(),
                            mesh.V.row(meshEI.second).transpose(),
                            Base::V.row(meshEJ.first).transpose(),
                            Base::V.row(meshEJ.second).transpose(),
                            mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                            mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                            Base::V.row(meshEJ.first).transpose(),
                            Base::V.row(meshEJ.second).transpose(),
                            CCDDistRatio * d_sqrt,
                            largestAlphasAS[cI])) {
                        if (largestAlphasAS[cI] < 1.0e-6) {
                            if (CTCD::edgeEdgeCTCD(mesh.V.row(meshEI.first).transpose(),
                                    mesh.V.row(meshEI.second).transpose(),
                                    Base::V.row(meshEJ.first).transpose(),
                                    Base::V.row(meshEJ.second).transpose(),
                                    mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                                    mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                                    Base::V.row(meshEJ.first).transpose(),
                                    Base::V.row(meshEJ.second).transpose(),
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
            }
#ifdef USE_TBB
        );
#endif
        stepSize = std::min(stepSize, largestAlphasAS.minCoeff());
    }
    return;
#endif

    largestFeasibleStepSize_CCD(mesh, sh, searchDir, slackness, stepSize);
}

#ifdef IPC_WITH_TIGHT_INCLUSION
template <int dim>
void MeshCO<dim>::largestFeasibleStepSize_TightInclusion(
    const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    const Eigen::VectorXd& searchDir,
    double tolerance,
    const std::vector<std::pair<int, int>>& constraintSet,
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
        if (constraintSet[cI].first < 0) {
            if (constraintSet[cI].second < 0) {
                // TP
                int vI = -constraintSet[cI].second - 1;
                const RowVector3i& sfVInd = mesh.SF.row(-constraintSet[cI].first - 1);

                double d_sqrt; // d is squared distance
                computePointTriD(Base::V.row(vI), mesh.V.row(sfVInd[0]),
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
                    Base::V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose(),
                    mesh.V.row(sfVInd[1]).transpose(),
                    mesh.V.row(sfVInd[2]).transpose(),
                    Base::V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                    mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                    mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                    /*err=*/tight_inclusion_vf_err,
                    /*ms=*/std::min(TIGHT_INCLUSION_DIST_P * d_sqrt, TIGHT_INCLUSION_MIN_DIST),
                    toi,
                    tolerance,
                    /*t_max=*/stepSize,
                    /*max_itr=*/TIGHT_INCLUSION_MAX_ITER,
                    output_tolerance,
                    /*CCD_TYPE=*/TIGHT_INCLUSION_CCD_TYPE,
                    /*no_zero_toi=*/TIGHT_INCLUSION_NO_ZERO_TOI);

                if (has_collision && toi < 1e-6) {
                    has_collision = inclusion_ccd::vertexFaceCCD_double(
                        Base::V.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose(),
                        mesh.V.row(sfVInd[1]).transpose(),
                        mesh.V.row(sfVInd[2]).transpose(),
                        Base::V.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                        mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                        mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                        /*err=*/tight_inclusion_vf_err,
                        /*ms=*/0,
                        toi,
                        tolerance,
                        /*max_t=*/stepSize,
                        /* max_itr=*/-1,
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
            else {
                // PT
                int vI = mesh.SVI[-constraintSet[cI].first - 1];
                const RowVector3i& sfVInd = Base::F.row(constraintSet[cI].second);

                double d_sqrt;
                computePointTriD(mesh.V.row(vI), Base::V.row(sfVInd[0]),
                    Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]), d_sqrt);
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
                    Base::V.row(sfVInd[0]).transpose(),
                    Base::V.row(sfVInd[1]).transpose(),
                    Base::V.row(sfVInd[2]).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                    Base::V.row(sfVInd[0]).transpose(),
                    Base::V.row(sfVInd[1]).transpose(),
                    Base::V.row(sfVInd[2]).transpose(),
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
                        Base::V.row(sfVInd[0]).transpose(),
                        Base::V.row(sfVInd[1]).transpose(),
                        Base::V.row(sfVInd[2]).transpose(),
                        mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                        Base::V.row(sfVInd[0]).transpose(),
                        Base::V.row(sfVInd[1]).transpose(),
                        Base::V.row(sfVInd[2]).transpose(),
                        /*err=*/tight_inclusion_vf_err,
                        /*ms=*/0,
                        toi,
                        tolerance,
                        /*max_t=*/stepSize,
                        /* max_itr=*/-1,
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
        else {
            // EE
            const auto& meshEI = mesh.SFEdges[constraintSet[cI].first];
            const auto& meshEJ = edges[constraintSet[cI].second];

            double d_sqrt;
            computeEdgeEdgeD(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), d_sqrt);
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
                mesh.V.row(meshEI.first).transpose(),
                mesh.V.row(meshEI.second).transpose(),
                Base::V.row(meshEJ.first).transpose(),
                Base::V.row(meshEJ.second).transpose(),
                mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                Base::V.row(meshEJ.first).transpose(),
                Base::V.row(meshEJ.second).transpose(),
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
                double output_tolerance;
                bool has_collision = inclusion_ccd::vertexFaceCCD_double(
                    mesh.V.row(meshEI.first).transpose(),
                    mesh.V.row(meshEI.second).transpose(),
                    Base::V.row(meshEJ.first).transpose(),
                    Base::V.row(meshEJ.second).transpose(),
                    mesh.V.row(meshEI.first).transpose() + searchDir.segment<dim>(meshEI.first * dim),
                    mesh.V.row(meshEI.second).transpose() + searchDir.segment<dim>(meshEI.second * dim),
                    Base::V.row(meshEJ.first).transpose(),
                    Base::V.row(meshEJ.second).transpose(),
                    /*err=*/tight_inclusion_ee_err,
                    /*ms=*/0,
                    toi,
                    tolerance,
                    /*max_t=*/stepSize,
                    /* max_itr=*/-1,
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
    largestFeasibleStepSize_CCD(mesh, sh, searchDir, tolerance, stepSize);
#endif
}
#endif // IPC_WITH_TIGHT_INCLUSION

template <int dim>
void MeshCO<dim>::largestFeasibleStepSize_exact(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    const Eigen::VectorXd& searchDir,
    const ccd::CCDMethod ccdMethod,
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
                largestAlphasAS[cI] = stepSize;
                if (constraintSet[cI].first < 0) {
                    if (constraintSet[cI].second < 0) {
                        // TP
                        int vI = -constraintSet[cI].second - 1;
                        const RowVector3i& sfVInd = mesh.SF.row(-constraintSet[cI].first - 1);

                        while (ccd::vertexFaceCCD(
                            Base::V.row(vI).transpose(),
                            mesh.V.row(sfVInd[0]).transpose(),
                            mesh.V.row(sfVInd[1]).transpose(),
                            mesh.V.row(sfVInd[2]).transpose(),
                            Base::V.row(vI).transpose(),
                            mesh.V.row(sfVInd[0]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(sfVInd[0] * dim),
                            mesh.V.row(sfVInd[1]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(sfVInd[1] * dim),
                            mesh.V.row(sfVInd[2]).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(sfVInd[2] * dim),
                            ccdMethod)) {
                            largestAlphasAS[cI] /= 2.0;
                        }
                    }
                    else {
                        // PT
                        int vI = mesh.SVI[-constraintSet[cI].first - 1];
                        const RowVector3i& sfVInd = Base::F.row(constraintSet[cI].second);

                        while (ccd::vertexFaceCCD(
                            mesh.V.row(vI).transpose(),
                            Base::V.row(sfVInd[0]).transpose(),
                            Base::V.row(sfVInd[1]).transpose(),
                            Base::V.row(sfVInd[2]).transpose(),
                            mesh.V.row(vI).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(vI * dim),
                            Base::V.row(sfVInd[0]).transpose(),
                            Base::V.row(sfVInd[1]).transpose(),
                            Base::V.row(sfVInd[2]).transpose(),
                            ccdMethod)) {
                            largestAlphasAS[cI] /= 2.0;
                        }
                    }
                }
                else {
                    // EE
                    const auto& meshEI = mesh.SFEdges[constraintSet[cI].first];
                    const auto& meshEJ = edges[constraintSet[cI].second];

                    while (ccd::edgeEdgeCCD(
                        mesh.V.row(meshEI.first).transpose(),
                        mesh.V.row(meshEI.second).transpose(),
                        Base::V.row(meshEJ.first).transpose(),
                        Base::V.row(meshEJ.second).transpose(),
                        mesh.V.row(meshEI.first).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(meshEI.first * dim),
                        mesh.V.row(meshEI.second).transpose() + largestAlphasAS[cI] * searchDir.segment<dim>(meshEI.second * dim),
                        Base::V.row(meshEJ.first).transpose(),
                        Base::V.row(meshEJ.second).transpose(),
                        ccdMethod)) {
                        largestAlphasAS[cI] /= 2.0;
                    }
                }
            }
#ifdef USE_TBB
        );
#endif
        stepSize = std::min(stepSize, largestAlphasAS.minCoeff());
    }
    return;
#endif

    largestFeasibleStepSize_CCD_exact(mesh, sh, searchDir, ccdMethod, stepSize);
}

template <int dim>
void MeshCO<dim>::largestFeasibleStepSize_CCD(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh, const Eigen::VectorXd& searchDir,
    double slackness, double& stepSize)
{
    const double CCDDistRatio = 1.0 - slackness;

    // point-triangle
    Eigen::VectorXd largestAlphaPT(Base::F.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::F.rows(), 1, [&](int sfI)
#else
    for (int sfI = 0; sfI < Base::F.rows(); ++sfI)
#endif
        {
            largestAlphaPT[sfI] = 1.0;
            const RowVector3i& sfVInd = Base::F.row(sfI);
#ifdef USE_SH_LFSS
            std::unordered_set<int> pointInds; // NOTE: different constraint order will result in numerically different results
            sh.queryTriangleForPoints(Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]),
                0.0, pointInds);
            for (const auto& svI : pointInds)
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
            {
                int vI = mesh.SVI[svI];

                double d_sqrt;
                computePointTriD(mesh.V.row(vI), Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlpha = 1.0;
                if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                        Base::V.row(sfVInd[0]).transpose(),
                        Base::V.row(sfVInd[1]).transpose(),
                        Base::V.row(sfVInd[2]).transpose(),
                        mesh.V.row(vI).transpose() + searchDir.template segment<dim>(vI * dim),
                        Base::V.row(sfVInd[0]).transpose(),
                        Base::V.row(sfVInd[1]).transpose(),
                        Base::V.row(sfVInd[2]).transpose(),
                        CCDDistRatio * d_sqrt,
                        largestAlpha)) {
                    if (largestAlpha < 1.0e-6) {
                        if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                                Base::V.row(sfVInd[0]).transpose(),
                                Base::V.row(sfVInd[1]).transpose(),
                                Base::V.row(sfVInd[2]).transpose(),
                                mesh.V.row(vI).transpose() + searchDir.template segment<dim>(vI * dim),
                                Base::V.row(sfVInd[0]).transpose(),
                                Base::V.row(sfVInd[1]).transpose(),
                                Base::V.row(sfVInd[2]).transpose(),
                                0.0, largestAlpha)) {
                            largestAlpha *= slackness;
                        }
                        else {
                            largestAlpha = 1.0;
                        }
                    }
                }
                largestAlphaPT[sfI] = std::min(largestAlpha, largestAlphaPT[sfI]);
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphaPT.minCoeff());

    // triangle-point, edge-point, point-point
    Eigen::VectorXd largestAlphaTP(Base::V.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::V.rows(), 1, [&](int vI)
#else
    for (int vI = 0; vI < Base::V.rows(); ++vI)
#endif
        {
            largestAlphaTP[vI] = 1.0;
#ifdef USE_SH_LFSS
            std::unordered_set<int> svInds, edgeInds, triInds; // NOTE: different constraint order will result in numerically different results
            sh.queryPointForPrimitives(Base::V.row(vI), Eigen::Matrix<double, 1, dim>::Zero(),
                svInds, edgeInds, triInds);

            // triangle-point
            for (const auto& sfI : triInds)
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
            {
                const RowVector3i& sfVInd = mesh.SF.row(sfI);

                double d_sqrt;
                computePointTriD(Base::V.row(vI), mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlpha = 1.0;
                if (CTCD::vertexFaceCTCD(Base::V.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose(),
                        mesh.V.row(sfVInd[1]).transpose(),
                        mesh.V.row(sfVInd[2]).transpose(),
                        Base::V.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose() + searchDir.template segment<dim>(sfVInd[0] * dim),
                        mesh.V.row(sfVInd[1]).transpose() + searchDir.template segment<dim>(sfVInd[1] * dim),
                        mesh.V.row(sfVInd[2]).transpose() + searchDir.template segment<dim>(sfVInd[2] * dim),
                        CCDDistRatio * d_sqrt,
                        largestAlpha)) {
                    if (largestAlpha < 1.0e-6) {
                        if (CTCD::vertexFaceCTCD(Base::V.row(vI).transpose(),
                                mesh.V.row(sfVInd[0]).transpose(),
                                mesh.V.row(sfVInd[1]).transpose(),
                                mesh.V.row(sfVInd[2]).transpose(),
                                Base::V.row(vI).transpose(),
                                mesh.V.row(sfVInd[0]).transpose() + searchDir.template segment<dim>(sfVInd[0] * dim),
                                mesh.V.row(sfVInd[1]).transpose() + searchDir.template segment<dim>(sfVInd[1] * dim),
                                mesh.V.row(sfVInd[2]).transpose() + searchDir.template segment<dim>(sfVInd[2] * dim),
                                0.0, largestAlpha)) {
                            largestAlpha *= slackness;
                        }
                        else {
                            largestAlpha = 1.0;
                        }
                    }
                }
                largestAlphaTP[vI] = std::min(largestAlpha, largestAlphaTP[vI]);
            }

        // point-point
#ifdef USE_SH_LFSS
            for (const auto& svJ : svInds)
#else
        for (int svJ = 0; svJ < mesh.SVI.size(); ++svJ)
#endif
            {
                int vJ = mesh.SVI[svJ];

                double d_sqrt = (Base::V.row(vI) - mesh.V.row(vJ)).norm();

                double largestAlpha = 1.0;
                if (CTCD::vertexVertexCTCD(Base::V.row(vI).transpose(),
                        mesh.V.row(vJ).transpose(),
                        Base::V.row(vI).transpose(),
                        mesh.V.row(vJ).transpose() + searchDir.template segment<dim>(vJ * dim),
                        CCDDistRatio * d_sqrt,
                        largestAlpha)) {
                    if (largestAlpha < 1.0e-6) {
                        if (CTCD::vertexVertexCTCD(Base::V.row(vI).transpose(),
                                mesh.V.row(vJ).transpose(),
                                Base::V.row(vI).transpose(),
                                mesh.V.row(vJ).transpose() + searchDir.template segment<dim>(vJ * dim),
                                0.0, largestAlpha)) {
                            largestAlpha *= slackness;
                        }
                        else {
                            largestAlpha = 1.0;
                        }
                    }
                }
                largestAlphaTP[vI] = std::min(largestAlpha, largestAlphaTP[vI]);
            }

        // edge-point
#ifdef USE_SH_LFSS
            for (const auto& eI : edgeInds) {
                const auto& meshEI = mesh.SFEdges[eI];
#else
        for (const auto& meshEI : mesh.SFEdges) {
#endif
                double d_sqrt;
                computePointEdgeD(Base::V.row(vI), mesh.V.row(meshEI.first),
                    mesh.V.row(meshEI.second), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlpha = 1.0;
                if (CTCD::vertexEdgeCTCD(Base::V.row(vI).transpose(),
                        mesh.V.row(meshEI.first).transpose(),
                        mesh.V.row(meshEI.second).transpose(),
                        Base::V.row(vI).transpose(),
                        mesh.V.row(meshEI.first).transpose() + searchDir.template segment<dim>(meshEI.first * dim),
                        mesh.V.row(meshEI.second).transpose() + searchDir.template segment<dim>(meshEI.second * dim),
                        CCDDistRatio * d_sqrt,
                        largestAlpha)) {
                    if (largestAlpha < 1.0e-6) {
                        if (CTCD::vertexEdgeCTCD(Base::V.row(vI).transpose(),
                                mesh.V.row(meshEI.first).transpose(),
                                mesh.V.row(meshEI.second).transpose(),
                                Base::V.row(vI).transpose(),
                                mesh.V.row(meshEI.first).transpose() + searchDir.template segment<dim>(meshEI.first * dim),
                                mesh.V.row(meshEI.second).transpose() + searchDir.template segment<dim>(meshEI.second * dim),
                                0.0, largestAlpha)) {
                            largestAlpha *= slackness;
                        }
                        else {
                            largestAlpha = 1.0;
                        }
                    }
                }
                largestAlphaTP[vI] = std::min(largestAlpha, largestAlphaTP[vI]);
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphaTP.minCoeff());

    // edge-edge, point-edge
    Eigen::VectorXd largestAlphaEE(edges.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)edges.size(), 1, [&](int eJ)
#else
    for (int eJ = 0; eJ < edges.size(); ++eJ)
#endif
        {
            largestAlphaEE[eJ] = 1.0;
            const auto& meshEJ = edges[eJ];
#ifdef USE_SH_LFSS
            std::vector<int> svInds, edgeInds; // NOTE: different constraint order will result in numerically different results
            sh.queryEdgeForPE(Base::V.row(meshEJ.first), Base::V.row(meshEJ.second),
                svInds, edgeInds);

            // edge-edge
            for (const auto& eI : edgeInds) {
                const auto& meshEI = mesh.SFEdges[eI];
#else
        for (const auto& meshEI : mesh.SFEdges) {
#endif
                double d_sqrt;
                computeEdgeEdgeD(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                    Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlphas = 1.0;
                if (CTCD::edgeEdgeCTCD(mesh.V.row(meshEI.first).transpose(),
                        mesh.V.row(meshEI.second).transpose(),
                        Base::V.row(meshEJ.first).transpose(),
                        Base::V.row(meshEJ.second).transpose(),
                        mesh.V.row(meshEI.first).transpose() + searchDir.template segment<dim>(meshEI.first * dim),
                        mesh.V.row(meshEI.second).transpose() + searchDir.template segment<dim>(meshEI.second * dim),
                        Base::V.row(meshEJ.first).transpose(),
                        Base::V.row(meshEJ.second).transpose(),
                        CCDDistRatio * d_sqrt,
                        largestAlphas)) {
                    if (largestAlphas < 1.0e-6) {
                        if (CTCD::edgeEdgeCTCD(mesh.V.row(meshEI.first).transpose(),
                                mesh.V.row(meshEI.second).transpose(),
                                Base::V.row(meshEJ.first).transpose(),
                                Base::V.row(meshEJ.second).transpose(),
                                mesh.V.row(meshEI.first).transpose() + searchDir.template segment<dim>(meshEI.first * dim),
                                mesh.V.row(meshEI.second).transpose() + searchDir.template segment<dim>(meshEI.second * dim),
                                Base::V.row(meshEJ.first).transpose(),
                                Base::V.row(meshEJ.second).transpose(),
                                0.0, largestAlphas)) {
                            if (largestAlphas == 0.0) {
                                // numerically parallel edge-edge causes CCD code to fail
                                largestAlphas = 1.0;
                            }
                            largestAlphas *= slackness;
                        }
                        else {
                            largestAlphas = 1.0;
                        }
                    }
                }
                largestAlphaEE[eJ] = std::min(largestAlphas, largestAlphaEE[eJ]);
            }

        // point-edge
#ifdef USE_SH_LFSS
            for (const auto& svI : svInds)
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
            {
                int vI = mesh.SVI[svI];

                double d_sqrt;
                computePointEdgeD(mesh.V.row(vI),
                    Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlphas = 1.0;
                if (CTCD::vertexEdgeCTCD(mesh.V.row(vI).transpose(),
                        Base::V.row(meshEJ.first).transpose(),
                        Base::V.row(meshEJ.second).transpose(),
                        mesh.V.row(vI).transpose() + searchDir.template segment<dim>(vI * dim),
                        Base::V.row(meshEJ.first).transpose(),
                        Base::V.row(meshEJ.second).transpose(),
                        CCDDistRatio * d_sqrt,
                        largestAlphas)) {
                    if (largestAlphas < 1.0e-6) {
                        if (CTCD::vertexEdgeCTCD(mesh.V.row(vI).transpose(),
                                Base::V.row(meshEJ.first).transpose(),
                                Base::V.row(meshEJ.second).transpose(),
                                mesh.V.row(vI).transpose() + searchDir.template segment<dim>(vI * dim),
                                Base::V.row(meshEJ.first).transpose(),
                                Base::V.row(meshEJ.second).transpose(),
                                0.0, largestAlphas)) {
                            if (largestAlphas == 0.0) {
                                largestAlphas = 1.0;
                            }
                            largestAlphas *= slackness;
                        }
                        else {
                            largestAlphas = 1.0;
                        }
                    }
                }
                largestAlphaEE[eJ] = std::min(largestAlphas, largestAlphaEE[eJ]);
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphaEE.minCoeff());
}

#ifdef IPC_WITH_TIGHT_INCLUSION
template <int dim>
void MeshCO<dim>::largestFeasibleStepSize_CCD_TightInclusion(
    const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    const Eigen::VectorXd& searchDir,
    double tolerance, double& stepSize)
{
    // point-triangle
#ifdef USE_TBB
    std::mutex stepSizeLock;
    tbb::parallel_for(0, (int)Base::F.rows(), 1, [&](int sfI) {
#else
    for (int sfI = 0; sfI < Base::F.rows(); ++sfI) {
#endif
        const RowVector3i& sfVInd = Base::F.row(sfI);
#ifdef USE_SH_LFSS
        std::unordered_set<int> pointInds; // NOTE: different constraint order will result in numerically different results
        sh.queryTriangleForPoints(Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]),
            0.0, pointInds);
        for (const auto& svI : pointInds) {
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
#endif
            int vI = mesh.SVI[svI];

            double d_sqrt;
            computePointTriD(mesh.V.row(vI), Base::V.row(sfVInd[0]),
                Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]), d_sqrt);
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
                Base::V.row(sfVInd[0]).transpose(),
                Base::V.row(sfVInd[1]).transpose(),
                Base::V.row(sfVInd[2]).transpose(),
                mesh.V.row(vI).transpose() + searchDir.template segment<dim>(vI * dim),
                Base::V.row(sfVInd[0]).transpose(),
                Base::V.row(sfVInd[1]).transpose(),
                Base::V.row(sfVInd[2]).transpose(),
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
                    Base::V.row(sfVInd[0]).transpose(),
                    Base::V.row(sfVInd[1]).transpose(),
                    Base::V.row(sfVInd[2]).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.template segment<dim>(vI * dim),
                    Base::V.row(sfVInd[0]).transpose(),
                    Base::V.row(sfVInd[1]).transpose(),
                    Base::V.row(sfVInd[2]).transpose(),
                    /*err=*/tight_inclusion_vf_err,
                    /*ms=*/0,
                    toi,
                    tolerance,
                    /*max_t=*/stepSize,
                    /* max_itr=*/-1,
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

// triangle-point
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::V.rows(), 1, [&](int vI) {
#else
    for (int vI = 0; vI < Base::V.rows(); ++vI) {
#endif
#ifdef USE_SH_LFSS
        std::unordered_set<int> svInds, edgeInds, triInds; // NOTE: different constraint order will result in numerically different results
        sh.queryPointForPrimitives(Base::V.row(vI), Eigen::Matrix<double, 1, dim>::Zero(),
            svInds, edgeInds, triInds);

        // triangle-point
        for (const auto& sfI : triInds) {
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI) {
#endif
            const RowVector3i& sfVInd = mesh.SF.row(sfI);

            double d_sqrt;
            computePointTriD(
                Base::V.row(vI),
                mesh.V.row(sfVInd[0]),
                mesh.V.row(sfVInd[1]),
                mesh.V.row(sfVInd[2]),
                d_sqrt);
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
                Base::V.row(vI).transpose(),
                mesh.V.row(sfVInd[0]).transpose(),
                mesh.V.row(sfVInd[1]).transpose(),
                mesh.V.row(sfVInd[2]).transpose(),
                Base::V.row(vI).transpose(),
                mesh.V.row(sfVInd[0]).transpose() + searchDir.template segment<dim>(sfVInd[0] * dim),
                mesh.V.row(sfVInd[1]).transpose() + searchDir.template segment<dim>(sfVInd[1] * dim),
                mesh.V.row(sfVInd[2]).transpose() + searchDir.template segment<dim>(sfVInd[2] * dim),
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
                    Base::V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose(),
                    mesh.V.row(sfVInd[1]).transpose(),
                    mesh.V.row(sfVInd[2]).transpose(),
                    Base::V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose() + searchDir.template segment<dim>(sfVInd[0] * dim),
                    mesh.V.row(sfVInd[1]).transpose() + searchDir.template segment<dim>(sfVInd[1] * dim),
                    mesh.V.row(sfVInd[2]).transpose() + searchDir.template segment<dim>(sfVInd[2] * dim),
                    /*err=*/tight_inclusion_vf_err,
                    /*ms=*/0,
                    toi,
                    tolerance,
                    /*max_t=*/stepSize,
                    /* max_itr=*/-1,
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

// edge-edge
#ifdef USE_TBB
    tbb::parallel_for(0, (int)edges.size(), 1, [&](int eJ) {
#else
    for (int eJ = 0; eJ < edges.size(); ++eJ) {
#endif
        const auto& meshEJ = edges[eJ];
#ifdef USE_SH_LFSS
        std::vector<int> svInds, edgeInds; // NOTE: different constraint order will result in numerically different results
        sh.queryEdgeForPE(Base::V.row(meshEJ.first), Base::V.row(meshEJ.second),
            svInds, edgeInds);

        // edge-edge
        for (const auto& eI : edgeInds) {
            const auto& meshEI = mesh.SFEdges[eI];
#else
        for (const auto& meshEI : mesh.SFEdges) {
#endif

            double d_sqrt;
            computeEdgeEdgeD(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), d_sqrt);
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
                mesh.V.row(meshEI.first).transpose(),
                mesh.V.row(meshEI.second).transpose(),
                Base::V.row(meshEJ.first).transpose(),
                Base::V.row(meshEJ.second).transpose(),
                mesh.V.row(meshEI.first).transpose() + searchDir.template segment<dim>(meshEI.first * dim),
                mesh.V.row(meshEI.second).transpose() + searchDir.template segment<dim>(meshEI.second * dim),
                Base::V.row(meshEJ.first).transpose(),
                Base::V.row(meshEJ.second).transpose(),
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
                has_collision = inclusion_ccd::vertexFaceCCD_double(
                    mesh.V.row(meshEI.first).transpose(),
                    mesh.V.row(meshEI.second).transpose(),
                    Base::V.row(meshEJ.first).transpose(),
                    Base::V.row(meshEJ.second).transpose(),
                    mesh.V.row(meshEI.first).transpose() + searchDir.template segment<dim>(meshEI.first * dim),
                    mesh.V.row(meshEI.second).transpose() + searchDir.template segment<dim>(meshEI.second * dim),
                    Base::V.row(meshEJ.first).transpose(),
                    Base::V.row(meshEJ.second).transpose(),
                    /*err=*/tight_inclusion_ee_err,
                    /*ms=*/0,
                    toi,
                    tolerance,
                    /*max_t=*/stepSize,
                    /* max_itr=*/-1,
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
}
#endif // IPC_WITH_TIGHT_INCLUSION

template <int dim>
void MeshCO<dim>::largestFeasibleStepSize_CCD_exact(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh, const Eigen::VectorXd& searchDir,
    const ccd::CCDMethod ccdMethod, double& stepSize)
{
    // point-triangle
    Eigen::VectorXd largestAlphaPT(Base::F.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::F.rows(), 1, [&](int sfI)
#else
    for (int sfI = 0; sfI < Base::F.rows(); ++sfI)
#endif
        {
            largestAlphaPT[sfI] = stepSize;
            const RowVector3i& sfVInd = Base::F.row(sfI);
#ifdef USE_SH_LFSS
            std::unordered_set<int> pointInds; // NOTE: different constraint order will result in numerically different results
            sh.queryTriangleForPoints(Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]),
                0.0, pointInds);
            for (const auto& svI : pointInds)
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
            {
                int vI = mesh.SVI[svI];

                while (ccd::vertexFaceCCD(
                    mesh.V.row(vI).transpose(),
                    Base::V.row(sfVInd[0]).transpose(),
                    Base::V.row(sfVInd[1]).transpose(),
                    Base::V.row(sfVInd[2]).transpose(),
                    mesh.V.row(vI).transpose() + largestAlphaPT[sfI] * searchDir.template segment<dim>(vI * dim),
                    Base::V.row(sfVInd[0]).transpose(),
                    Base::V.row(sfVInd[1]).transpose(),
                    Base::V.row(sfVInd[2]).transpose(),
                    ccdMethod)) {
                    largestAlphaPT[sfI] /= 2.0;
                }
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphaPT.minCoeff());

    // triangle-point
    Eigen::VectorXd largestAlphaTP(Base::V.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::V.rows(), 1, [&](int vI)
#else
    for (int vI = 0; vI < Base::V.rows(); ++vI)
#endif
        {
            largestAlphaTP[vI] = stepSize;
#ifdef USE_SH_LFSS
            std::unordered_set<int> triInds; // NOTE: different constraint order will result in numerically different results
            sh.queryPointForTriangles(Base::V.row(vI), 0.0, triInds);

            // triangle-point
            for (const auto& sfI : triInds)
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
            {
                const RowVector3i& sfVInd = mesh.SF.row(sfI);

                while (ccd::vertexFaceCCD(
                    Base::V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose(),
                    mesh.V.row(sfVInd[1]).transpose(),
                    mesh.V.row(sfVInd[2]).transpose(),
                    Base::V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose() + largestAlphaTP[vI] * searchDir.template segment<dim>(sfVInd[0] * dim),
                    mesh.V.row(sfVInd[1]).transpose() + largestAlphaTP[vI] * searchDir.template segment<dim>(sfVInd[1] * dim),
                    mesh.V.row(sfVInd[2]).transpose() + largestAlphaTP[vI] * searchDir.template segment<dim>(sfVInd[2] * dim),
                    ccdMethod)) {
                    largestAlphaTP[vI] /= 2.0;
                }
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphaTP.minCoeff());

    // edge-edge
    Eigen::VectorXd largestAlphaEE(edges.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)edges.size(), 1, [&](int eJ)
#else
    for (int eJ = 0; eJ < edges.size(); ++eJ)
#endif
        {
            largestAlphaEE[eJ] = stepSize;
            const auto& meshEJ = edges[eJ];
#ifdef USE_SH_LFSS
            std::vector<int> edgeInds; // NOTE: different constraint order will result in numerically different results
            sh.queryEdgeForEdges(Base::V.row(meshEJ.first), Base::V.row(meshEJ.second),
                0.0, edgeInds);

            // edge-edge
            for (const auto& eI : edgeInds) {
                const auto& meshEI = mesh.SFEdges[eI];
#else
        for (const auto& meshEI : mesh.SFEdges) {
#endif
                while (ccd::edgeEdgeCCD(
                    mesh.V.row(meshEI.first).transpose(),
                    mesh.V.row(meshEI.second).transpose(),
                    Base::V.row(meshEJ.first).transpose(),
                    Base::V.row(meshEJ.second).transpose(),
                    mesh.V.row(meshEI.first).transpose() + largestAlphaEE[eJ] * searchDir.template segment<dim>(meshEI.first * dim),
                    mesh.V.row(meshEI.second).transpose() + largestAlphaEE[eJ] * searchDir.template segment<dim>(meshEI.second * dim),
                    Base::V.row(meshEJ.first).transpose(),
                    Base::V.row(meshEJ.second).transpose(),
                    ccdMethod)) {
                    largestAlphaEE[eJ] /= 2.0;
                }
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphaEE.minCoeff());
}

template <int dim>
void MeshCO<dim>::computeConstraintSet(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh, double dHat,
    std::vector<MMCVID>& constraintSet,
    std::vector<MMCVID>& paraEEMMCVIDSet,
    std::vector<std::pair<int, int>>& paraEEeIeJSet,
    bool getPTEE, std::vector<std::pair<int, int>>& cs_PTEE) const
{
    // point-triangle
    timer_temp3.start(8);
    std::vector<std::vector<MMCVID>> constraintSetPT(Base::F.rows());
    std::vector<std::vector<int>> cs_PT;
    if (getPTEE) {
        cs_PT.resize(Base::F.rows());
    }
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::F.rows(), 1, [&](int sfI)
#else
    for (int sfI = 0; sfI < Base::F.rows(); ++sfI)
#endif
        {
            double dHatFI = dHat;
            const RowVector3i& sfVInd = Base::F.row(sfI);
#ifdef USE_SH_CCS
            std::unordered_set<int> pointInds; // NOTE: different constraint order will result in numerically different results
            sh.queryTriangleForPoints(Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]),
                std::sqrt(dHatFI), pointInds);
            for (const auto& svI : pointInds)
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
            {
                int vI = mesh.SVI[svI];

                int dtype = dType_PT(mesh.V.row(vI), Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]));
                double d;
                switch (dtype) {
                case 0: {
                    d_PP(mesh.V.row(vI), Base::V.row(sfVInd[0]), d);
                    if (d < dHatFI) {
                        constraintSetPT[sfI].emplace_back(-vI - 1, sfVInd[0], -1, -1);
                    }
                    break;
                }

                case 1: {
                    d_PP(mesh.V.row(vI), Base::V.row(sfVInd[1]), d);
                    if (d < dHatFI) {
                        constraintSetPT[sfI].emplace_back(-vI - 1, sfVInd[1], -1, -1);
                    }
                    break;
                }

                case 2: {
                    d_PP(mesh.V.row(vI), Base::V.row(sfVInd[2]), d);
                    if (d < dHatFI) {
                        constraintSetPT[sfI].emplace_back(-vI - 1, sfVInd[2], -1, -1);
                    }
                    break;
                }

                case 3: {
                    d_PE(mesh.V.row(vI), Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), d);
                    if (d < dHatFI) {
                        constraintSetPT[sfI].emplace_back(-vI - 1, sfVInd[0], sfVInd[1], -1);
                    }
                    break;
                }

                case 4: {
                    d_PE(mesh.V.row(vI), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]), d);
                    if (d < dHatFI) {
                        constraintSetPT[sfI].emplace_back(-vI - 1, sfVInd[1], sfVInd[2], -1);
                    }
                    break;
                }

                case 5: {
                    d_PE(mesh.V.row(vI), Base::V.row(sfVInd[2]), Base::V.row(sfVInd[0]), d);
                    if (d < dHatFI) {
                        constraintSetPT[sfI].emplace_back(-vI - 1, sfVInd[2], sfVInd[0], -1);
                    }
                    break;
                }

                case 6: {
                    d_PT(mesh.V.row(vI), Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]), d);
                    if (d < dHatFI) {
                        constraintSetPT[sfI].emplace_back(-vI - 1, sfVInd[0], sfVInd[1], sfVInd[2]);
                    }
                    break;
                }

                default:
                    break;
                }

                if (getPTEE && d < dHatFI) {
                    cs_PT[sfI].emplace_back(svI);
                }
            }
        }
#ifdef USE_TBB
    );
#endif

    // triangle-point
    std::vector<std::vector<MMCVID>> constraintSetTP(Base::V.rows());
    std::vector<std::vector<int>> cs_TP;
    if (getPTEE) {
        cs_TP.resize(Base::V.rows());
    }
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::V.rows(), 1, [&](int vI)
#else
    for (int vI = 0; vI < Base::V.rows(); ++vI)
#endif
        {
            double dHatVI = dHat;
#ifdef USE_SH_CCS
            std::unordered_set<int> triInds; // NOTE: different constraint order will result in numerically different results
            sh.queryPointForTriangles(Base::V.row(vI), std::sqrt(dHatVI), triInds);
            for (const auto& sfI : triInds)
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
            {
                const RowVector3i& sfVInd = mesh.SF.row(sfI);

                int dtype = dType_PT(Base::V.row(vI), mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]));
                double d;
                switch (dtype) {
                case 0: {
                    d_PP(Base::V.row(vI), mesh.V.row(sfVInd[0]), d);
                    if (d < dHatVI) {
                        constraintSetTP[vI].emplace_back(-sfVInd[0] - 1, vI, -1, -1);
                    }
                    break;
                }

                case 1: {
                    d_PP(Base::V.row(vI), mesh.V.row(sfVInd[1]), d);
                    if (d < dHatVI) {
                        constraintSetTP[vI].emplace_back(-sfVInd[1] - 1, vI, -1, -1);
                    }
                    break;
                }

                case 2: {
                    d_PP(Base::V.row(vI), mesh.V.row(sfVInd[2]), d);
                    if (d < dHatVI) {
                        constraintSetTP[vI].emplace_back(-sfVInd[2] - 1, vI, -1, -1);
                    }
                    break;
                }

                case 3: {
                    d_PE(Base::V.row(vI), mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), d);
                    if (d < dHatVI) {
                        constraintSetTP[vI].emplace_back(-sfVInd[0] - 1, -sfVInd[1] - 1, vI, -1);
                    }
                    break;
                }

                case 4: {
                    d_PE(Base::V.row(vI), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), d);
                    if (d < dHatVI) {
                        constraintSetTP[vI].emplace_back(-sfVInd[1] - 1, -sfVInd[2] - 1, vI, -1);
                    }
                    break;
                }

                case 5: {
                    d_PE(Base::V.row(vI), mesh.V.row(sfVInd[2]), mesh.V.row(sfVInd[0]), d);
                    if (d < dHatVI) {
                        constraintSetTP[vI].emplace_back(-sfVInd[2] - 1, -sfVInd[0] - 1, vI, -1);
                    }
                    break;
                }

                case 6: {
                    d_PT(Base::V.row(vI), mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), d);
                    if (d < dHatVI) {
                        constraintSetTP[vI].emplace_back(-sfVInd[0] - 1, -sfVInd[1] - 1, -sfVInd[2] - 1, vI);
                    }
                    break;
                }

                default:
                    break;
                }

                if (getPTEE && d < dHatVI) {
                    cs_TP[vI].emplace_back(sfI);
                }
            }
        }
#ifdef USE_TBB
    );
#endif
    timer_temp3.stop();

    // edge-edge
    timer_temp3.start(9);
    std::vector<std::vector<MMCVID>> constraintSetEE(edges.size());
    std::vector<std::vector<int>> cs_EE;
    if (getPTEE) {
        cs_EE.resize(edges.size());
    }
#ifdef USE_TBB
    tbb::parallel_for(0, (int)edges.size(), 1, [&](int eJ)
#else
    for (int eJ = 0; eJ < edges.size(); ++eJ)
#endif
        {
            double dHatEJ = dHat;
            timer_mt.start(6);
            const auto& meshEJ = edges[eJ];
            timer_mt.stop();
#ifdef USE_SH_CCS
            std::vector<int> edgeInds; // NOTE: different constraint order will result in numerically different results
            // timer_mt.start(0);
            sh.queryEdgeForEdges(Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), std::sqrt(dHatEJ), edgeInds);
            // timer_mt.stop();
            timer_mt.start(23);
            for (const auto& eI : edgeInds) {
#else
        for (int eI = 0; eI < mesh.SFEdges.size(); eI++) {
#endif
                timer_mt.start(6);
                const auto& meshEI = mesh.SFEdges[eI];
                timer_mt.stop();
                timer_mt.start(1);
                int dtype = dType_EE(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), Base::V.row(meshEJ.first), Base::V.row(meshEJ.second));
                timer_mt.stop();

                timer_mt.start(22);
                double EECrossSqNorm, eps_x;
                computeEECrossSqNorm(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), EECrossSqNorm);
                compute_eps_x(mesh, Base::V, meshEI.first, meshEI.second, meshEJ.first, meshEJ.second, eps_x);
                // int add_e = -1;
                int add_e = (EECrossSqNorm < eps_x) ? -eI - 2 : -1;
                // == -1: regular,
                // <= -2 && >= -mesh.SFEdges.size()-1: nearly parallel PP or PE,
                // <= -mesh.SFEdges.size()-2: nearly parallel EE

                double d;
                timer_mt.start(2);
                switch (dtype) {
                case 0: {
                    d_PP(mesh.V.row(meshEI.first), Base::V.row(meshEJ.first), d);
                    if (d < dHatEJ) {
                        constraintSetEE[eJ].emplace_back(-meshEI.first - 1, meshEJ.first, -1, add_e);
                    }
                    break;
                }

                case 1: {
                    d_PP(mesh.V.row(meshEI.first), Base::V.row(meshEJ.second), d);
                    if (d < dHatEJ) {
                        constraintSetEE[eJ].emplace_back(-meshEI.first - 1, meshEJ.second, -1, add_e);
                    }
                    break;
                }

                case 2: {
                    d_PE(mesh.V.row(meshEI.first), Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), d);
                    if (d < dHatEJ) {
                        constraintSetEE[eJ].emplace_back(-meshEI.first - 1, meshEJ.first, meshEJ.second, add_e);
                    }
                    break;
                }

                case 3: {
                    d_PP(mesh.V.row(meshEI.second), Base::V.row(meshEJ.first), d);
                    if (d < dHatEJ) {
                        constraintSetEE[eJ].emplace_back(-meshEI.second - 1, meshEJ.first, -1, add_e);
                    }
                    break;
                }

                case 4: {
                    d_PP(mesh.V.row(meshEI.second), Base::V.row(meshEJ.second), d);
                    if (d < dHatEJ) {
                        constraintSetEE[eJ].emplace_back(-meshEI.second - 1, meshEJ.second, -1, add_e);
                    }
                    break;
                }

                case 5: {
                    d_PE(mesh.V.row(meshEI.second), Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), d);
                    if (d < dHatEJ) {
                        constraintSetEE[eJ].emplace_back(-meshEI.second - 1, meshEJ.first, meshEJ.second, add_e);
                    }
                    break;
                }

                case 6: {
                    d_PE(Base::V.row(meshEJ.first), mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), d);
                    if (d < dHatEJ) {
                        constraintSetEE[eJ].emplace_back(-meshEI.first - 1, -meshEI.second - 1, meshEJ.first, add_e);
                    }
                    break;
                }

                case 7: {
                    d_PE(Base::V.row(meshEJ.second), mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), d);
                    if (d < dHatEJ) {
                        constraintSetEE[eJ].emplace_back(-meshEI.first - 1, -meshEI.second - 1, meshEJ.second, add_e);
                    }
                    break;
                }

                case 8: {
                    d_EE(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second), Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), d);
                    if (d < dHatEJ) {
                        if (add_e <= -2) {
                            constraintSetEE[eJ].emplace_back(meshEI.first, meshEI.second, meshEJ.first, -meshEJ.second - mesh.SFEdges.size() - 2);
                        }
                        else {
                            constraintSetEE[eJ].emplace_back(meshEI.first, meshEI.second, meshEJ.first, meshEJ.second);
                        }
                    }
                    break;
                }

                default:
                    break;
                }

                if (getPTEE && d < dHatEJ) {
                    cs_EE[eJ].emplace_back(eI);
                }
                timer_mt.stop();
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
        cs_PTEE.reserve(cs_PT.size() + cs_TP.size() + cs_EE.size());
        for (int sfI = 0; sfI < cs_PT.size(); ++sfI) {
            for (const auto& svI : cs_PT[sfI]) {
                cs_PTEE.emplace_back(-svI - 1, sfI);
            }
        }
        for (int vI = 0; vI < cs_TP.size(); ++vI) {
            for (const auto& sfI : cs_TP[vI]) {
                cs_PTEE.emplace_back(-sfI - 1, -vI - 1);
            }
        }
        for (int eI = 0; eI < cs_EE.size(); ++eI) {
            for (const auto& eJ : cs_EE[eI]) {
                cs_PTEE.emplace_back(eJ, eI);
            }
        }
    }

    constraintSet.resize(0);
    constraintSet.reserve(constraintSetPT.size() + constraintSetTP.size() + constraintSetEE.size());
    // for (const auto& csI : constraintSetPT) {
    //     constraintSet.insert(constraintSet.end(), csI.begin(), csI.end());
    // }
    // for (const auto& csI : constraintSetTP) {
    //     constraintSet.insert(constraintSet.end(), csI.begin(), csI.end());
    // }
    // for (const auto& csI : constraintSetEE) {
    //     constraintSet.insert(constraintSet.end(), csI.begin(), csI.end());
    // }
    std::map<MMCVID, int> constraintCounter;
    for (const auto& csI : constraintSetPT) {
        for (const auto& cI : csI) {
            if (cI[3] < 0) {
                // PP or PE or EP
                ++constraintCounter[cI];
            }
            else {
                constraintSet.emplace_back(cI);
            }
        }
    }
    for (const auto& csI : constraintSetTP) {
        for (const auto& cI : csI) {
            if (cI[3] < 0) {
                // PP or PE or EP
                ++constraintCounter[cI];
            }
            else {
                constraintSet.emplace_back(cI);
            }
        }
    }
    paraEEMMCVIDSet.resize(0);
    paraEEeIeJSet.resize(0);
    int eJ = 0;
    for (const auto& csI : constraintSetEE) {
        for (const auto& cI : csI) {
            if (cI[3] >= 0) {
                // regular EE
                constraintSet.emplace_back(cI);
            }
            else if (cI[3] == -1) {
                // regular PP or PE or EP
                ++constraintCounter[cI];
            }
            else if (cI[3] >= -mesh.SFEdges.size() - 1) {
                // nearly parallel PP or PE or EP
                paraEEMMCVIDSet.emplace_back(cI[0], cI[1], cI[2], -1);
                paraEEeIeJSet.emplace_back(-cI[3] - 2, eJ);
            }
            else {
                // nearly parallel EE
                paraEEMMCVIDSet.emplace_back(cI[0], cI[1], cI[2], -cI[3] - mesh.SFEdges.size() - 2);
                paraEEeIeJSet.emplace_back(-1, -1);
            }
        }
        ++eJ;
    }

    constraintSet.reserve(constraintSet.size() + constraintCounter.size());
    for (const auto& ccI : constraintCounter) {
        constraintSet.emplace_back(MMCVID(ccI.first[0], ccI.first[1], ccI.first[2], -ccI.second));
    }
    timer_temp3.stop();
}

template <int dim>
void MeshCO<dim>::augmentParaEEEnergy(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& paraEEMMCVIDSet,
    const std::vector<std::pair<int, int>>& paraEEeIeJSet,
    Eigen::VectorXd& d_inc, Eigen::VectorXd& energyVec_inc,
    double dHat, double coef) const
{
    int startCI = d_inc.size();
    Base::evaluateConstraints(mesh, paraEEMMCVIDSet, d_inc);

    energyVec_inc.conservativeResize(d_inc.size());
    for (int cI = startCI; cI < energyVec_inc.size(); ++cI) {
        if (d_inc[cI] <= 0.0) {
            std::cout << "constraintVal = " << d_inc[cI] << " when evaluating constraint ";
            exit(0);
        }
        else {
            const MMCVID& MMCVIDI = paraEEMMCVIDSet[cI - startCI];
            double eps_x, e;
            if (MMCVIDI[3] >= 0) {
                // EE
                compute_eps_x(mesh, Base::V, MMCVIDI[0], MMCVIDI[1], MMCVIDI[2], MMCVIDI[3], eps_x);
                compute_e(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]),
                    Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), eps_x, e);
            }
            else {
                // PP or PE
                const std::pair<int, int>& eIeJ = paraEEeIeJSet[cI - startCI];
                const std::pair<int, int>& eI = mesh.SFEdges[eIeJ.first];
                const std::pair<int, int>& eJ = edges[eIeJ.second];
                compute_eps_x(mesh, Base::V, eI.first, eI.second, eJ.first, eJ.second, eps_x);
                compute_e(mesh.V.row(eI.first), mesh.V.row(eI.second),
                    Base::V.row(eJ.first), Base::V.row(eJ.second), eps_x, e);
            }
            compute_b(d_inc[cI], dHat, energyVec_inc[cI]);
            energyVec_inc[cI] *= coef * e;
        }
    }
}

template <int dim>
void MeshCO<dim>::augmentParaEEGradient(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& paraEEMMCVIDSet,
    const std::vector<std::pair<int, int>>& paraEEeIeJSet,
    Eigen::VectorXd& grad_inc, double dHat, double coef) const
{
    Eigen::VectorXd e_db_div_dd;
    Base::evaluateConstraints(mesh, paraEEMMCVIDSet, e_db_div_dd);

    for (int cI = 0; cI < e_db_div_dd.size(); ++cI) {
        const MMCVID& MMCVIDI = paraEEMMCVIDSet[cI];

        double b;
        compute_b(e_db_div_dd[cI], dHat, b);
        compute_g_b(e_db_div_dd[cI], dHat, e_db_div_dd[cI]);

        double eps_x, e;
        double coef_b = coef * b;
        if (MMCVIDI[3] >= 0) {
            // EE
            compute_eps_x(mesh, Base::V, MMCVIDI[0], MMCVIDI[1], MMCVIDI[2], MMCVIDI[3], eps_x);
            compute_e(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), eps_x, e);

            Eigen::Matrix<double, 12, 1> e_g;
            compute_e_g(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), eps_x, e_g);

            grad_inc.template segment<dim>(MMCVIDI[0] * dim) += coef_b * e_g.template segment<dim>(0);
            grad_inc.template segment<dim>(MMCVIDI[1] * dim) += coef_b * e_g.template segment<dim>(dim);
        }
        else {
            // PP or PE
            const std::pair<int, int>& eIeJ = paraEEeIeJSet[cI];
            const std::pair<int, int>& eI = mesh.SFEdges[eIeJ.first];
            const std::pair<int, int>& eJ = edges[eIeJ.second];
            compute_eps_x(mesh, Base::V, eI.first, eI.second, eJ.first, eJ.second, eps_x);
            compute_e(mesh.V.row(eI.first), mesh.V.row(eI.second), Base::V.row(eJ.first), Base::V.row(eJ.second), eps_x, e);

            Eigen::Matrix<double, 12, 1> e_g;
            compute_e_g(mesh.V.row(eI.first), mesh.V.row(eI.second), Base::V.row(eJ.first), Base::V.row(eJ.second), eps_x, e_g);

            grad_inc.template segment<dim>(eI.first * dim) += coef_b * e_g.template segment<dim>(0);
            grad_inc.template segment<dim>(eI.second * dim) += coef_b * e_g.template segment<dim>(dim);
        }
        e_db_div_dd[cI] *= e;
    }
    leftMultiplyConstraintJacobianT(mesh, paraEEMMCVIDSet, e_db_div_dd, grad_inc, coef);
}

template <int dim>
void MeshCO<dim>::augmentParaEEHessian(const Mesh<dim>& mesh,
    const std::vector<MMCVID>& paraEEMMCVIDSet,
    const std::vector<std::pair<int, int>>& paraEEeIeJSet,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* H_inc,
    double dHat, double coef, bool projectDBC) const
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

                    compute_eps_x(mesh, Base::V, MMCVIDI[0], MMCVIDI[1], MMCVIDI[2], MMCVIDI[3], eps_x);
                    compute_e(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), eps_x, e);
                    compute_e_g(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), eps_x, e_g);
                    compute_e_H(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), eps_x, e_H);

                    g_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), grad_d);
                    H_EE(mesh.V.row(MMCVIDI[0]), mesh.V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), Base::V.row(MMCVIDI[3]), H_d);

                    rowIStart[cI][0] = mesh.isProjectDBCVertex(MMCVIDI[0], projectDBC) ? -1 : (MMCVIDI[0] * dim);
                    rowIStart[cI][1] = mesh.isProjectDBCVertex(MMCVIDI[1], projectDBC) ? -1 : (MMCVIDI[1] * dim);
                    rowIStart[cI][2] = -1;
                    rowIStart[cI][3] = -1;
                }
                else {
                    // PP or PE
                    compute_b(d, dHat, b);
                    compute_g_b(d, dHat, g_b);
                    compute_H_b(d, dHat, H_b);

                    const std::pair<int, int>& eIeJ = paraEEeIeJSet[cI];
                    const std::pair<int, int>& eI = mesh.SFEdges[eIeJ.first];
                    const std::pair<int, int>& eJ = edges[eIeJ.second];
                    compute_eps_x(mesh, Base::V, eI.first, eI.second, eJ.first, eJ.second, eps_x);
                    compute_e(mesh.V.row(eI.first), mesh.V.row(eI.second), Base::V.row(eJ.first), Base::V.row(eJ.second), eps_x, e);
                    compute_e_g(mesh.V.row(eI.first), mesh.V.row(eI.second), Base::V.row(eJ.first), Base::V.row(eJ.second), eps_x, e_g);
                    compute_e_H(mesh.V.row(eI.first), mesh.V.row(eI.second), Base::V.row(eJ.first), Base::V.row(eJ.second), eps_x, e_H);

                    rowIStart[cI][0] = mesh.isProjectDBCVertex(eI.first, projectDBC) ? -1 : (eI.first * dim);
                    rowIStart[cI][1] = mesh.isProjectDBCVertex(eI.second, projectDBC) ? -1 : (eI.second * dim);
                    rowIStart[cI][2] = -1;
                    rowIStart[cI][3] = -1;

                    int v0I = -MMCVIDI[0] - 1;
                    if (MMCVIDI[2] >= 0) {
                        if (MMCVIDI[1] >= 0) {
                            // PE
                            Eigen::Matrix<double, 9, 1> grad_d_PE;
                            g_PE(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), grad_d_PE);
                            Eigen::Matrix<double, 9, 9> H_d_PE;
                            H_PE(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), Base::V.row(MMCVIDI[2]), H_d_PE);

                            // fill in
                            int ind[4] = { eI.first, eI.second, eJ.first, eJ.second };
                            int indMap[3];
                            for (int i = 0; i < 2; ++i) {
                                if (v0I == ind[i]) {
                                    indMap[0] = i;
                                    break;
                                }
                            }
                            for (int i = 2; i < 4; ++i) {
                                if (MMCVIDI[1] == ind[i]) {
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
                            // EP
                            int v1I = -MMCVIDI[1] - 1;
                            Eigen::Matrix<double, 9, 1> grad_d_PE;
                            g_PE(Base::V.row(MMCVIDI[2]), mesh.V.row(v0I), mesh.V.row(v1I), grad_d_PE);
                            Eigen::Matrix<double, 9, 9> H_d_PE;
                            H_PE(Base::V.row(MMCVIDI[2]), mesh.V.row(v0I), mesh.V.row(v1I), H_d_PE);

                            // fill in
                            int ind[4] = { eI.first, eI.second, eJ.first, eJ.second };
                            int indMap[3];
                            for (int i = 0; i < 2; ++i) {
                                if (v0I == ind[i]) {
                                    indMap[1] = i;
                                }
                                else if (v1I == ind[i]) {
                                    indMap[2] = i;
                                }
                            }
                            for (int i = 2; i < 4; ++i) {
                                if (MMCVIDI[2] == ind[i]) {
                                    indMap[0] = i;
                                    break;
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
                    }
                    else {
                        // PP
                        Eigen::Matrix<double, 6, 1> grad_d_PP;
                        g_PP(mesh.V.row(v0I), Base::V.row(MMCVIDI[1]), grad_d_PP);
                        Eigen::Matrix<double, 6, 6> H_d_PP;
                        H_PP(H_d_PP);

                        int ind[4] = { eI.first, eI.second, eJ.first, eJ.second };
                        int indMap[2];
                        for (int i = 0; i < 2; ++i) {
                            if (v0I == ind[i]) {
                                indMap[0] = i;
                                break;
                            }
                        }
                        for (int i = 2; i < 4; ++i) {
                            if (MMCVIDI[1] == ind[i]) {
                                indMap[1] = i;
                                break;
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
bool MeshCO<dim>::checkEdgeTriIntersection(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh)
{
    std::vector<std::vector<int>> intersectedEdges(Base::F.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::F.rows(), 1, [&](int cofI)
#else
    for (int cofI = 0; cofI < Base::F.rows(); ++cofI)
#endif
        {
            const Eigen::RowVector3i& cofVInd = Base::F.row(cofI);
#ifdef USE_SH_CCS
            std::unordered_set<int> sEdgeInds;
            sh.queryTriangleForEdges(Base::V.row(cofVInd[0]), Base::V.row(cofVInd[1]), Base::V.row(cofVInd[2]), 0.0, sEdgeInds);
            for (const auto& eI : sEdgeInds)
#else
        for (int eI = 0; eI < mesh.SFEdges.size(); ++eI)
#endif
            {
                const auto& meshEI = mesh.SFEdges[eI];
                if (IglUtils::segTriIntersect(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                        Base::V.row(cofVInd[0]), Base::V.row(cofVInd[1]), Base::V.row(cofVInd[2]))) {
                    intersectedEdges[cofI].emplace_back(eI);
                }
            }
        }
#ifdef USE_TBB
    );
#endif

    std::vector<std::vector<int>> intersectedTris(edges.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)edges.size(), 1, [&](int eI)
#else
    for (int eI = 0; eI < edges.size(); ++eI)
#endif
        {
            const auto& coEI = edges[eI];
#ifdef USE_SH_CCS
            std::unordered_set<int> sTriInds;
            sh.queryEdgeForTriangles(Base::V.row(coEI.first), Base::V.row(coEI.second), 0.0, sTriInds);
            for (const auto& sfI : sTriInds)
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
            {
                const Eigen::RowVector3i& sfVInd = mesh.SF.row(sfI);
                if (IglUtils::segTriIntersect(Base::V.row(coEI.first), Base::V.row(coEI.second),
                        mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]))) {
                    intersectedTris[eI].emplace_back(sfI);
                }
            }
        }
#ifdef USE_TBB
    );
#endif

    bool intersected = false;
    for (int cofI = 0; cofI < intersectedEdges.size(); ++cofI) {
        for (const auto& eI : intersectedEdges[cofI]) {
            intersected = true;
            std::cout << "elastic edge - MCO triangle intersection detected" << std::endl;
            std::cout << mesh.SFEdges[eI].first << " " << mesh.SFEdges[eI].second << std::endl;
            std::cout << Base::F.row(cofI) << std::endl;
            std::cout << mesh.V.row(mesh.SFEdges[eI].first) << std::endl;
            std::cout << mesh.V.row(mesh.SFEdges[eI].second) << std::endl;
            std::cout << Base::V.row(Base::F(cofI, 0)) << std::endl;
            std::cout << Base::V.row(Base::F(cofI, 1)) << std::endl;
            std::cout << Base::V.row(Base::F(cofI, 2)) << std::endl;
        }
    }
    for (int eI = 0; eI < intersectedTris.size(); ++eI) {
        for (const auto& sfI : intersectedTris[eI]) {
            intersected = true;
            std::cout << "MCO edge - elastic triangle intersection detected" << std::endl;
            std::cout << edges[eI].first << " " << edges[eI].second << std::endl;
            std::cout << mesh.SF.row(sfI) << std::endl;
            std::cout << Base::V.row(edges[eI].first) << std::endl;
            std::cout << Base::V.row(edges[eI].second) << std::endl;
            std::cout << mesh.V.row(mesh.SF(sfI, 0)) << std::endl;
            std::cout << mesh.V.row(mesh.SF(sfI, 1)) << std::endl;
            std::cout << mesh.V.row(mesh.SF(sfI, 2)) << std::endl;
        }
    }
    return !intersected;
}

template <int dim>
bool MeshCO<dim>::checkEdgeTriIntersectionIfAny(const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh)
{
    Eigen::ArrayXi intersected(Base::F.rows());
    intersected.setZero();
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::F.rows(), 1, [&](int cofI)
#else
    for (int cofI = 0; cofI < Base::F.rows(); ++cofI)
#endif
        {
            const Eigen::RowVector3i& cofVInd = Base::F.row(cofI);
#ifdef USE_SH_CCS
            std::unordered_set<int> sEdgeInds;
            sh.queryTriangleForEdges(Base::V.row(cofVInd[0]), Base::V.row(cofVInd[1]), Base::V.row(cofVInd[2]), 0.0, sEdgeInds);
            for (const auto& eI : sEdgeInds)
#else
        for (int eI = 0; eI < mesh.SFEdges.size(); ++eI)
#endif
            {
                const auto& meshEI = mesh.SFEdges[eI];
                if (IglUtils::segTriIntersect(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                        Base::V.row(cofVInd[0]), Base::V.row(cofVInd[1]), Base::V.row(cofVInd[2]))) {
                    intersected[cofI] = 1;
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

    intersected.setZero(edges.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)edges.size(), 1, [&](int eI)
#else
    for (int eI = 0; eI < edges.size(); ++eI)
#endif
        {
            const auto& coEI = edges[eI];
#ifdef USE_SH_CCS
            std::unordered_set<int> sTriInds;
            sh.queryEdgeForTriangles(Base::V.row(coEI.first), Base::V.row(coEI.second), 0.0, sTriInds);
            for (const auto& sfI : sTriInds)
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
            {
                const Eigen::RowVector3i& sfVInd = mesh.SF.row(sfI);
                if (IglUtils::segTriIntersect(Base::V.row(coEI.first), Base::V.row(coEI.second),
                        mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]))) {
                    intersected[eI] = 1;
                    break;
                }
            }
        }
#ifdef USE_TBB
    );
#endif

    return (intersected == 0).all();
}

template <int dim>
void MeshCO<dim>::updateConstraints_QP(
    const Mesh<dim>& mesh,
    const std::vector<MMCVID>& activeSet,
    const CollisionConstraintType constraintType,
    std::vector<Eigen::Triplet<double>>& A_triplet,
    Eigen::VectorXd& l) const
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
                spdlog::error("where=\"MeshCO:evaluateConstraintsQP\" msg=\"Cannot compute the constraint value of MMCVID without toi!\"");
                spdlog::debug("MMCVIDI={:s}", MMCVIDI.str());
                spdlog::debug("mmcvid_to_toi:");
                for (auto const& [mmcvid, toi] : mmcvid_to_toi) {
                    spdlog::debug("{:s}: {:g}", mmcvid.str(), toi);
                }
                throw "missing toi for constraint";
            }

            if (MMCVIDI[0] >= 0) {
                // edge-edge
                compute_collision_constraint_gradient(
                    mesh.V_prev.row(MMCVIDI[0]),
                    mesh.V_prev.row(MMCVIDI[1]),
                    Base::V.row(MMCVIDI[2]),
                    Base::V.row(MMCVIDI[3]),
                    mesh.V.row(MMCVIDI[0]),
                    mesh.V.row(MMCVIDI[1]),
                    Base::V.row(MMCVIDI[2]),
                    Base::V.row(MMCVIDI[3]),
                    constraintType, /*is_edge_edge=*/true, toi,
                    gradc);

                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < dim; j++) {
                        A_triplet.emplace_back(
                            constraintI,
                            MMCVIDI[i] * dim + j,
                            gradc[i * dim + j]);
                    }
                }
            }
            else {
                // point and triangle
                if (MMCVIDI[1] >= 0) {
                    // point-triangle
                    compute_collision_constraint_gradient(
                        mesh.V_prev.row(-MMCVIDI[0] - 1),
                        Base::V.row(MMCVIDI[1]),
                        Base::V.row(MMCVIDI[2]),
                        Base::V.row(MMCVIDI[3]),
                        mesh.V.row(-MMCVIDI[0] - 1),
                        Base::V.row(MMCVIDI[1]),
                        Base::V.row(MMCVIDI[2]),
                        Base::V.row(MMCVIDI[3]),
                        constraintType, /*is_edge_edge=*/false, toi,
                        gradc);

                    for (int j = 0; j < dim; j++) {
                        A_triplet.emplace_back(
                            constraintI,
                            (-MMCVIDI[0] - 1) * dim + j,
                            gradc[j]);
                    }
                }
                else {
                    // triangle-point
                    compute_collision_constraint_gradient(
                        Base::V.row(MMCVIDI[3]),
                        mesh.V.row(-MMCVIDI[0] - 1),
                        mesh.V.row(-MMCVIDI[1] - 1),
                        mesh.V.row(-MMCVIDI[2] - 1),
                        Base::V.row(MMCVIDI[3]),
                        mesh.V_prev.row(-MMCVIDI[0] - 1),
                        mesh.V_prev.row(-MMCVIDI[1] - 1),
                        mesh.V_prev.row(-MMCVIDI[2] - 1),
                        constraintType, /*is_edge_edge=*/false, toi,
                        gradc);

                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < dim; j++) {
                            A_triplet.emplace_back(
                                constraintI,
                                (-MMCVIDI[i] - 1) * dim + j,
                                gradc[(i + 1) * dim + j]);
                        }
                    }
                }
            }

            ++constraintI;
        }

        // Use coef=-1 to move the linear component to the right side of the
        // inequality:  ∇g(x₀) * x + g(x₀) ≥ ϵ ≡ ∇g(x₀) * x ≥ -g(x₀) + ϵ
        this->evaluateConstraintsQP(
            mesh, activeSet, constraintType, l, /*coef=*/-1.0);
    }
    else {
        // TODO: 2D collisions
        spdlog::error("IPC does not handle 2D collisions.");
    }
}

void addValidConstraints(
    const Eigen::VectorXd& tois,
    const std::vector<MMCVID>& mmcvids,
    std::vector<MMCVID>& activeSet)
{
    // TODO: Parallelize
    for (int i = 0; i < tois.size(); i++) {
        if (isfinite(tois[i]) && tois[i] >= 0) {
            activeSet.push_back(mmcvids[i]);
        }
    }
}

// Update the active set using CCD
template <int dim>
bool MeshCO<dim>::updateActiveSet_QP(
    const Mesh<dim>& mesh, // FEM mesh
    const Eigen::VectorXd& searchDir, // Linear displacement of the mesh vertices
    const CollisionConstraintType constraintType, // Update method
    std::vector<MMCVID>& activeSet, // Constraint active set to update
    const ccd::CCDMethod ccdMethod,
    const double eta, // Surface thickness
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

    // Only save the earliest impact
    Eigen::VectorXd min_tois = Eigen::VectorXd::Constant(
        mesh.SVI.size(), std::numeric_limits<double>::infinity());
    std::vector<MMCVID> mmcvids(min_tois.size(), MMCVID());

    // point-triangle
    // Loop over mesh collision object faces
    for (int MCOFI = 0; MCOFI < Base::F.rows(); ++MCOFI) {
        const RowVector3i& MCTriVInd = Base::F.row(MCOFI);

        // Loop over FEM mesh surface vertices
#ifdef USE_SH_CCS
        std::unordered_set<int> pointInds; // NOTE: different constraint order will result in numerically different results
        sh.queryTriangleForPoints(
            this->V.row(MCTriVInd[0]),
            this->V.row(MCTriVInd[1]),
            this->V.row(MCTriVInd[2]),
            eta, pointInds);
        for (const int& svI : pointInds) {
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
#endif
            int vI = mesh.SVI[svI];

            double toi;
            bool intersects;
            switch (ccdMethod) {
            case ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER: {
                intersects = CTCD::vertexFaceCTCD(
                    mesh.V_prev.row(vI).transpose(),
                    Base::V.row(MCTriVInd[0]).transpose(),
                    Base::V.row(MCTriVInd[1]).transpose(),
                    Base::V.row(MCTriVInd[2]).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                    Base::V.row(MCTriVInd[0]).transpose(),
                    Base::V.row(MCTriVInd[1]).transpose(),
                    Base::V.row(MCTriVInd[2]).transpose(),
                    eta,
                    toi);
                break;
            }
            case ccd::CCDMethod::TIGHT_INCLUSION: {
#ifdef IPC_WITH_TIGHT_INCLUSION
                double output_tolerance;
                intersects = inclusion_ccd::vertexFaceCCD_double(
                    mesh.V_prev.row(vI).transpose(),
                    Base::V.row(MCTriVInd[0]).transpose(),
                    Base::V.row(MCTriVInd[1]).transpose(),
                    Base::V.row(MCTriVInd[2]).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                    Base::V.row(MCTriVInd[0]).transpose(),
                    Base::V.row(MCTriVInd[1]).transpose(),
                    Base::V.row(MCTriVInd[2]).transpose(),
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
                    Base::V.row(MCTriVInd[0]).transpose(),
                    Base::V.row(MCTriVInd[1]).transpose(),
                    Base::V.row(MCTriVInd[2]).transpose(),
                    mesh.V.row(vI).transpose() + searchDir.segment<dim>(vI * dim),
                    Base::V.row(MCTriVInd[0]).transpose(),
                    Base::V.row(MCTriVInd[1]).transpose(),
                    Base::V.row(MCTriVInd[2]).transpose(),
                    ccdMethod,
                    toi);
            }

            MMCVID mmcvid = MMCVID(-vI - 1, // mesh point
                MCTriVInd[0], MCTriVInd[2], MCTriVInd[1]); // CO triangle
            mmcvid_to_toi.insert_or_assign(mmcvid, intersects ? toi : std::numeric_limits<double>::infinity());
            bool isNewConstraint = prevActiveSet.find(mmcvid) == prevActiveSet.end();
            if (intersects && toi < min_tois[svI] && (wasActiveSetCleared || isNewConstraint)) {
                newConstraintsAdded |= isNewConstraint;
                min_tois[svI] = toi;
                mmcvids[svI] = mmcvid;
            }
        }
    }

    addValidConstraints(min_tois, mmcvids, activeSet);
    min_tois = Eigen::VectorXd::Constant(
        mesh.SF.rows(), std::numeric_limits<double>::infinity());
    mmcvids = std::vector<MMCVID>(min_tois.size(), MMCVID());

    // triangle-point
    // Loop over mesh collision object vertices
    for (int vI = 0; vI < Base::V.rows(); ++vI) {
        // Loop over FEM mesh surface faces
#ifdef USE_SH_CCS
        std::unordered_set<int> triInds; // NOTE: different constraint order will result in numerically different results
        sh.queryPointForTriangles(Base::V.row(vI), eta, triInds);
        for (const int& sfI : triInds) {
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI) {
#endif
            const RowVector3i& sfVInd = mesh.SF.row(sfI);

            double toi;
            bool intersects;
            switch (ccdMethod) {
            case ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER:
                intersects = CTCD::vertexFaceCTCD(
                    Base::V.row(vI).transpose(),
                    mesh.V_prev.row(sfVInd[0]).transpose(),
                    mesh.V_prev.row(sfVInd[1]).transpose(),
                    mesh.V_prev.row(sfVInd[2]).transpose(),
                    Base::V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                    mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                    mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                    eta,
                    toi);
                break;
            case ccd::CCDMethod::TIGHT_INCLUSION: {
#ifdef IPC_WITH_TIGHT_INCLUSION
                double output_tolerance;
                intersects = inclusion_ccd::vertexFaceCCD_double(
                    Base::V.row(vI).transpose(),
                    mesh.V_prev.row(sfVInd[0]).transpose(),
                    mesh.V_prev.row(sfVInd[1]).transpose(),
                    mesh.V_prev.row(sfVInd[2]).transpose(),
                    Base::V.row(vI).transpose(),
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
                    Base::V.row(vI).transpose(),
                    mesh.V_prev.row(sfVInd[0]).transpose(),
                    mesh.V_prev.row(sfVInd[1]).transpose(),
                    mesh.V_prev.row(sfVInd[2]).transpose(),
                    Base::V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose() + searchDir.segment<dim>(sfVInd[0] * dim),
                    mesh.V.row(sfVInd[1]).transpose() + searchDir.segment<dim>(sfVInd[1] * dim),
                    mesh.V.row(sfVInd[2]).transpose() + searchDir.segment<dim>(sfVInd[2] * dim),
                    ccdMethod,
                    toi);
            }

            MMCVID mmcvid = MMCVID(
                -sfVInd[0] - 1, -sfVInd[2] - 1, -sfVInd[1] - 1, // mesh triangle
                vI); // CO Point
            mmcvid_to_toi.insert_or_assign(mmcvid, intersects ? toi : std::numeric_limits<double>::infinity());

            bool isNewConstraint = prevActiveSet.find(mmcvid) == prevActiveSet.end();
            if (intersects && toi < min_tois[sfI] && (wasActiveSetCleared || isNewConstraint)) {
                newConstraintsAdded |= isNewConstraint;
                min_tois[sfI] = toi;
                mmcvids[sfI] = mmcvid;
            }
        }
    }

    addValidConstraints(min_tois, mmcvids, activeSet);
    min_tois = Eigen::VectorXd::Constant(
        mesh.SFEdges.size(), std::numeric_limits<double>::infinity());
    mmcvids = std::vector<MMCVID>(min_tois.size(), MMCVID());

    // edge-edge
    // Loop over meshCO edges
    // Get the mesh CO edge as a pair of indiced into Base::V
    for (const std::pair<int, int>& co_edge : this->edges) {
        Eigen::Matrix<double, dim, 1> co_edge_dir = (this->V.row(co_edge.second) - this->V.row(co_edge.first)).normalized();

        // Loop over FEM mesh edges
#ifdef USE_SH_CCS
        std::vector<int> sEdgeInds;
        sh.queryEdgeForEdges(
            this->V.row(co_edge.first), this->V.row(co_edge.second),
            eta, sEdgeInds);
        for (const int& seI : sEdgeInds) {
#else
        for (int seI = 0; seI < mesh.SFEdges.size(); ++seI) {
#endif
            // Get the FEM mesh edge as a pair of indiced into mesh.V
            const std::pair<int, int>& mesh_edge = mesh.SFEdges[seI];
            Eigen::Matrix<double, dim, 1> mesh_edge_dir = (mesh.V.row(mesh_edge.second) - mesh.V.row(mesh_edge.first)).normalized();

            // Skip parallel edges
            const double tol = 1e-12;
            if (abs(mesh_edge_dir.dot(co_edge_dir)) > 1 - tol) {
                continue;
            }

            // Check the orientation of the edge.
            MMCVID mmcvid;
            double test_volume;
            compute_c<dim>(
                mesh.V_prev.row(mesh_edge.first),
                mesh.V_prev.row(mesh_edge.second),
                Base::V.row(co_edge.first),
                Base::V.row(co_edge.second),
                test_volume, 1.0);
            if (test_volume >= 0) {
                mmcvid = MMCVID(
                    mesh_edge.first, mesh_edge.second,
                    co_edge.first, co_edge.second);
            }
            else {
                mmcvid = MMCVID(
                    mesh_edge.second, mesh_edge.first,
                    co_edge.first, co_edge.second);
            }

            double toi;
            bool intersects;
            switch (ccdMethod) {
            case ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER:
                intersects = CTCD::edgeEdgeCTCD(
                    mesh.V_prev.row(mesh_edge.first).transpose(),
                    mesh.V_prev.row(mesh_edge.second).transpose(),
                    Base::V.row(co_edge.first).transpose(),
                    Base::V.row(co_edge.second).transpose(),
                    mesh.V.row(mesh_edge.first).transpose() + searchDir.segment<dim>(mesh_edge.first * dim),
                    mesh.V.row(mesh_edge.second).transpose() + searchDir.segment<dim>(mesh_edge.second * dim),
                    Base::V.row(co_edge.first).transpose(),
                    Base::V.row(co_edge.second).transpose(),
                    eta,
                    toi);
                break;
            case ccd::CCDMethod::TIGHT_INCLUSION: {
#ifdef IPC_WITH_TIGHT_INCLUSION
                double output_tolerance;
                intersects = inclusion_ccd::edgeEdgeCCD_double(
                    mesh.V_prev.row(mesh_edge.first).transpose(),
                    mesh.V_prev.row(mesh_edge.second).transpose(),
                    Base::V.row(co_edge.first).transpose(),
                    Base::V.row(co_edge.second).transpose(),
                    mesh.V.row(mesh_edge.first).transpose() + searchDir.segment<dim>(mesh_edge.first * dim),
                    mesh.V.row(mesh_edge.second).transpose() + searchDir.segment<dim>(mesh_edge.second * dim),
                    Base::V.row(co_edge.first).transpose(),
                    Base::V.row(co_edge.second).transpose(),
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
                    mesh.V_prev.row(mesh_edge.first).transpose(),
                    mesh.V_prev.row(mesh_edge.second).transpose(),
                    Base::V.row(co_edge.first).transpose(),
                    Base::V.row(co_edge.second).transpose(),
                    mesh.V.row(mesh_edge.first).transpose() + searchDir.segment<dim>(mesh_edge.first * dim),
                    mesh.V.row(mesh_edge.second).transpose() + searchDir.segment<dim>(mesh_edge.second * dim),
                    Base::V.row(co_edge.first).transpose(),
                    Base::V.row(co_edge.second).transpose(),
                    ccdMethod,
                    toi);
            }

            if (intersects && constraintType == CollisionConstraintType::VERSCHOOR) {
                Eigen::Matrix<double, dim, 1> mesh_v0_toi = (mesh.V.row(mesh_edge.first) - mesh.V_prev.row(mesh_edge.first)) * toi + mesh.V_prev.row(mesh_edge.first);
                Eigen::Matrix<double, dim, 1> mesh_v1_toi = (mesh.V.row(mesh_edge.second) - mesh.V_prev.row(mesh_edge.second)) * toi + mesh.V_prev.row(mesh_edge.second);
                mesh_edge_dir = (mesh_v1_toi - mesh_v0_toi).normalized();
                if (abs(mesh_edge_dir.dot(co_edge_dir)) > 1 - tol) {
                    continue;
                }
            }

            mmcvid_to_toi.insert_or_assign(mmcvid, intersects ? toi : std::numeric_limits<double>::infinity());
            bool isNewConstraint = prevActiveSet.find(mmcvid) == prevActiveSet.end();
            if (intersects && toi < min_tois[seI] && (wasActiveSetCleared || isNewConstraint)) {
                newConstraintsAdded |= isNewConstraint;
                min_tois[seI] = toi;
                mmcvids[seI] = mmcvid;
            }
        }
    }

    addValidConstraints(min_tois, mmcvids, activeSet);

    return newConstraintsAdded;
}

// Determine if the mesh intersects with this collision object.
template <int dim>
bool MeshCO<dim>::isIntersected(
    const Mesh<dim>& mesh,
    const Eigen::MatrixXd& V0,
    const ccd::CCDMethod ccdMethod) const
{
#ifdef USE_EXACT_CCD
    if (ccdMethod == ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER) {
        return false;
    }

    // TODO: spatial hashing
#ifdef USE_SH_INTERSECTED
    Eigen::MatrixXd searchDir = mesh.V - V0;
    searchDir = Map<VectorXd>(searchDir.data(), searchDir.size());
    SpatialHash<dim> sh(mesh, -searchDir, 1, mesh.avgEdgeLen / 3.0);
#endif

    // point-triangle
    // Loop over mesh collision object faces
    for (int MCOFI = 0; MCOFI < this->F.rows(); ++MCOFI) {
        const RowVector3i& MCTriVInd = this->F.row(MCOFI);
        // Loop over FEM mesh surface vertices
#ifdef USE_SH_INTERSECTED
        std::unordered_set<int> pointInds; // NOTE: different constraint order will result in numerically different results
        sh.queryTriangleForPoints(
            this->V.row(MCTriVInd[0]),
            this->V.row(MCTriVInd[1]),
            this->V.row(MCTriVInd[2]),
            /*radius=*/0, pointInds);
        for (const int& svI : pointInds) {
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
#endif
            int vI = mesh.SVI[svI];
            if (ccd::vertexFaceCCD(
                    V0.row(vI).transpose(),
                    this->V.row(MCTriVInd[0]).transpose(),
                    this->V.row(MCTriVInd[1]).transpose(),
                    this->V.row(MCTriVInd[2]).transpose(),
                    mesh.V.row(vI).transpose(),
                    this->V.row(MCTriVInd[0]).transpose(),
                    this->V.row(MCTriVInd[1]).transpose(),
                    this->V.row(MCTriVInd[2]).transpose(),
                    ccdMethod)) {
                return true;
            }
        }
    }

    // triangle-point
    // Loop over mesh collision object vertices
    for (int vI = 0; vI < this->V.rows(); ++vI) {
        // Loop over FEM mesh surface faces
#ifdef USE_SH_INTERSECTED
        std::unordered_set<int> triInds; // NOTE: different constraint order will result in numerically different results
        sh.queryPointForTriangles(Base::V.row(vI), /*radius=*/0, triInds);
        for (const int& sfI : triInds) {
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI) {
#endif
            const RowVector3i& sfVInd = mesh.SF.row(sfI);
            if (ccd::vertexFaceCCD(
                    this->V.row(vI).transpose(),
                    V0.row(sfVInd[0]).transpose(),
                    V0.row(sfVInd[1]).transpose(),
                    V0.row(sfVInd[2]).transpose(),
                    this->V.row(vI).transpose(),
                    mesh.V.row(sfVInd[0]).transpose(),
                    mesh.V.row(sfVInd[1]).transpose(),
                    mesh.V.row(sfVInd[2]).transpose(),
                    ccdMethod)) {
                return true;
            }
        }
    }

    // edge-edge
    // Loop over meshCO edges
    // Get the mesh CO edge as a pair of indiced into this->V
    for (const std::pair<int, int>& co_edge : edges) {
        // Loop over FEM mesh edges
        // Get the FEM mesh edge as a pair of indiced into mesh.V
#ifdef USE_SH_INTERSECTED
        std::vector<int> sEdgeInds;
        sh.queryEdgeForEdges(
            this->V.row(co_edge.first), this->V.row(co_edge.second),
            /*radius=*/0, sEdgeInds);
        for (const int& seI : sEdgeInds) {
#else
        for (int seI = 0; seI < mesh.SFEdges.size(); ++seI) {
#endif
            const std::pair<int, int>& mesh_edge = mesh.SFEdges[seI];
            if (ccd::edgeEdgeCCD(
                    V0.row(mesh_edge.first).transpose(),
                    V0.row(mesh_edge.second).transpose(),
                    this->V.row(co_edge.first).transpose(),
                    this->V.row(co_edge.second).transpose(),
                    mesh.V.row(mesh_edge.first).transpose(),
                    mesh.V.row(mesh_edge.second).transpose(),
                    this->V.row(co_edge.first).transpose(),
                    this->V.row(co_edge.second).transpose(),
                    ccdMethod)) {
                return true;
            }
        }
    }
#endif
    return false;
}

template <int dim>
void MeshCO<dim>::filterSearchDir_QP(const Mesh<dim>& mesh,
    Eigen::VectorXd& searchDir,
    std::vector<MMCVID>& activeSet_next)
{
    // TODO: spatial hashing
    // TODO: parallelize

    // point-triangle
    double stepsize = 1.0;
    for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
        int vI = mesh.SVI[svI];
        for (int MCOFI = 0; MCOFI < Base::F.rows(); ++MCOFI) {
            const RowVector3i& MCTriVInd = Base::F.row(MCOFI);
            double c;
            compute_c<dim>(mesh.V.row(vI),
                Base::V.row(MCTriVInd[0]),
                Base::V.row(MCTriVInd[2]),
                Base::V.row(MCTriVInd[1]),
                c, 1.0);
            if (c <= 0.0) {
                continue;
            }

            double largestAlpha = 1.0;
#if (DIM == 3)
            computeInjectiveStepSize_3d(mesh.V.row(vI),
                Base::V.row(MCTriVInd[0]),
                Base::V.row(MCTriVInd[2]),
                Base::V.row(MCTriVInd[1]),
                searchDir.segment<dim>(vI * dim),
                Eigen::Vector3d::Zero(),
                Eigen::Vector3d::Zero(),
                Eigen::Vector3d::Zero(),
                1.0e-8, largestAlpha);
#else

#endif

            if (largestAlpha < 1.0) {
                // check intersection
                if (pointInTriangleSweep<dim>(
                        mesh.V.row(vI) + largestAlpha * searchDir.segment<dim>(vI * dim).transpose(),
                        Base::V.row(MCTriVInd[0]),
                        Base::V.row(MCTriVInd[2]),
                        Base::V.row(MCTriVInd[1]))) {
                    activeSet_next.emplace_back(-vI - 1, MCTriVInd[0], MCTriVInd[2], MCTriVInd[1]);

                    searchDir.segment<dim>(vI * dim) *= largestAlpha; // NOTE: will possibly result in non-descent direction
                }
                else {
                    // std::cout << "non-intersecting constraint 'violation' not counted for largest step size" << std::endl;
                }
            }
        }
    }

    // TODO: triangle-point

    // TODO: edge-edge
}

template <int dim>
void MeshCO<dim>::move(const Eigen::Matrix<double, dim, 1>& deltaX,
    const Mesh<dim>& mesh, const SpatialHash<dim>& sh,
    double slackness, double& stepSizeLeft)
{
    // movement of Base::origin specified in deltaX
    // movement of each vertex specified in V_target,
    // which is consistent with the movement of Base::origin
    assert(Base::V_target.rows() == Base::V.rows());
    Vdist2 = (Base::V_target - Base::V).rowwise().squaredNorm();

    double stepSize = 1.0;
    // point-triangle
    largestAlphaPT.resize(Base::F.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::F.rows(), 1, [&](int sfI)
#else
    for (int sfI = 0; sfI < Base::F.rows(); ++sfI)
#endif
        {
            largestAlphaPT[sfI] = 1.0;
            const RowVector3i& sfVInd = Base::F.row(sfI);
#ifdef USE_SH_CCS
            std::unordered_set<int> pointInds; // NOTE: different constraint order will result in numerically different results
            sh.queryTriangleForPoints(Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]),
                Base::V_target.row(sfVInd[0]) - Base::V.row(sfVInd[0]),
                Base::V_target.row(sfVInd[1]) - Base::V.row(sfVInd[1]),
                Base::V_target.row(sfVInd[2]) - Base::V.row(sfVInd[2]), pointInds);
            for (const auto& svI : pointInds)
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
            {
                int vI = mesh.SVI[svI];

                double d_sqrt;
                computePointTriD(mesh.V.row(vI), Base::V.row(sfVInd[0]), Base::V.row(sfVInd[1]), Base::V.row(sfVInd[2]), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlpha = 1.0;
                if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                        Base::V.row(sfVInd[0]).transpose(),
                        Base::V.row(sfVInd[1]).transpose(),
                        Base::V.row(sfVInd[2]).transpose(),
                        mesh.V.row(vI).transpose(),
                        Base::V_target.row(sfVInd[0]).transpose(),
                        Base::V_target.row(sfVInd[1]).transpose(),
                        Base::V_target.row(sfVInd[2]).transpose(),
                        0.5 * d_sqrt,
                        largestAlpha)) {
                    if (largestAlpha < 1.0e-6) {
                        if (CTCD::vertexFaceCTCD(mesh.V.row(vI).transpose(),
                                Base::V.row(sfVInd[0]).transpose(),
                                Base::V.row(sfVInd[1]).transpose(),
                                Base::V.row(sfVInd[2]).transpose(),
                                mesh.V.row(vI).transpose(),
                                Base::V_target.row(sfVInd[0]).transpose(),
                                Base::V_target.row(sfVInd[1]).transpose(),
                                Base::V_target.row(sfVInd[2]).transpose(),
                                0.0, largestAlpha)) {
                            largestAlpha *= slackness;
                        }
                        else {
                            largestAlpha = 1.0;
                        }
                    }
                }
                largestAlphaPT[sfI] = std::min(largestAlpha, largestAlphaPT[sfI]);
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphaPT.minCoeff());

    // triangle-point
    largestAlphaTP.resize(Base::V.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)Base::V.rows(), 1, [&](int vI)
#else
    for (int vI = 0; vI < Base::V.rows(); ++vI)
#endif
        {
            largestAlphaTP[vI] = 1.0;
#ifdef USE_SH_CCS
            std::unordered_set<int> triInds; // NOTE: different constraint order will result in numerically different results
            sh.queryPointForTriangles(Base::V.row(vI), Base::V_target.row(vI) - Base::V.row(vI), triInds);
            for (const auto& sfI : triInds)
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
            {
                const RowVector3i& sfVInd = mesh.SF.row(sfI);

                double d_sqrt;
                computePointTriD(Base::V.row(vI), mesh.V.row(sfVInd[0]), mesh.V.row(sfVInd[1]), mesh.V.row(sfVInd[2]), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlpha = 1.0;
                if (CTCD::vertexFaceCTCD(Base::V.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose(),
                        mesh.V.row(sfVInd[1]).transpose(),
                        mesh.V.row(sfVInd[2]).transpose(),
                        Base::V_target.row(vI).transpose(),
                        mesh.V.row(sfVInd[0]).transpose(),
                        mesh.V.row(sfVInd[1]).transpose(),
                        mesh.V.row(sfVInd[2]).transpose(),
                        0.5 * d_sqrt,
                        largestAlpha)) {
                    if (largestAlpha < 1.0e-6) {
                        if (CTCD::vertexFaceCTCD(Base::V.row(vI).transpose(),
                                mesh.V.row(sfVInd[0]).transpose(),
                                mesh.V.row(sfVInd[1]).transpose(),
                                mesh.V.row(sfVInd[2]).transpose(),
                                Base::V_target.row(vI).transpose(),
                                mesh.V.row(sfVInd[0]).transpose(),
                                mesh.V.row(sfVInd[1]).transpose(),
                                mesh.V.row(sfVInd[2]).transpose(),
                                0.0, largestAlpha)) {
                            largestAlpha *= slackness;
                        }
                        else {
                            largestAlpha = 1.0;
                        }
                    }
                }
                largestAlphaTP[vI] = std::min(largestAlpha, largestAlphaTP[vI]);
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphaTP.minCoeff());

    // edge-edge
    largestAlphaEE.resize(edges.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)edges.size(), 1, [&](int eJ)
#else
    for (int eJ = 0; eJ < edges.size(); ++eJ)
#endif
        {
            largestAlphaEE[eJ] = 1.0;
            const auto& meshEJ = edges[eJ];
#ifdef USE_SH_CCS
            std::vector<int> edgeInds; // NOTE: different constraint order will result in numerically different results
            sh.queryEdgeForEdges(Base::V.row(meshEJ.first), Base::V.row(meshEJ.second),
                Base::V_target.row(meshEJ.first) - Base::V.row(meshEJ.first),
                Base::V_target.row(meshEJ.second) - Base::V.row(meshEJ.second), edgeInds);
            for (const auto& eI : edgeInds) {
                const auto& meshEI = mesh.SFEdges[eI];
#else
        for (const auto& meshEI : mesh.SFEdges) {
#endif
                double d_sqrt;
                computeEdgeEdgeD(mesh.V.row(meshEI.first), mesh.V.row(meshEI.second),
                    Base::V.row(meshEJ.first), Base::V.row(meshEJ.second), d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlphas = 1.0;
                if (CTCD::edgeEdgeCTCD(mesh.V.row(meshEI.first).transpose(),
                        mesh.V.row(meshEI.second).transpose(),
                        Base::V.row(meshEJ.first).transpose(),
                        Base::V.row(meshEJ.second).transpose(),
                        mesh.V.row(meshEI.first).transpose(),
                        mesh.V.row(meshEI.second).transpose(),
                        Base::V_target.row(meshEJ.first).transpose(),
                        Base::V_target.row(meshEJ.second).transpose(),
                        0.5 * d_sqrt,
                        largestAlphas)) {
                    if (largestAlphas < 1.0e-6) {
                        if (CTCD::edgeEdgeCTCD(mesh.V.row(meshEI.first).transpose(),
                                mesh.V.row(meshEI.second).transpose(),
                                Base::V.row(meshEJ.first).transpose(),
                                Base::V.row(meshEJ.second).transpose(),
                                mesh.V.row(meshEI.first).transpose(),
                                mesh.V.row(meshEI.second).transpose(),
                                Base::V_target.row(meshEJ.first).transpose(),
                                Base::V_target.row(meshEJ.second).transpose(),
                                0.0, largestAlphas)) {
                            if (largestAlphas == 0.0) {
                                // numerically parallel edge-edge causes CCD code to fail
                                largestAlphas = 1.0;
                            }
                            largestAlphas *= slackness;
                        }
                        else {
                            largestAlphas = 1.0;
                        }
                    }
                }
                largestAlphaEE[eJ] = std::min(largestAlphas, largestAlphaEE[eJ]);
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, largestAlphaEE.minCoeff());

    Base::origin += deltaX * stepSize;
    Base::V += (Base::V_target - Base::V) * stepSize;
    stepSizeLeft = 1.0 - stepSize;
}

template class MeshCO<DIM>;

} // namespace IPC
