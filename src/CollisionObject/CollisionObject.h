//
//  CollisionObject.h
//  IPC
//
//  Created by Minchen Li on 10/22/18.
//

#ifndef CollisionObject_h
#define CollisionObject_h

#include "Energy.hpp"
#include "Mesh.hpp"
#include "LinSysSolver.hpp"
#include "MeshCollisionUtils.hpp"
#include "CollisionConstraints.hpp"
#include "SpatialHash.hpp"

#include <ccd.hpp>

#include <igl/writeOBJ.h>

#include <Eigen/Eigen>

#include <vector>

namespace IPC {

template <int dim>
class CollisionObject {
public:
    Eigen::Matrix<double, dim, 1> origin;
    Eigen::Matrix<double, dim, 1> velocitydt;
    double friction;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd V_target; // for moving meshCO

public:
    virtual ~CollisionObject(void){};

public:
    virtual void updateConstraints_QP(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        std::vector<Eigen::Triplet<double>>& A_triplet,
        Eigen::VectorXd& l) const = 0;

    virtual void evaluateConstraint(const Mesh<dim>& mesh,
        int vI, double& val, double coef = 1.0) const = 0;
    virtual void evaluateConstraints(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        Eigen::VectorXd& val, double coef = 1.0) const
    {
        int constraintStartInd = val.size();
        val.conservativeResize(constraintStartInd + activeSet.size());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)activeSet.size(), 1, [&](int cI)
#else
        for (int cI = 0; cI < activeSet.size(); ++cI)
#endif
            {
                int vI = activeSet[cI];
                evaluateConstraint(mesh, vI, val[constraintStartInd + cI], coef);
            }
#ifdef USE_TBB
        );
#endif
    }

    virtual void evaluateConstraints_all(const Mesh<dim>& mesh,
        Eigen::VectorXd& val, double coef = 1.0) const
    {
        int constraintStartInd = val.size();
        val.conservativeResize(constraintStartInd + mesh.V.rows());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.V.rows(), 1, [&](int vI)
#else
        for (int vI = 0; vI < mesh.V.rows(); vI++)
#endif
            {
                evaluateConstraint(mesh, vI, val[constraintStartInd + vI], coef);
            }
#ifdef USE_TBB
        );
#endif
    }

    virtual void evaluateConstraintQP(const Mesh<dim>& mesh,
        int vI, double& val, double coef = 1.0) const = 0;
    virtual void evaluateConstraintsQP(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        Eigen::VectorXd& val, double coef = 1.0) const
    {
        int constraintStartInd = val.size();
        val.conservativeResize(constraintStartInd + activeSet.size());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)activeSet.size(), 1, [&](int cI)
#else
        for (int cI = 0; cI < activeSet.size(); ++cI)
#endif
            {
                int vI = activeSet[cI];
                evaluateConstraintQP(mesh, vI, val[constraintStartInd + cI], coef);
            }
#ifdef USE_TBB
        );
#endif
    }
    virtual void evaluateConstraintsQP_all(const Mesh<dim>& mesh,
        Eigen::VectorXd& val, double coef = 1.0) const
    {
        int constraintStartInd = val.size();
        val.conservativeResize(constraintStartInd + mesh.V.rows());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.V.rows(), 1, [&](int vI)
#else
        for (int vI = 0; vI < mesh.V.rows(); vI++)
#endif
            {
                evaluateConstraintQP(mesh, vI, val[constraintStartInd + vI], coef);
            }
#ifdef USE_TBB
        );
#endif
    }

    virtual void leftMultiplyConstraintJacobianT(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output,
        double coef = 1.0) const = 0;
    virtual void leftMultiplyConstraintJacobianTQP(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output,
        double coef = 1.0) const = 0;
    virtual void augmentIPHessian(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* mtr_incremental,
        double dHat, double coef = 1.0, bool projectDBC = true) const = 0;

    /**
     * @brief Update the active set of constraints for the QP/SQP solver.
     *
     * @param[in]  mesh            FEM mesh colliding with this mesh collision object.
     * @param[in]  searchDir       Linear displacement of the mesh vertices.
     * @param[in]  constraintType  Method of updating the active set.
     * @param[out] activeSet       Constraint active set to update.
     * @param[in]  eta             Collision offset distance for early activation.
     */
    virtual bool updateActiveSet_QP(
        const Mesh<dim>& mesh,
        const Eigen::VectorXd& searchDir,
        const CollisionConstraintType constraintType,
        std::vector<MMCVID>& activeSet,
        const ccd::CCDMethod ccdMethod,
        const double eta = 0,
        const double ccd_tol = 1e-6) { return false; };

    virtual void filterSearchDir_QP(const Mesh<dim>& mesh,
        Eigen::VectorXd& searchDir,
        std::vector<int>& activeSet_next)
        = 0;
    virtual void largestFeasibleStepSize(const Mesh<dim>& mesh,
        const Eigen::VectorXd& searchDir,
        double slackness,
        std::vector<int>& activeSet_next,
        double& stepSize)
        = 0;

    // for mesh collision objects
    virtual void evaluateConstraint(const Mesh<dim>& mesh,
        const MMCVID& MMCVIDI, double& val, double coef = 1.0) const {}

    virtual void evaluateConstraints(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        Eigen::VectorXd& val, double coef = 1.0) const
    {
        int constraintStartInd = val.size();
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

    virtual void leftMultiplyConstraintJacobianT(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output_incremental,
        double coef = 1.0) const {}

    /**
     * @brief Compute the collision constraints used for SQP method.
     *
     * @param[in]  mesh            FEM mesh colliding with this object.
     * @param[in]  MMCVIDI         Constraint pair.
     * @param[in]  constraintType  Type of constraints to use.
     * @param[out] val             Computed constraint value.
     * @param[in]  coef            Coefficent to multiply the value.
     */
    virtual void evaluateConstraintQP(
        const Mesh<dim>& mesh, const MMCVID& MMCVIDI,
        const CollisionConstraintType constraintType, double& val,
        double coef) const {}

    /**
     * @brief Compute the collision constraints used for SQP method.
     *
     * @param[in]  mesh            FEM mesh colliding with this object.
     * @param[in]  activeSet       Active set of constraint pair.
     * @param[in]  constraintType  Type of constraints to use.
     * @param[out] val             Computed constraint values.
     * @param[in]  coef            Coefficent to multiply the values.
     */
    virtual void evaluateConstraintsQP(
        const Mesh<dim>& mesh, const std::vector<MMCVID>& activeSet,
        const CollisionConstraintType constraintType, Eigen::VectorXd& val,
        double coef = 1.0) const
    {
        int constraintStartInd = val.size();
        val.conservativeResize(constraintStartInd + activeSet.size());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)activeSet.size(), 1, [&](int cI)
#else
        for (int cI = 0; cI < activeSet.size(); ++cI)
#endif
            {
                evaluateConstraintQP(
                    mesh, activeSet[cI], constraintType,
                    val[constraintStartInd + cI], coef);
            }
#ifdef USE_TBB
        );
#endif
    }

    /**
     * @brief Evaluate the Jacobian of the constraints for the SQP method.
     *
     * @param[in]  mesh                FEM mesh colliding with this object.
     * @param[in]  activeSet           Set of constraint identifier.
     * @param[in]  input               Input value.
     * @param[in]  constraintType      Type of constraints to use.
     * @param[out] output_incremental  Computed output of Jacobian times input.
     * @param[in]  coef                Ceofficient of the constraints.
     */
    virtual void leftMultiplyConstraintJacobianTQP(
        const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        const Eigen::VectorXd& input,
        const CollisionConstraintType constraintType,
        Eigen::VectorXd& output_incremental,
        double coef = 1.0) const {}

    virtual void augmentIPHessian(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* mtr_incremental,
        double dHat, double coef = 1.0, bool projectDBC = true) const {}

    virtual void largestFeasibleStepSize(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        double slackness,
        const std::vector<std::pair<int, int>>& constraintSet,
        double& stepSize) {}
    virtual void largestFeasibleStepSize_CCD(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh, const Eigen::VectorXd& searchDir,
        double slackness, double& stepSize) {}

#ifdef IPC_WITH_TIGHT_INCLUSION
    virtual void largestFeasibleStepSize_TightInclusion(
        const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        double tolerance,
        const std::vector<std::pair<int, int>>& constraintSet,
        double& stepSize)
    {
    }
    virtual void largestFeasibleStepSize_CCD_TightInclusion(
        const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        double tolerance,
        double& stepSize)
    {
    }
#endif

    virtual void largestFeasibleStepSize_exact(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        ccd::CCDMethod method,
        const std::vector<std::pair<int, int>>& constraintSet,
        double& stepSize)
    {
    }
    virtual void largestFeasibleStepSize_CCD_exact(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh, const Eigen::VectorXd& searchDir,
        ccd::CCDMethod method, double& stepSize) {}

    // for meshCO when SQP is used
    virtual void updateConstraints_QP(
        const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        const CollisionConstraintType constraintType,
        std::vector<Eigen::Triplet<double>>& A_triplet,
        Eigen::VectorXd& l) const {}

    virtual void filterSearchDir_QP(const Mesh<dim>& mesh,
        Eigen::VectorXd& searchDir,
        std::vector<MMCVID>& activeSet_next) {}

    // for distance-based IP
    virtual void computeConstraintSet(const Mesh<dim>& mesh,
        double dHat, std::vector<int>& constraintSet) const
    {
        std::vector<int> isActive(mesh.SVI.size(), 0);
#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.SVI.size(), 1, [&](int svI)
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
            {
                int vI = mesh.SVI[svI];
                if (!mesh.isDBCVertex(vI) && mesh.vICoDim(vI) == 3) {
                    double d = 0;
                    evaluateConstraint(mesh, vI, d);
                    if (d < dHat) {
                        isActive[svI] = 1;
                    }
                }
            }
#ifdef USE_TBB
        );
#endif

        constraintSet.resize(0);
        for (int svI = 0; svI < isActive.size(); ++svI) {
            if (isActive[svI]) {
                constraintSet.emplace_back(mesh.SVI[svI]);
            }
        }
    }
    virtual void computeConstraintSet(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh, double dHat,
        std::vector<MMCVID>& constraintSet,
        std::vector<MMCVID>& paraEEMMCVIDSet,
        std::vector<std::pair<int, int>>& paraEEeIeJSet,
        bool getPTEE, std::vector<std::pair<int, int>>& cs_PTEE) const {}

    virtual void augmentParaEEEnergy(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& paraEEMMCVIDSet,
        const std::vector<std::pair<int, int>>& paraEEeIeJSet,
        Eigen::VectorXd& d_inc, Eigen::VectorXd& energyVec_inc,
        double dHat, double coef) const {}

    virtual void augmentParaEEGradient(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& paraEEMMCVIDSet,
        const std::vector<std::pair<int, int>>& paraEEeIeJSet,
        Eigen::VectorXd& grad_inc, double dHat, double coef) const {}

    virtual void augmentParaEEHessian(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& paraEEMMCVIDSet,
        const std::vector<std::pair<int, int>>& paraEEeIeJSet,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* H_inc,
        double dHat, double coef, bool projectDBC) const {}

    /**
     * @brief Determine if the mesh intersects with this collision object.
     *
     * @param[in]  mesh    FEM mesh colliding with this collision object.
     * @param[in]  V0      FEM mesh vertices at previous iteration.
     * @param[in]  method  NOT USED
     *
     * @returns  True if the mesh intersect with this collision object.
     */
    virtual bool isIntersected(
        const Mesh<dim>& mesh,
        const Eigen::MatrixXd& V0,
        const ccd::CCDMethod method = ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER) const
    {
        Eigen::VectorXd constraint_vals;
        this->evaluateConstraints_all(mesh, constraint_vals, 1.0);
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.vICoDim(vI) == dim && !mesh.isDBCVertex(vI)) {
                if (constraint_vals[vI] <= 0.0) {
                    return true;
                }
            }
        }
        return false;
    }

    virtual void computeFrictionEnergy(const Eigen::MatrixXd& V,
        const Eigen::MatrixXd& Vt, const std::vector<int>& activeSet,
        const Eigen::VectorXd& multipliers,
        double& Ef, double eps2, double coef) const
    {
        throw "computeFrictionEnergy not implemented!";
    }
    virtual void augmentFrictionGradient(const Eigen::MatrixXd& V,
        const Eigen::MatrixXd& Vt, const std::vector<int>& activeSet,
        const Eigen::VectorXd& multipliers,
        Eigen::VectorXd& grad_inc, double eps2, double coef) const
    {
        throw "augmentFrictionGradient not implemented!";
    }
    virtual void augmentFrictionHessian(const Mesh<dim>& mesh,
        const Eigen::MatrixXd& Vt, const std::vector<int>& activeSet,
        const Eigen::VectorXd& multipliers,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* H_inc,
        double eps2, double coef, bool projectDBC = true) const
    {
        throw "augmentFrictionHessian not implemented!";
    }

    virtual void computeRelDX(const Eigen::RowVector3d& relDX3D, Eigen::Vector2d& relDX) const
    {
        throw "computeRelDX not implemented!";
    }

    virtual bool checkEdgeTriIntersection(
        const Mesh<dim>& mesh, const SpatialHash<dim>& sh)
    {
        throw "checkEdgeTriIntersection not implemented!";
    }
    virtual bool checkEdgeTriIntersectionIfAny(
        const Mesh<dim>& mesh, const SpatialHash<dim>& sh)
    {
        throw "checkEdgeTriIntersectionIfAny not implemented!";
    }

    virtual void move(const Eigen::Matrix<double, dim, 1>& deltaX,
        const Mesh<dim>& mesh, const SpatialHash<dim>& sh,
        double slackness, double& stepSizeLeft)
    {
        throw "move not implemented!";
    }

    virtual void initRenderingData(double extensionScale = 1.0) = 0;

    virtual void draw(Eigen::MatrixXd& p_V,
        Eigen::MatrixXi& p_F,
        Eigen::MatrixXd& color,
        double extensionScale = 1.0) const
    {
        int oldVSize = p_V.rows();
        p_V.conservativeResize(oldVSize + V.rows(), 3);
        p_V.bottomRows(V.rows()) = V;

        int oldFSize = p_F.rows();
        p_F.conservativeResize(oldFSize + F.rows(), 3);
        p_F.bottomRows(F.rows()) = F;
        p_F.bottomRows(F.rows()).array() += oldVSize;

        color.conservativeResize(color.rows() + F.rows(), 3);
        color.bottomRows(F.rows()).setConstant(0.9);
    }

    virtual void saveMesh(const std::string& filePath, const Eigen::Vector3d& shift = Eigen::Vector3d::Zero()) const
    {
        igl::writeOBJ(filePath, V.rowwise() + shift.transpose(), F);
    }
};

} // namespace IPC

#endif /* CollisionObject_h */
