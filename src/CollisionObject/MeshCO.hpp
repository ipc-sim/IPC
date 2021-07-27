//
//  MeshCO.hpp
//  IPC
//
//  Created by Minchen Li on 4/25/19.
//

#ifndef MeshCO_hpp
#define MeshCO_hpp

#include "CollisionObject.h"

#include "Energy.hpp"

namespace IPC {

template <int dim>
class MeshCO : public CollisionObject<dim> {
    typedef CollisionObject<dim> Base;

public:
    std::string objFilePath;
    std::vector<std::pair<int, int>> edges;
    double scale;

    Eigen::VectorXd largestAlphaPT, largestAlphaEE, largestAlphaTP, Vdist2;
    // information on how large the dHat needs to be

public:
    MeshCO(const char* objFilePath,
        const Eigen::Matrix<double, dim, 1>& p_position,
        const Eigen::Matrix3d& p_rotation,
        double scale, double p_friction);

protected:
    std::unordered_map<MMCVID, double, MMCVIDHash> mmcvid_to_toi;

public:
    virtual void updateConstraints_QP(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        std::vector<Eigen::Triplet<double>>& A_triplet,
        Eigen::VectorXd& l) const override {}

    virtual void evaluateConstraint(const Mesh<dim>& mesh,
        int vI, double& val, double coef = 1.0) const override {}
    virtual void evaluateConstraintQP(const Mesh<dim>& mesh,
        int vI, double& val, double coef = 1.0) const override {}

    virtual void leftMultiplyConstraintJacobianT(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output_incremental,
        double coef = 1.0) const override {}
    virtual void leftMultiplyConstraintJacobianTQP(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output_incremental,
        double coef = 1.0) const override {}
    virtual void augmentIPHessian(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* mtr_incremental,
        double dHat, double coef = 1.0, bool projectDBC = true) const override {}

    virtual void filterSearchDir_QP(const Mesh<dim>& mesh,
        Eigen::VectorXd& searchDir,
        std::vector<int>& activeSet_next) override {}
    virtual void largestFeasibleStepSize(const Mesh<dim>& mesh,
        const Eigen::VectorXd& searchDir,
        double slackness,
        std::vector<int>& activeSet_next,
        double& stepSize) override {}

    void evaluateConstraint(const Mesh<dim>& mesh,
        const MMCVID& MMCVIDI, double& val, double coef = 1.0) const override;

    void leftMultiplyConstraintJacobianT(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output_incremental,
        double coef = 1.0) const override;

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
        double coef) const override;

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
        double coef = 1.0) const override;

    void augmentIPHessian(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* mtr_incremental,
        double dHat, double coef = 1.0, bool projectDBC = true) const override;

    void largestFeasibleStepSize(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        double slackness,
        const std::vector<std::pair<int, int>>& constraintSet,
        double& stepSize) override;
    void largestFeasibleStepSize_CCD(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh, const Eigen::VectorXd& searchDir,
        double slackness, double& stepSize) override;

#ifdef IPC_WITH_TIGHT_INCLUSION
    void largestFeasibleStepSize_TightInclusion(
        const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        double tolerance,
        const std::vector<std::pair<int, int>>& constraintSet,
        double& stepSize) override;
    void largestFeasibleStepSize_CCD_TightInclusion(
        const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        double tolerance,
        double& stepSize) override;
#endif

    void largestFeasibleStepSize_exact(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        const ccd::CCDMethod method,
        const std::vector<std::pair<int, int>>& constraintSet,
        double& stepSize) override;
    void largestFeasibleStepSize_CCD_exact(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh, const Eigen::VectorXd& searchDir,
        const ccd::CCDMethod method, double& stepSize) override;

    void computeConstraintSet(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh, double dHat,
        std::vector<MMCVID>& constraintSet,
        std::vector<MMCVID>& paraEEMMCVIDSet,
        std::vector<std::pair<int, int>>& paraEEeIeJSet,
        bool getPTEE, std::vector<std::pair<int, int>>& cs_PTEE) const override;

    void augmentParaEEEnergy(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& paraEEMMCVIDSet,
        const std::vector<std::pair<int, int>>& paraEEeIeJSet,
        Eigen::VectorXd& d_inc, Eigen::VectorXd& energyVec_inc,
        double dHat, double coef) const override;

    void augmentParaEEGradient(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& paraEEMMCVIDSet,
        const std::vector<std::pair<int, int>>& paraEEeIeJSet,
        Eigen::VectorXd& grad_inc, double dHat, double coef) const override;

    void augmentParaEEHessian(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& paraEEMMCVIDSet,
        const std::vector<std::pair<int, int>>& paraEEeIeJSet,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* H_inc,
        double dHat, double coef, bool projectDBC) const override;

    bool checkEdgeTriIntersection(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh) override;
    bool checkEdgeTriIntersectionIfAny(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh) override;

    // when SQP is used
    virtual void updateConstraints_QP(
        const Mesh<dim>& mesh, const std::vector<MMCVID>& activeSet,
        const CollisionConstraintType constraintType,
        std::vector<Eigen::Triplet<double>>& A_triplet, Eigen::VectorXd& l) const override;

    /**
     * @brief Update the active set of constraints
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
        const double ccd_tol = 1e-6) override;

    virtual void filterSearchDir_QP(
        const Mesh<dim>& mesh, Eigen::VectorXd& searchDir,
        std::vector<MMCVID>& activeSet_next) override;

    /**
     * @brief Determine if the mesh intersects with this collision object.
     *
     * Check if the mesh intersects with the collision object along the linear
     * trajectories from V0 to mesh.V
     *
     * @param[in]  mesh    FEM mesh colliding with this collision object.
     * @param[in]  V0      FEM mesh vertices at previous iteration.
     * @param[in]  method  Method of exact CCD to use.
     *
     * @returns  True if the mesh intersect with this collision object.
     */
    virtual bool isIntersected(
        const Mesh<dim>& mesh,
        const Eigen::MatrixXd& V0,
        ccd::CCDMethod method) const override;

    virtual void move(const Eigen::Matrix<double, dim, 1>& deltaX,
        const Mesh<dim>& mesh, const SpatialHash<dim>& sh,
        double slackness, double& stepSizeLeft) override;

    virtual void initRenderingData(double extensionScale = 1.0) override{};
};

} // namespace IPC

#endif /* MeshCO_hpp */
