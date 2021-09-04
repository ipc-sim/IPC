//
//  SelfCollisionHandler.hpp
//  IPC
//
//  Created by Minchen Li on 6/03/19.
//

#ifndef SelfCollisionHandler_hpp
#define SelfCollisionHandler_hpp

#include "MeshCollisionUtils.hpp"
#include "CollisionConstraints.hpp"
#include "SpatialHash.hpp"
#include "LinSysSolver.hpp"
#include "Mesh.hpp"
#include <ccd.hpp>

namespace IPC {

template <int dim>
class SelfCollisionHandler {
public:
    static void evaluateConstraint(const Mesh<dim>& mesh,
        const MMCVID& MMCVIDI, double& val, double coef = 1.0);
    static void evaluateConstraints(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        Eigen::VectorXd& val, double coef = 1.0);

    static void leftMultiplyConstraintJacobianT(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output_incremental,
        double coef = 1.0);

    /**
     * @brief Evaluate the constraint for the SQP method.
     *
     * @param[in]  mesh         Mesh to compute the self collision constraints.
     * @param[in]  MMCVIDI         Constraint identifier.
     * @param[in]  constraintType  Type of constraints to use.
     * @param[in]  toi             Constraint's toi used by Verschoor.
     * @param[out] val             Computed value of the constraint.
     * @param[in]  coef            Ceofficient of the constraint.
     */
    static void evaluateConstraintQP(
        const Mesh<dim>& mesh,
        const MMCVID& MMCVIDI,
        const CollisionConstraintType constraintType,
        double toi,
        double& val,
        double coef = 1.0);

    /**
     * @brief Evaluate the constraints for the SQP method.
     *
     * @param[in]  mesh         Mesh to compute the self collision constraints.
     * @param[in]  activeSet       Set of constraint identifier.
     * @param[in]  constraintType  Type of constraints to use.
     * @param[in]  mmcvid_to_toi   Map of constraints' toi used by Verschoor.
     * @param[out] val             Computed value of the constraints.
     * @param[in]  coef            Ceofficient of the constraints.
     */
    static void evaluateConstraintsQP(
        const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        const CollisionConstraintType constraintType,
        const std::unordered_map<MMCVID, double, MMCVIDHash>& mmcvid_to_toi,
        Eigen::VectorXd& val,
        double coef = 1.0);

    /**
     * @brief Evaluate the Jacobian of the constraints for the SQP method.
     *
     * @param[in]  mesh  Mesh to compute the self collision constraints.
     * @param[in]  activeSet           Set of constraint identifier.
     * @param[in]  input               Input value.
     * @param[in]  constraintType      Type of constraints to use.
     * @param[in]  mmcvid_to_toi       Map of constraints' toi used by Verschoor.
     * @param[out] output_incremental  Computed output of Jacobian times input.
     * @param[in]  coef                Ceofficient of the constraints.
     */
    static void leftMultiplyConstraintJacobianTQP(
        const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        const Eigen::VectorXd& input,
        const CollisionConstraintType constraintType,
        const std::unordered_map<MMCVID, double, MMCVIDHash>& mmcvid_to_toi,
        Eigen::VectorXd& output_incremental,
        double coef = 1.0);

    static void augmentConnectivity(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        std::vector<std::set<int>>& vNeighbor);
    static void augmentConnectivity(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& paraEEMMCVIDSet,
        const std::vector<std::pair<int, int>>& paraEEeIeJSet,
        std::vector<std::set<int>>& vNeighbor);

    static void augmentIPHessian(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* mtr_incremental,
        double dHat, double coef = 1.0, bool projectDBC = true);

    // CFL based
    static void largestFeasibleStepSize(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        double slackness,
        const std::vector<std::pair<int, int>>& constraintSet,
        std::vector<std::pair<int, int>>& candidates,
        double& stepSize);
    // CCD based
    static void largestFeasibleStepSize_CCD(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh, const Eigen::VectorXd& searchDir,
        double slackness, std::vector<std::pair<int, int>>& candidates,
        double& stepSize);

#ifdef IPC_WITH_TIGHT_INCLUSION
    static void largestFeasibleStepSize_TightInclusion(
        const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        double tolerance,
        const std::vector<std::pair<int, int>>& constraintSet,
        std::vector<std::pair<int, int>>& candidates,
        double& stepSize);
    static void largestFeasibleStepSize_CCD_TightInclusion(
        const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        double tolerance,
        std::vector<std::pair<int, int>>& candidates,
        double& stepSize);
#endif

    static void largestFeasibleStepSize_exact(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        const Eigen::VectorXd& searchDir,
        const ccd::CCDMethod method,
        const std::vector<std::pair<int, int>>& constraintSet,
        double& stepSize);
    static void largestFeasibleStepSize_CCD_exact(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh, const Eigen::VectorXd& searchDir,
        const ccd::CCDMethod method, double& stepSize);

    static void updateConstraints_QP(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& activeSet,
        const CollisionConstraintType constraintType,
        const std::unordered_map<MMCVID, double, MMCVIDHash>& mmcvid_to_toi,
        std::vector<Eigen::Triplet<double>>& A_triplet, Eigen::VectorXd& l);

    static void filterSearchDir_QP(const Mesh<dim>& mesh,
        Eigen::VectorXd& searchDir,
        std::vector<MMCVID>& activeSet_next);

    /**
     * @brief Update the active set of constraints
     *
     * @param[in]  mesh            FEM mesh colliding with this mesh collision object.
     * @param[in]  searchDir       Linear displacement of the mesh vertices.
     * @param[in]  constraintType  Method of updating the active set.
     * @param[out] activeSet       Constraint active set to update.
     * @param[out] mmcvid_to_toi   Map of constrain toi's used by Verschoor.
     * @param[in]  eta             Collision offset distance for early activation.
     */
    static bool updateActiveSet_QP(
        const Mesh<dim>& mesh, const Eigen::VectorXd& searchDir,
        const CollisionConstraintType constraintType,
        std::vector<MMCVID>& activeSet,
        std::unordered_map<MMCVID, double, MMCVIDHash>& mmcvid_to_toi,
        const ccd::CCDMethod ccdMethod,
        const double eta = 0,
        const double ccd_tol = 1e-6);

    static void computeConstraintSet(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        double dHat, std::vector<MMCVID>& constraintSet,
        std::vector<MMCVID>& paraEEMMCVIDSet,
        std::vector<std::pair<int, int>>& paraEEeIeJSet,
        bool getPTEE, std::vector<std::pair<int, int>>& cs_PTEE);

    static void computeDistCoordAndTanBasis(
        const Mesh<dim>& mesh,
        const std::vector<MMCVID>& constraintSet,
        std::vector<Eigen::Vector2d>& MMDistCoord,
        std::vector<Eigen::Matrix<double, 3, 2>>& MMTanBasis);

    static void computeFrictionEnergy(const Eigen::MatrixXd& V,
        const Eigen::MatrixXd& Vt, const std::vector<MMCVID>& constraintSet,
        const Eigen::VectorXd& multipliers,
        const std::vector<Eigen::Vector2d>& MMDistCoord,
        const std::vector<Eigen::Matrix<double, 3, 2>>& MMTanBasis,
        double& Ef, double eps2, double coef);

    static void augmentFrictionGradient(const Eigen::MatrixXd& V,
        const Eigen::MatrixXd& Vt, const std::vector<MMCVID>& constraintSet,
        const Eigen::VectorXd& multipliers,
        const std::vector<Eigen::Vector2d>& MMDistCoord,
        const std::vector<Eigen::Matrix<double, 3, 2>>& MMTanBasis,
        Eigen::VectorXd& grad_inc, double eps2, double coef);

    static void augmentFrictionHessian(const Mesh<dim>& mesh,
        const Eigen::MatrixXd& Vt, const std::vector<MMCVID>& constraintSet,
        const Eigen::VectorXd& multipliers,
        const std::vector<Eigen::Vector2d>& MMDistCoord,
        const std::vector<Eigen::Matrix<double, 3, 2>>& MMTanBasis,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* H_inc,
        double eps2, double coef, bool projectDBC);

    static void augmentFrictionHessian(const Mesh<dim>& mesh,
        const Eigen::MatrixXd& Vt, const std::vector<MMCVID>& constraintSet,
        const Eigen::VectorXd& multipliers,
        const std::vector<Eigen::Vector2d>& MMDistCoord,
        const std::vector<Eigen::Matrix<double, 3, 2>>& MMTanBasis,
        const std::function<void(size_t, size_t, double)>& addCoeff,
        double eps2, double coef, bool projectDBC);

    static void augmentParaEEGradient(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& paraEEMMCVIDSet,
        const std::vector<std::pair<int, int>>& paraEEeIeJSet,
        Eigen::VectorXd& grad_inc, double dHat, double coef);

    static void augmentParaEEHessian(const Mesh<dim>& mesh,
        const std::vector<MMCVID>& paraEEMMCVIDSet,
        const std::vector<std::pair<int, int>>& paraEEeIeJSet,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* H_inc,
        double dHat, double coef, bool projectDBC);

    static bool checkEdgeTriIntersection(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh);
    static bool checkEdgeTriIntersectionIfAny(const Mesh<dim>& mesh,
        const SpatialHash<dim>& sh);

    /**
     * @brief Determine if the mesh is in an intersected state.
     *
     * Check if the mesh enters an intersected state along the linear
     * trajectories from V0 to mesh.V.
     *
     * @param[in]  mesh    Mesh to check for intersections.
     * @param[in]  V0      Mesh vertices at previous iteration.
     * @param[in]  method  Method of exact CCD to use.
     *
     * @returns  True if the mesh is intersected.
     */
    static bool isIntersected(
        const Mesh<dim>& mesh,
        const Eigen::MatrixXd& V0,
        const ccd::CCDMethod method);
};

} // namespace IPC

#endif /* SelfCollisionHandler_hpp */
