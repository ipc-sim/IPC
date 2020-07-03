//
//  HalfSpace.hpp
//  IPC
//
//  Created by Minchen Li on 10/22/18.
//

#ifndef HalfSpace_hpp
#define HalfSpace_hpp

#include "CollisionObject.h"

#include "Energy.hpp"

namespace IPC {

template <int dim>
class HalfSpace : public CollisionObject<dim> {
    typedef CollisionObject<dim> Base;

protected:
    Eigen::Matrix<double, dim, 1> normal;
    double D;
    // normal[0] x + normal[1] y + normal[2] z + D = 0
    Eigen::Matrix<double, dim, dim> nnT;
    Eigen::Matrix<double, 3, 2> basis;

    Eigen::Matrix<double, dim, dim> rotMtr; // for visualization

public:
    HalfSpace(const Eigen::Matrix<double, dim, 1>& origin,
        const Eigen::Matrix<double, dim, 1>& normal,
        const Eigen::Matrix<double, dim, 1>& p_velocitydt,
        double p_friction);

    HalfSpace(double p_Y, double p_friction);

    void init(const Eigen::Matrix<double, dim, 1>& origin,
        const Eigen::Matrix<double, dim, 1>& normal,
        const Eigen::Matrix<double, dim, 1>& p_velocitydt,
        double p_friction);

public:
    virtual void updateConstraints_QP(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        std::vector<Eigen::Triplet<double>>& A_triplet, Eigen::VectorXd& l) const override;

    virtual void evaluateConstraint(const Mesh<dim>& mesh,
        int vI, double& val, double coef = 1.0) const override;
    virtual void evaluateConstraintQP(const Mesh<dim>& mesh,
        int vI, double& val, double coef = 1.0) const override;

    virtual void leftMultiplyConstraintJacobianT(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output_incremental,
        double coef = 1.0) const override;
    virtual void leftMultiplyConstraintJacobianTQP(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        const Eigen::VectorXd& input,
        Eigen::VectorXd& output,
        double coef = 1.0) const override;
    virtual void augmentIPHessian(const Mesh<dim>& mesh,
        const std::vector<int>& activeSet,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* mtr_incremental,
        double dHat, double coef = 1.0, bool projectDBC = true) const override;

    virtual void filterSearchDir_QP(const Mesh<dim>& mesh,
        Eigen::VectorXd& searchDir,
        std::vector<int>& activeSet_next) override;
    virtual void largestFeasibleStepSize(const Mesh<dim>& mesh,
        const Eigen::VectorXd& searchDir,
        double slackness,
        std::vector<int>& activeSet_next,
        double& stepSize) override;

    virtual void computeFrictionEnergy(const Eigen::MatrixXd& V,
        const Eigen::MatrixXd& Vt, const std::vector<int>& activeSet,
        const Eigen::VectorXd& multipliers,
        double& Ef, double eps2, double coef) const override;
    virtual void augmentFrictionGradient(const Eigen::MatrixXd& V,
        const Eigen::MatrixXd& Vt, const std::vector<int>& activeSet,
        const Eigen::VectorXd& multipliers,
        Eigen::VectorXd& grad_inc, double eps2, double coef) const override;
    virtual void augmentFrictionHessian(const Mesh<dim>& mesh,
        const Eigen::MatrixXd& Vt, const std::vector<int>& activeSet,
        const Eigen::VectorXd& multipliers,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* H_inc,
        double eps2, double coef, bool projectDBC = true) const override;

    virtual void computeRelDX(const Eigen::RowVector3d& relDX3D, Eigen::Vector2d& relDX) const override;

    virtual void move(const Eigen::Matrix<double, dim, 1>& deltaX,
        const Mesh<dim>& mesh, const SpatialHash<dim>& sh,
        double slackness, double& stepSizeLeft) override;

    virtual void initRenderingData(double extensionScale = 1.0) override;
};

} // namespace IPC

#endif /* HalfSpace_hpp */
