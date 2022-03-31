//
//  HalfSpace.cpp
//  IPC
//
//  Created by Minchen Li on 10/22/18.
//

#include "HalfSpace.hpp"

#include "BarrierFunctions.hpp"
#include "IglUtils.hpp"

#ifdef USE_TBB
#include <tbb/parallel_for.h>
#endif

namespace IPC {

template <int dim>
HalfSpace<dim>::HalfSpace(const Eigen::Matrix<double, dim, 1>& p_origin,
    const Eigen::Matrix<double, dim, 1>& p_normal,
    const Eigen::Matrix<double, dim, 1>& p_velocitydt,
    double p_friction)
{
    init(p_origin, p_normal, p_velocitydt, p_friction);
}

template <int dim>
HalfSpace<dim>::HalfSpace(double p_Y, double p_friction)
{
    Eigen::Matrix<double, dim, 1> p_origin, p_normal, p_velocitydt;
    p_origin.setZero();
    p_origin[1] = p_Y;
    p_normal.setZero();
    p_normal[1] = 1.0;
    p_velocitydt.setZero();

    init(p_origin, p_normal, p_velocitydt, p_friction);
}

template <int dim>
void HalfSpace<dim>::init(const Eigen::Matrix<double, dim, 1>& p_origin,
    const Eigen::Matrix<double, dim, 1>& p_normal,
    const Eigen::Matrix<double, dim, 1>& p_velocitydt,
    double p_friction)
{
    Base::origin = p_origin;
    normal = p_normal;
    normal.normalize();
    Base::velocitydt = p_velocitydt;
    D = -normal.dot(Base::origin);
    nnT = normal * normal.transpose();
    if constexpr (dim == 3) {
        Eigen::RowVector3d xCross = Eigen::RowVector3d::UnitX().cross(normal);
        Eigen::RowVector3d yCross = Eigen::RowVector3d::UnitY().cross(normal);
        if (xCross.squaredNorm() > yCross.squaredNorm()) {
            basis.col(0) = xCross.normalized().transpose();
            basis.col(1) = normal.cross(xCross).normalized().transpose();
        }
        else {
            basis.col(0) = yCross.normalized().transpose();
            basis.col(1) = normal.cross(yCross).normalized().transpose();
        }
    }

    Eigen::Matrix<double, dim, 1> defaultN;
    defaultN.setZero();
    defaultN[1] = 1.0;
    double rotAngle = std::acos(std::max(-1.0, std::min(1.0, normal.dot(defaultN))));
    if (rotAngle < 1.0e-3) {
        rotMtr.setIdentity();
    }
    else {
#if (DIM == 3)
        rotMtr = Eigen::AngleAxisd(rotAngle, (defaultN.cross(normal)).normalized());
#else
        rotMtr = Eigen::Rotation2Dd(rotAngle);
#endif
    }

    Base::friction = p_friction;

    initRenderingData();
}

template <int dim>
void HalfSpace<dim>::updateConstraints_QP(const Mesh<dim>& mesh,
    const std::vector<int>& activeSet,
    std::vector<Eigen::Triplet<double>>& A_triplet, Eigen::VectorXd& l) const
{
    A_triplet.reserve(A_triplet.size() + activeSet.size() * dim);
    int constraintI = l.size();
    for (const auto& vI : activeSet) {
        A_triplet.emplace_back(constraintI, vI * dim, normal[0]);
        A_triplet.emplace_back(constraintI, vI * dim + 1, normal[1]);
        if constexpr (dim == 3) {
            A_triplet.emplace_back(constraintI, vI * dim + 2, normal[2]);
        }
        ++constraintI;
    }

    Base::evaluateConstraintsQP(mesh, activeSet, l, -1.0);
}

template <int dim>
void HalfSpace<dim>::evaluateConstraint(const Mesh<dim>& mesh,
    int vI, double& val, double coef) const
{
    double dist = normal.transpose().dot(mesh.V.row(vI)) + D;
    val = coef * dist * dist;
}
template <int dim>
void HalfSpace<dim>::evaluateConstraintQP(const Mesh<dim>& mesh,
    int vI, double& val, double coef) const
{
    double dist = normal.transpose().dot(mesh.V.row(vI)) + D;
    val = coef * dist;
}

template <int dim>
void HalfSpace<dim>::leftMultiplyConstraintJacobianT(const Mesh<dim>& mesh,
    const std::vector<int>& activeSet,
    const Eigen::VectorXd& input,
    Eigen::VectorXd& output_incremental,
    double coef) const
{
    assert(input.size() == activeSet.size());
    assert(output_incremental.size() == mesh.V.rows() * dim);

#ifdef USE_TBB
    tbb::parallel_for(0, (int)activeSet.size(), 1, [&](int cI)
#else
    for (int cI = 0; cI < activeSet.size(); ++cI)
#endif
        {
            int vI = activeSet[cI];
            double dist = normal.transpose().dot(mesh.V.row(vI)) + D;
            output_incremental.segment<dim>(vI * dim) += coef * input[cI] * 2.0 * dist * normal;
        }
#ifdef USE_TBB
    );
#endif
}
template <int dim>
void HalfSpace<dim>::leftMultiplyConstraintJacobianTQP(const Mesh<dim>& mesh,
    const std::vector<int>& activeSet,
    const Eigen::VectorXd& input,
    Eigen::VectorXd& output_incremental,
    double coef) const
{
    assert(input.size() == activeSet.size());
    assert(output_incremental.size() == mesh.V.rows() * dim);

#ifdef USE_TBB
    tbb::parallel_for(0, (int)activeSet.size(), 1, [&](int cI)
#else
    for (int cI = 0; cI < activeSet.size(); ++cI)
#endif
        {
            int vI = activeSet[cI];
            output_incremental.segment<dim>(vI * dim) += coef * input[cI] * normal;
        }
#ifdef USE_TBB
    );
#endif
}

template <int dim>
void HalfSpace<dim>::augmentIPHessian(const Mesh<dim>& mesh,
    const std::vector<int>& activeSet,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* mtr_incremental,
    double dHat, double coef, bool projectDBC) const
{
#ifdef USE_TBB
    tbb::parallel_for(0, (int)activeSet.size(), 1, [&](int constraintI)
#else
    for (int constraintI = 0; constraintI < activeSet.size(); ++constraintI)
#endif
        {
            int vI = activeSet[constraintI];
            if (!mesh.isDBCVertex(vI) || !projectDBC) {
                double constraintVal = 0;
                evaluateConstraint(mesh, vI, constraintVal);
                double g_b = 0;
                compute_g_b(constraintVal, dHat, g_b);
                double H_b = 0;
                compute_H_b(constraintVal, dHat, H_b);
                double param = 4.0 * H_b * constraintVal + 2.0 * g_b;
                if (param > 0.0) {
                    Eigen::Matrix<double, dim, dim> H = coef * param * nnT;
                    // signed distance formulation for QP:
                    // Eigen::Matrix<double, dim, dim> H = coef * H_b * nnT;

                    int rowStartI = vI * dim;
                    mtr_incremental->addCoeff(rowStartI, rowStartI, H(0, 0));
                    mtr_incremental->addCoeff(rowStartI, rowStartI + 1, H(0, 1));
                    mtr_incremental->addCoeff(rowStartI + 1, rowStartI, H(1, 0));
                    mtr_incremental->addCoeff(rowStartI + 1, rowStartI + 1, H(1, 1));
                    if constexpr (dim == 3) {
                        mtr_incremental->addCoeff(rowStartI, rowStartI + 2, H(0, 2));
                        mtr_incremental->addCoeff(rowStartI + 1, rowStartI + 2, H(1, 2));

                        mtr_incremental->addCoeff(rowStartI + 2, rowStartI, H(2, 0));
                        mtr_incremental->addCoeff(rowStartI + 2, rowStartI + 1, H(2, 1));
                        mtr_incremental->addCoeff(rowStartI + 2, rowStartI + 2, H(2, 2));
                    }
                }
            }
        }
#ifdef USE_TBB
    );
#endif
}

template <int dim>
void HalfSpace<dim>::filterSearchDir_QP(const Mesh<dim>& mesh,
    Eigen::VectorXd& searchDir,
    std::vector<int>& activeSet_next)
{
    Eigen::VectorXd constraintVal;
    Base::evaluateConstraintsQP_all(mesh, constraintVal, -1.0);
    for (int vI = 0; vI < mesh.V.rows(); vI++) {
        double coef = normal.dot(searchDir.segment<dim>(vI * dim));
        if (coef < 0.0) { // if going towards the halfSpace
            double maxStepSize = constraintVal[vI] / coef;
            if (maxStepSize < 1.0) { // includes points inside halfSpace (maxStepSize <= 0)
                if (mesh.m_isBoundaryVert[vI]) {
                    activeSet_next.emplace_back(vI);
                }
                // project out normal search direction so that stepSize won't be limited:
                searchDir.segment<dim>(vI * dim) -= (1.0 - maxStepSize) * coef * normal;
                // NOTE: clamp search direction colinearly will fail the line search
            }
        }
    }

    std::sort(activeSet_next.begin(), activeSet_next.end());
    activeSet_next.erase(std::unique(activeSet_next.begin(), activeSet_next.end()), activeSet_next.end());
}

template <int dim>
void HalfSpace<dim>::largestFeasibleStepSize(const Mesh<dim>& mesh,
    const Eigen::VectorXd& searchDir,
    double slackness,
    std::vector<int>& activeSet_next,
    double& stepSize)
{
    Eigen::VectorXd maxStepSizes(mesh.SVI.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SVI.size(), 1, [&](int svI)
#else
    for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
        {
            maxStepSizes[svI] = 1.0;
            int vI = mesh.SVI[svI];
            if (!mesh.isDBCVertex(vI)) {
                double coef = normal.dot(searchDir.segment<dim>(vI * dim));
                if (coef < 0.0) { // if going towards the halfSpace
                    double dist = normal.transpose().dot(mesh.V.row(vI)) + D;
                    maxStepSizes[svI] = -dist / coef * slackness;
                }
            }
        }
#ifdef USE_TBB
    );
#endif
    stepSize = std::min(stepSize, maxStepSizes.minCoeff());
}

template <int dim>
void HalfSpace<dim>::computeFrictionEnergy(const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& Vt, const std::vector<int>& activeSet,
    const Eigen::VectorXd& multipliers,
    double& Ef, double eps2, double coef) const
{
    assert(multipliers.size() == activeSet.size());
    double eps = std::sqrt(eps2);

    if constexpr (dim == 3) {
        // TODO: parallelize
        Ef = 0.0;
        int contactPairI = 0;
        for (const auto& vI : activeSet) {
            Eigen::Matrix<double, dim, 1> VDiff = (V.row(vI) - Vt.row(vI)).transpose();
            VDiff -= Base::velocitydt;
            Eigen::Matrix<double, dim, 1> VProj = VDiff - VDiff.dot(normal) * normal;
            double VProjMag2 = VProj.squaredNorm();
            if (VProjMag2 > eps2) {
                Ef += Base::friction * multipliers[contactPairI] * (std::sqrt(VProjMag2) - eps * 0.5);
            }
            else {
                Ef += Base::friction * multipliers[contactPairI] * VProjMag2 / eps * 0.5;
            }
            ++contactPairI;
        }
        Ef *= coef;
    }
}
template <int dim>
void HalfSpace<dim>::augmentFrictionGradient(const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& Vt, const std::vector<int>& activeSet,
    const Eigen::VectorXd& multipliers,
    Eigen::VectorXd& grad_inc, double eps2, double coef) const
{
    assert(multipliers.size() == activeSet.size());
    double eps = std::sqrt(eps2);

    if constexpr (dim == 3) {
        // TODO: parallelize
        int contactPairI = 0;
        for (const auto& vI : activeSet) {
            Eigen::Matrix<double, dim, 1> VDiff = (V.row(vI) - Vt.row(vI)).transpose();
            VDiff -= Base::velocitydt;
            Eigen::Matrix<double, dim, 1> VProj = VDiff - VDiff.dot(normal) * normal;
            double VProjMag2 = VProj.squaredNorm();
            if (VProjMag2 > eps2) {
                grad_inc.template segment<dim>(vI * dim) += coef * Base::friction * multipliers[contactPairI] / std::sqrt(VProjMag2) * VProj;
            }
            else {
                grad_inc.template segment<dim>(vI * dim) += coef * Base::friction * multipliers[contactPairI] / eps * VProj;
            }
            ++contactPairI;
        }
    }
}
template <int dim>
void HalfSpace<dim>::augmentFrictionHessian(const Mesh<dim>& mesh,
    const Eigen::MatrixXd& Vt, const std::vector<int>& activeSet,
    const Eigen::VectorXd& multipliers,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* H_inc,
    double eps2, double coef, bool projectDBC) const
{
    assert(multipliers.size() == activeSet.size());
    double eps = std::sqrt(eps2);

    // TODO: parallelize
    int contactPairI = 0;
    for (const auto& vI : activeSet) {
        if (projectDBC && mesh.isDBCVertex(vI)) {
            continue;
        }

        double multiplier_vI = coef * Base::friction * multipliers[contactPairI];
        Eigen::Matrix<double, dim, dim> H_vI;

        Eigen::Matrix<double, dim, 1> VDiff = (mesh.V.row(vI) - Vt.row(vI)).transpose();
        VDiff -= Base::velocitydt;
        Eigen::Matrix<double, dim, 1> VProj = VDiff - VDiff.dot(normal) * normal;
        double VProjMag2 = VProj.squaredNorm();
        if (VProjMag2 > eps2) {
            double VProjMag = std::sqrt(VProjMag2);

            H_vI = (VProj * (-multiplier_vI / VProjMag2 / VProjMag)) * VProj.transpose();
            H_vI += (Eigen::Matrix<double, dim, dim>::Identity() - normal * normal.transpose()) * (multiplier_vI / VProjMag);

            IglUtils::makePD(H_vI);
        }
        else {
            H_vI = (Eigen::Matrix<double, dim, dim>::Identity() - normal * normal.transpose()) * (multiplier_vI / eps);
            // already SPD
        }

        int startInd = vI * dim;
        H_inc->addCoeff(startInd, startInd, H_vI(0, 0));
        H_inc->addCoeff(startInd, startInd + 1, H_vI(0, 1));
        H_inc->addCoeff(startInd + 1, startInd, H_vI(1, 0));
        H_inc->addCoeff(startInd + 1, startInd + 1, H_vI(1, 1));
        if constexpr (dim == 3) {
            H_inc->addCoeff(startInd, startInd + 2, H_vI(0, 2));
            H_inc->addCoeff(startInd + 1, startInd + 2, H_vI(1, 2));

            H_inc->addCoeff(startInd + 2, startInd, H_vI(2, 0));
            H_inc->addCoeff(startInd + 2, startInd + 1, H_vI(2, 1));
            H_inc->addCoeff(startInd + 2, startInd + 2, H_vI(2, 2));
        }

        ++contactPairI;
    }
}

template <int dim>
void HalfSpace<dim>::computeRelDX(const Eigen::RowVector3d& relDX3D, Eigen::Vector2d& relDX) const
{
    relDX = (relDX3D * basis).transpose();
}

template <int dim>
void HalfSpace<dim>::move(const Eigen::Matrix<double, dim, 1>& deltaX,
    const Mesh<dim>& mesh, const SpatialHash<dim>& sh,
    double slackness, double& stepSizeLeft)
{
    Eigen::VectorXd maxStepSizes(mesh.SVI.size());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.SVI.size(), 1, [&](int svI)
#else
    for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
        {
            maxStepSizes[svI] = 1.0;
            int vI = mesh.SVI[svI];
            double coef = normal.dot(-deltaX);
            if (coef < 0.0) { // if going towards the object
                double dist = normal.transpose().dot(mesh.V.row(vI)) + D;
                maxStepSizes[svI] = -dist / coef * slackness;
            }
        }
#ifdef USE_TBB
    );
#endif
    double stepSize = std::min(1.0, maxStepSizes.minCoeff());

    init(Base::origin + stepSize * deltaX, normal, Base::velocitydt, Base::friction);

    stepSizeLeft = 1.0 - stepSize;
}

template <int dim>
void HalfSpace<dim>::initRenderingData(double extensionScale)
{
    double size = extensionScale * 10.0;
    int elemAmt = 200;
    double spacing = size / std::sqrt(elemAmt / 2.0);
    assert(size >= spacing);
    int gridSize = static_cast<int>(size / spacing) + 1;
    spacing = size / (gridSize - 1);

    Base::V.resize(gridSize * gridSize, 3);
    for (int rowI = 0; rowI < gridSize; rowI++) {
        for (int colI = 0; colI < gridSize; colI++) {
            int vI = rowI * gridSize + colI;
            Eigen::Matrix<double, 1, dim> coord;
            coord[0] = spacing * colI - size / 2.0;
            coord[1] = 0.0;
            if constexpr (dim == 3) {
                coord[2] = spacing * rowI - size / 2.0;
            }
            Base::V.row(vI).leftCols(dim) = (rotMtr * coord.transpose() + Base::origin).transpose();
            if constexpr (dim == 2) {
                // TODO: debug
                Base::V(vI, 2) = 0.0;
            }
        }
    }

    Base::F.resize((gridSize - 1) * (gridSize - 1) * 2, 3);
    for (int rowI = 0; rowI < gridSize - 1; rowI++) {
        for (int colI = 0; colI < gridSize - 1; colI++) {
            int squareI = rowI * (gridSize - 1) + colI;
            Base::F.row(squareI * 2) = Eigen::Vector3i(rowI * gridSize + colI,
                (rowI + 1) * gridSize + colI,
                (rowI + 1) * gridSize + colI + 1);
            Base::F.row(squareI * 2 + 1) = Eigen::Vector3i(rowI * gridSize + colI,
                (rowI + 1) * gridSize + colI + 1,
                rowI * gridSize + colI + 1);
        }
    }
}

template class HalfSpace<DIM>;
} // namespace IPC
