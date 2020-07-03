//
//  FrictionUtils.hpp
//  IPC
//
//  Created by Minchen Li on 9/26/19.
//

#ifndef FrictionUtils_hpp
#define FrictionUtils_hpp

#include "Types.hpp"

#include "MeshCollisionUtils.hpp"

#include <Eigen/Eigen>

#include <iostream>
#include <array>

namespace IPC {

// Point - Triangle

inline void computeTangentBasis_PT(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    Eigen::Matrix<double, 3, 2>& basis)
{
    Eigen::Vector3d v12 = (v2 - v1).transpose();
    basis.col(0) = v12.normalized();
    basis.col(1) = v12.cross((v3 - v1).transpose()).cross(v12).normalized();
}

inline void computeClosestPoint_PT(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    Eigen::Vector2d& beta)
{
    Eigen::Matrix<double, 2, 3> basis;
    basis.row(0) = v2 - v1;
    basis.row(1) = v3 - v1;
    beta = (basis * basis.transpose()).ldlt().solve(basis * (v0 - v1).transpose());
}

inline void computeRelDX_PT(
    const Eigen::RowVector3d& dx0,
    const Eigen::RowVector3d& dx1,
    const Eigen::RowVector3d& dx2,
    const Eigen::RowVector3d& dx3,
    double beta1, double beta2,
    Eigen::RowVector3d& relDX)
{
    relDX = dx0 - (dx1 + beta1 * (dx2 - dx1) + beta2 * (dx3 - dx1));
}

inline void liftRelDXTanToMesh_PT(
    const Eigen::Vector2d& relDXTan,
    const Eigen::Matrix<double, 3, 2>& basis,
    double beta1, double beta2,
    Eigen::Matrix<double, 12, 1>& TTTDX)
{
    TTTDX.template segment<3>(0) = basis * relDXTan;
    TTTDX.template segment<3>(3) = (-1 + beta1 + beta2) * TTTDX.template segment<3>(0);
    TTTDX.template segment<3>(6) = -beta1 * TTTDX.template segment<3>(0);
    TTTDX.template segment<3>(9) = -beta2 * TTTDX.template segment<3>(0);
}

inline void computeTTT_PT(
    const Eigen::Matrix<double, 3, 2>& basis,
    double beta1, double beta2,
    Eigen::Matrix<double, 12, 12>& TTT)
{
    Eigen::Matrix<double, 2, 12> TT;
    TT.template block<2, 3>(0, 0) = basis.transpose();
    TT.template block<2, 3>(0, 3) = (-1 + beta1 + beta2) * basis.transpose();
    TT.template block<2, 3>(0, 6) = -beta1 * basis.transpose();
    TT.template block<2, 3>(0, 9) = -beta2 * basis.transpose();
    TTT = TT.transpose() * TT;
}

// Edge - Edge

inline void computeTangentBasis_EE(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    Eigen::Matrix<double, 3, 2>& basis)
{
    Eigen::Vector3d v01 = (v1 - v0).transpose();
    basis.col(0) = v01.normalized();
    basis.col(1) = v01.cross((v3 - v2).transpose()).cross(v01).normalized();
}

inline void computeClosestPoint_EE(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    Eigen::Vector2d& gamma)
{
    Eigen::RowVector3d e20 = v0 - v2;
    Eigen::RowVector3d e01 = v1 - v0;
    Eigen::RowVector3d e23 = v3 - v2;

    Eigen::Matrix2d coefMtr;
    coefMtr(0, 0) = e01.squaredNorm();
    coefMtr(0, 1) = coefMtr(1, 0) = -e23.dot(e01);
    coefMtr(1, 1) = e23.squaredNorm();

    Eigen::Vector2d rhs;
    rhs[0] = -e20.dot(e01);
    rhs[1] = e20.dot(e23);

    gamma = coefMtr.ldlt().solve(rhs);
}

inline void computeRelDX_EE(
    const Eigen::RowVector3d& dx0,
    const Eigen::RowVector3d& dx1,
    const Eigen::RowVector3d& dx2,
    const Eigen::RowVector3d& dx3,
    double gamma1, double gamma2,
    Eigen::RowVector3d& relDX)
{
    relDX = dx0 + gamma1 * (dx1 - dx0) - (dx2 + gamma2 * (dx3 - dx2));
}

inline void liftRelDXTanToMesh_EE(
    const Eigen::Vector2d& relDXTan,
    const Eigen::Matrix<double, 3, 2>& basis,
    double gamma1, double gamma2,
    Eigen::Matrix<double, 12, 1>& TTTDX)
{
    Eigen::Vector3d relDXTan3D = basis * relDXTan;
    TTTDX.template segment<3>(0) = (1.0 - gamma1) * relDXTan3D;
    TTTDX.template segment<3>(3) = gamma1 * relDXTan3D;
    TTTDX.template segment<3>(6) = (gamma2 - 1.0) * relDXTan3D;
    TTTDX.template segment<3>(9) = -gamma2 * relDXTan3D;
}

inline void computeTTT_EE(
    const Eigen::Matrix<double, 3, 2>& basis,
    double gamma1, double gamma2,
    Eigen::Matrix<double, 12, 12>& TTT)
{
    Eigen::Matrix<double, 2, 12> TT;
    TT.template block<2, 3>(0, 0) = (1.0 - gamma1) * basis.transpose();
    TT.template block<2, 3>(0, 3) = gamma1 * basis.transpose();
    TT.template block<2, 3>(0, 6) = (gamma2 - 1.0) * basis.transpose();
    TT.template block<2, 3>(0, 9) = -gamma2 * basis.transpose();
    TTT = TT.transpose() * TT;
}

// Point - Edge

inline void computeTangentBasis_PE(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    Eigen::Matrix<double, 3, 2>& basis)
{
    Eigen::Vector3d v12 = (v2 - v1).transpose();
    basis.col(0) = v12.normalized();
    basis.col(1) = v12.cross((v0 - v1).transpose()).normalized();
}

inline void computeClosestPoint_PE(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    double& yita)
{
    Eigen::RowVector3d e12 = v2 - v1;
    yita = (v0 - v1).dot(e12) / e12.squaredNorm();
}

inline void computeRelDX_PE(
    const Eigen::RowVector3d& dx0,
    const Eigen::RowVector3d& dx1,
    const Eigen::RowVector3d& dx2,
    double yita,
    Eigen::RowVector3d& relDX)
{
    relDX = dx0 - (dx1 + yita * (dx2 - dx1));
}

inline void liftRelDXTanToMesh_PE(
    const Eigen::Vector2d& relDXTan,
    const Eigen::Matrix<double, 3, 2>& basis,
    double yita,
    Eigen::Matrix<double, 9, 1>& TTTDX)
{
    TTTDX.template segment<3>(0) = basis * relDXTan;
    TTTDX.template segment<3>(3) = (yita - 1.0) * TTTDX.template segment<3>(0);
    TTTDX.template segment<3>(6) = -yita * TTTDX.template segment<3>(0);
}

inline void computeTTT_PE(
    const Eigen::Matrix<double, 3, 2>& basis,
    double yita,
    Eigen::Matrix<double, 9, 9>& TTT)
{
    Eigen::Matrix<double, 2, 9> TT;
    TT.template block<2, 3>(0, 0) = basis.transpose();
    TT.template block<2, 3>(0, 3) = (yita - 1.0) * basis.transpose();
    TT.template block<2, 3>(0, 6) = -yita * basis.transpose();
    TTT = TT.transpose() * TT;
}

// Point - Point

inline void computeTangentBasis_PP(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    Eigen::Matrix<double, 3, 2>& basis)
{
    Eigen::RowVector3d v01 = v1 - v0;
    Eigen::RowVector3d xCross = Eigen::RowVector3d::UnitX().cross(v01);
    Eigen::RowVector3d yCross = Eigen::RowVector3d::UnitY().cross(v01);
    if (xCross.squaredNorm() > yCross.squaredNorm()) {
        basis.col(0) = xCross.normalized().transpose();
        basis.col(1) = v01.cross(xCross).normalized().transpose();
    }
    else {
        basis.col(0) = yCross.normalized().transpose();
        basis.col(1) = v01.cross(yCross).normalized().transpose();
    }
}

inline void computeRelDX_PP(
    const Eigen::RowVector3d& dx0,
    const Eigen::RowVector3d& dx1,
    Eigen::RowVector3d& relDX)
{
    relDX = dx0 - dx1;
}

inline void liftRelDXTanToMesh_PP(
    const Eigen::Vector2d& relDXTan,
    const Eigen::Matrix<double, 3, 2>& basis,
    Eigen::Matrix<double, 6, 1>& TTTDX)
{
    TTTDX.template segment<3>(0) = basis * relDXTan;
    TTTDX.template segment<3>(3) = -TTTDX.template segment<3>(0);
}

inline void computeTTT_PP(
    const Eigen::Matrix<double, 3, 2>& basis,
    Eigen::Matrix<double, 6, 6>& TTT)
{
    Eigen::Matrix<double, 2, 6> TT;
    TT.template block<2, 3>(0, 0) = basis.transpose();
    TT.template block<2, 3>(0, 3) = -basis.transpose();
    TTT = TT.transpose() * TT;
}

// static friction clamping model
// C0 clamping
inline void f0_SF_C0(double x2, double eps_f, double& f0)
{
    f0 = x2 / (2.0 * eps_f) + eps_f / 2.0;
}

inline void f1_SF_div_relDXNorm_C0(double eps_f, double& result)
{
    result = 1.0 / eps_f;
}

inline void f2_SF_C0(double eps_f, double& f2)
{
    f2 = 1.0 / eps_f;
}

// C1 clamping
inline void f0_SF_C1(double x2, double eps_f, double& f0)
{
    f0 = x2 * (-std::sqrt(x2) / 3.0 + eps_f) / (eps_f * eps_f) + eps_f / 3.0;
}

inline void f1_SF_div_relDXNorm_C1(double x2, double eps_f, double& result)
{
    result = (-std::sqrt(x2) + 2.0 * eps_f) / (eps_f * eps_f);
}

inline void f2_SF_C1(double x2, double eps_f, double& f2)
{
    f2 = 2.0 * (eps_f - std::sqrt(x2)) / (eps_f * eps_f);
}

// C2 clamping
inline void f0_SF_C2(double x2, double eps_f, double& f0)
{
    f0 = x2 * (0.25 * x2 - (std::sqrt(x2) - 1.5 * eps_f) * eps_f) / (eps_f * eps_f * eps_f) + eps_f / 4.0;
}

inline void f1_SF_div_relDXNorm_C2(double x2, double eps_f, double& result)
{
    result = (x2 - (3.0 * std::sqrt(x2) - 3.0 * eps_f) * eps_f) / (eps_f * eps_f * eps_f);
}

inline void f2_SF_C2(double x2, double eps_f, double& f2)
{
    f2 = 3.0 * (x2 - (2.0 * std::sqrt(x2) - eps_f) * eps_f) / (eps_f * eps_f * eps_f);
}

// interfaces
inline void f0_SF(double relDXSqNorm, double eps_f, double& f0)
{
#if (SFCLAMPING_ORDER == 0)
    f0_SF_C0(relDXSqNorm, eps_f, f0);
#elif (SFCLAMPING_ORDER == 1)
    f0_SF_C1(relDXSqNorm, eps_f, f0);
#elif (SFCLAMPING_ORDER == 2)
    f0_SF_C2(relDXSqNorm, eps_f, f0);
#endif
}

inline void f1_SF_div_relDXNorm(double relDXSqNorm, double eps_f, double& result)
{
#if (SFCLAMPING_ORDER == 0)
    f1_SF_div_relDXNorm_C0(eps_f, result);
#elif (SFCLAMPING_ORDER == 1)
    f1_SF_div_relDXNorm_C1(relDXSqNorm, eps_f, result);
#elif (SFCLAMPING_ORDER == 2)
    f1_SF_div_relDXNorm_C2(relDXSqNorm, eps_f, result);
#endif
}

inline void f2_SF(double relDXSqNorm, double eps_f, double& f2)
{
#if (SFCLAMPING_ORDER == 0)
    f2_SF_C0(eps_f, f2);
#elif (SFCLAMPING_ORDER == 1)
    f2_SF_C1(relDXSqNorm, eps_f, f2);
#elif (SFCLAMPING_ORDER == 2)
    f2_SF_C2(relDXSqNorm, eps_f, f2);
#endif
}

} // namespace IPC

#endif /* FrictionUtils_hpp */
