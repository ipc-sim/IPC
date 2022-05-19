//
//  HyperFoamEnergy.cpp
//  IPC
//
//  Created by Ehsan Haghighat 5/18/22.
//

#include "ipc/Energy/Physics_Elasticity/HyperFoamEnergy.hpp"
#include "ipc/Utils/IglUtils.hpp"

namespace IPC {

template <int dim>
void HyperFoamEnergy<dim>::computeEnergyVal(const Mesh<dim>& data, int redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    double& energyVal) const
{
    Base::computeEnergyValBySVD(data, redoSVD, svd, F, coef, energyVal);
}

template <int dim>
void HyperFoamEnergy<dim>::computeGradient(const Mesh<dim>& data, bool redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    Eigen::VectorXd& gradient,
    bool projectDBC) const
{
    Base::computeGradientByPK(data, redoSVD, svd, F, coef, gradient, projectDBC);
}

template <int dim>
void HyperFoamEnergy<dim>::computeHessian(const Mesh<dim>& data, bool redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver,
    bool projectSPD,
    bool projectDBC) const
{
    Base::computeHessianByPK(data, redoSVD, svd, F, coef, linSysSolver, projectSPD, projectDBC);
}

template <int dim>
void HyperFoamEnergy<dim>::getEnergyValPerElem(const Mesh<dim>& data,
    Eigen::VectorXd& energyValPerElem,
    bool uniformWeight) const
{
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>> svd(data.F.rows());
    std::vector<Eigen::Matrix<double, dim, dim>> F(data.F.rows());
    Energy<dim>::getEnergyValPerElemBySVD(data, true, svd, F, energyValPerElem, uniformWeight);
}

template <int dim>
void HyperFoamEnergy<dim>::compute_E(const Eigen::Matrix<double, dim, 1>& singularValues,
    const MaterialProps& mat,
    double& E) const
{
    if (mat.u == 0.0 && mat.lambda == 0.0) {
        E = 0.;
        return;
    }

    const double sigmaASum = singularValues.array().pow(mat.alpha).sum();
    const double sigmaProd = singularValues.prod();
    const double sigmaProd_ab = std::pow(sigmaProd, -mat.alpha*mat.beta);
    
    const double coef = 2.0*mat.u / std::pow(mat.alpha, 2);
    E = coef * ((sigmaASum - dim) + (sigmaProd_ab - 1.) / mat.beta);
}

template <int dim>
void HyperFoamEnergy<dim>::compute_dE_div_dsigma(const Eigen::Matrix<double, dim, 1>& singularValues,
    const MaterialProps& mat,
    Eigen::Matrix<double, dim, 1>& dE_div_dsigma) const
{
    if (mat.u == 0.0 && mat.lambda == 0.0) {
        dE_div_dsigma.setZero();
        return;
    }

    const double sigmaProd_ab = std::pow(singularValues.prod(), -mat.alpha*mat.beta);
    const double coef = 2.0*mat.u/mat.alpha; 

    const auto invSingularValues = singularValues.array().inverse();
    const auto singulatValues_powa = singularValues.array().pow(mat.alpha);
    dE_div_dsigma = coef * invSingularValues * (singulatValues_powa - sigmaProd_ab);
}

template <int dim>
void HyperFoamEnergy<dim>::compute_d2E_div_dsigma2(const Eigen::Matrix<double, dim, 1>& singularValues,
    const MaterialProps& mat,
    Eigen::Matrix<double, dim, dim>& d2E_div_dsigma2) const
{
    if (mat.u == 0.0 && mat.lambda == 0.0) {
        d2E_div_dsigma2.setZero();
        return;
    }

    const double sigmaProd_ab = std::pow(singularValues.prod(), -mat.alpha*mat.beta);

    const auto singularValuesSq = singularValues * singularValues.transpose();

    const double coef = 2.0*mat.u*mat.beta*sigmaProd_ab;
    d2E_div_dsigma2 = singularValuesSq.array().inverse() * coef;

    const double coef2 = 2.0*mat.u/mat.alpha;
    const auto term2 = singularValues.array().pow(-2.0)*sigmaProd_ab;
    const auto term3 = singularValues.array().pow(mat.alpha-2.0)*(mat.alpha-1.0);
    const auto term23 = (term2 + term3)*coef2;

    d2E_div_dsigma2(0, 0) += term23(0);
    d2E_div_dsigma2(1, 1) += term23(1);
    d2E_div_dsigma2(2, 2) += term23(2);
}

template <int dim>
void HyperFoamEnergy<dim>::compute_BLeftCoef(const Eigen::Matrix<double, dim, 1>& singularValues,
    const MaterialProps& mat,
    Eigen::Matrix<double, dim*(dim - 1) / 2, 1>& BLeftCoef) const
{
    if (mat.u == 0.0 && mat.lambda == 0.0) {
        BLeftCoef.setZero();
        return;
    }

    const double sigmaProd = singularValues.prod();
    for (int cI = 0; cI < dim; cI++) {
        int cJ = (cI + 1) % dim;
        // term1 = J^(-alpha*beta) / (Si * Sj)
        double term1 = std::pow(sigmaProd, -mat.alpha*mat.beta) / singularValues[cI] / singularValues[cJ];
        // term2 = (Si^(alpha-1) - Sj^(alpha-1)) / (Si - Sj)
        double term2 = std::pow(singularValues[cI], mat.alpha-1.) - std::pow(singularValues[cJ], mat.alpha-1.);
        const double eps = 1.0e-6;
        if (std::abs(singularValues[cI] - singularValues[cJ]) < eps){
            term2 /= eps;
        } else {
            term2 /= (singularValues[cI] - singularValues[cJ]);
        }
        // Bij = 0.5*(Pi - Pj)/(Si - Sj) = (mu / alpha)*(term1 + term2)
        BLeftCoef[cI] = mat.u/mat.alpha*(term1 + term2);
    }
}
template <int dim>
void HyperFoamEnergy<dim>::compute_dE_div_dF(const Eigen::Matrix<double, dim, dim>& F,
    const AutoFlipSVD<Eigen::Matrix<double, dim, dim>>& svd,
    const MaterialProps& mat,
    Eigen::Matrix<double, dim, dim>& dE_div_dF) const
{
    if (mat.u == 0.0 && mat.lambda == 0.0) {
        dE_div_dF.setZero();
        return;
    }
    Eigen::Matrix<double, dim, 1> singularValues = svd.singularValues();
    const double J = singularValues.prod();
    // P_hat = (2u/alpha)*(Si^alpha - J^(-alpha*beta)) / Si
    const double coef = 2*mat.u / mat.alpha;
    Eigen::Matrix<double, dim, 1> dE_div_dsigma = coef * 
        (singularValues.array().pow(mat.alpha) - std::pow(J, -mat.alpha*mat.beta)) / singularValues.array();
    // P = U * P_hat * V'
    dE_div_dF = svd.matrixU() * dE_div_dsigma.asDiagonal() * svd.matrixV().transpose();
}

template <int dim>
void HyperFoamEnergy<dim>::checkEnergyVal(const Mesh<dim>& data) const // check with isometric case
{
    std::cout << "check energyVal computation..." << std::endl;

    double err = 0.0;
    for (int triI = 0; triI < data.F.rows(); triI++) {
        AutoFlipSVD<Eigen::Matrix<double, dim, dim>> svd(Eigen::Matrix<double, dim, dim>::Identity());

        double energyVal;
        compute_E(svd.singularValues(), data.matProps[triI], energyVal);
        err += data.triArea[triI] * energyVal;
    }

    std::cout << "energyVal computation error = " << err << std::endl;
}

template <int dim>
HyperFoamEnergy<dim>::HyperFoamEnergy(void)
    : Energy<dim>(true)
{
}

template class HyperFoamEnergy<DIM>;

} // namespace IPC
