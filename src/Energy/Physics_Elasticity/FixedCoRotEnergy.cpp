//
//  FixedCoRotEnergy.cpp
//  IPC
//
//  Created by Minchen Li on 6/20/18.
//

#include "FixedCoRotEnergy.hpp"
#include "IglUtils.hpp"
#include "Timer.hpp"

#ifdef USE_TBB
#include <tbb/parallel_for.h>
#endif

extern Timer timer_temp;

namespace IPC {

template <int dim>
void FixedCoRotEnergy<dim>::computeEnergyVal(const Mesh<dim>& data, int redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    double& energyVal) const
{
    Base::computeEnergyValBySVD(data, redoSVD, svd, F, coef, energyVal);
}
template <int dim>
void FixedCoRotEnergy<dim>::computeGradient(const Mesh<dim>& data, bool redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    Eigen::VectorXd& gradient,
    bool projectDBC) const
{
    Base::computeGradientByPK(data, redoSVD, svd, F, coef, gradient, projectDBC);
}
template <int dim>
void FixedCoRotEnergy<dim>::computeHessian(const Mesh<dim>& data, bool redoSVD,
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
void FixedCoRotEnergy<dim>::getEnergyValPerElem(const Mesh<dim>& data,
    Eigen::VectorXd& energyValPerElem,
    bool uniformWeight) const
{
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>> svd(data.F.rows());
    std::vector<Eigen::Matrix<double, dim, dim>> F(data.F.rows());
    Energy<dim>::getEnergyValPerElemBySVD(data, true, svd, F, energyValPerElem, uniformWeight);
}

template <int dim>
void FixedCoRotEnergy<dim>::compute_E(const Eigen::Matrix<double, dim, 1>& singularValues,
    double u, double lambda,
    double& E) const
{
    const double sigmam12Sum = (singularValues - Eigen::Matrix<double, dim, 1>::Ones()).squaredNorm();
    const double sigmaProdm1 = singularValues.prod() - 1.0;

    E = u * sigmam12Sum + lambda / 2.0 * sigmaProdm1 * sigmaProdm1;
}
template <int dim>
void FixedCoRotEnergy<dim>::compute_dE_div_dsigma(const Eigen::Matrix<double, dim, 1>& singularValues,
    double u, double lambda,
    Eigen::Matrix<double, dim, 1>& dE_div_dsigma) const
{
    const double sigmaProdm1lambda = lambda * (singularValues.prod() - 1.0);
    Eigen::Matrix<double, dim, 1> sigmaProd_noI;
    if constexpr (dim == 2) {
        sigmaProd_noI[0] = singularValues[1];
        sigmaProd_noI[1] = singularValues[0];
    }
    else {
        sigmaProd_noI[0] = singularValues[1] * singularValues[2];
        sigmaProd_noI[1] = singularValues[2] * singularValues[0];
        sigmaProd_noI[2] = singularValues[0] * singularValues[1];
    }

    double _2u = u * 2;
    dE_div_dsigma[0] = (_2u * (singularValues[0] - 1.0) + sigmaProd_noI[0] * sigmaProdm1lambda);
    dE_div_dsigma[1] = (_2u * (singularValues[1] - 1.0) + sigmaProd_noI[1] * sigmaProdm1lambda);
    if constexpr (dim == 3) {
        dE_div_dsigma[2] = (_2u * (singularValues[2] - 1.0) + sigmaProd_noI[2] * sigmaProdm1lambda);
    }
}
template <int dim>
void FixedCoRotEnergy<dim>::compute_d2E_div_dsigma2(const Eigen::Matrix<double, dim, 1>& singularValues,
    double u, double lambda,
    Eigen::Matrix<double, dim, dim>& d2E_div_dsigma2) const
{
    const double sigmaProd = singularValues.prod();
    Eigen::Matrix<double, dim, 1> sigmaProd_noI;
    if constexpr (dim == 2) {
        sigmaProd_noI[0] = singularValues[1];
        sigmaProd_noI[1] = singularValues[0];
    }
    else {
        sigmaProd_noI[0] = singularValues[1] * singularValues[2];
        sigmaProd_noI[1] = singularValues[2] * singularValues[0];
        sigmaProd_noI[2] = singularValues[0] * singularValues[1];
    }

    double _2u = u * 2;
    d2E_div_dsigma2(0, 0) = _2u + lambda * sigmaProd_noI[0] * sigmaProd_noI[0];
    d2E_div_dsigma2(1, 1) = _2u + lambda * sigmaProd_noI[1] * sigmaProd_noI[1];
    if constexpr (dim == 3) {
        d2E_div_dsigma2(2, 2) = _2u + lambda * sigmaProd_noI[2] * sigmaProd_noI[2];
    }

    if constexpr (dim == 2) {
        d2E_div_dsigma2(0, 1) = d2E_div_dsigma2(1, 0) = lambda * ((sigmaProd - 1.0) + sigmaProd_noI[0] * sigmaProd_noI[1]);
    }
    else {
        d2E_div_dsigma2(0, 1) = d2E_div_dsigma2(1, 0) = lambda * (singularValues[2] * (sigmaProd - 1.0) + sigmaProd_noI[0] * sigmaProd_noI[1]);
        d2E_div_dsigma2(0, 2) = d2E_div_dsigma2(2, 0) = lambda * (singularValues[1] * (sigmaProd - 1.0) + sigmaProd_noI[0] * sigmaProd_noI[2]);
        d2E_div_dsigma2(2, 1) = d2E_div_dsigma2(1, 2) = lambda * (singularValues[0] * (sigmaProd - 1.0) + sigmaProd_noI[2] * sigmaProd_noI[1]);
    }
}
template <int dim>
void FixedCoRotEnergy<dim>::compute_BLeftCoef(const Eigen::Matrix<double, dim, 1>& singularValues,
    double u, double lambda,
    Eigen::Matrix<double, dim*(dim - 1) / 2, 1>& BLeftCoef) const
{
    const double sigmaProd = singularValues.prod();
    const double halfLambda = lambda / 2.0;
    if constexpr (dim == 2) {
        BLeftCoef[0] = u - halfLambda * (sigmaProd - 1);
    }
    else {
        BLeftCoef[0] = u - halfLambda * singularValues[2] * (sigmaProd - 1);
        BLeftCoef[1] = u - halfLambda * singularValues[0] * (sigmaProd - 1);
        BLeftCoef[2] = u - halfLambda * singularValues[1] * (sigmaProd - 1);
    }
}
template <int dim>
void FixedCoRotEnergy<dim>::compute_dE_div_dF(const Eigen::Matrix<double, dim, dim>& F,
    const AutoFlipSVD<Eigen::Matrix<double, dim, dim>>& svd,
    double u, double lambda,
    Eigen::Matrix<double, dim, dim>& dE_div_dF) const
{
    Eigen::Matrix<double, dim, dim> JFInvT;
    IglUtils::computeCofactorMtr(F, JFInvT);
    dE_div_dF = (u * 2 * (F - svd.matrixU() * svd.matrixV().transpose()) + lambda * (svd.singularValues().prod() - 1) * JFInvT);
}

template <int dim>
void FixedCoRotEnergy<dim>::checkEnergyVal(const Mesh<dim>& data) const // check with isometric case
{
    std::cout << "check energyVal computation..." << std::endl;

    double err = 0.0;
    for (int triI = 0; triI < data.F.rows(); triI++) {
        AutoFlipSVD<Eigen::Matrix<double, dim, dim>> svd(Eigen::Matrix<double, dim, dim>::Identity());

        double energyVal;
        compute_E(svd.singularValues(), data.u[triI], data.lambda[triI], energyVal);
        err += data.triArea[triI] * energyVal;
    }

    std::cout << "energyVal computation error = " << err << std::endl;
}

template <int dim>
FixedCoRotEnergy<dim>::FixedCoRotEnergy(void)
    : Energy<dim>(false)
{
}

template class FixedCoRotEnergy<DIM>;

} // namespace IPC
