//
//  FixedCoRotEnergy.hpp
//  IPC
//
//  Created by Minchen Li on 6/20/18.
//

#ifndef FixedCoRotEnergy_hpp
#define FixedCoRotEnergy_hpp

#include "Energy.hpp"

namespace IPC {

template <int dim>
class FixedCoRotEnergy : public Energy<dim> {
    typedef Energy<dim> Base;

public:
    virtual void computeEnergyVal(const Mesh<dim>& data, int redoSVD,
        std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
        std::vector<Eigen::Matrix<double, dim, dim>>& F,
        double coef,
        double& energyVal) const;
    virtual void computeGradient(const Mesh<dim>& data, bool redoSVD,
        std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
        std::vector<Eigen::Matrix<double, dim, dim>>& F,
        double coef,
        Eigen::VectorXd& gradient,
        bool projectDBC = true) const;
    virtual void computeHessian(const Mesh<dim>& data, bool redoSVD,
        std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
        std::vector<Eigen::Matrix<double, dim, dim>>& F,
        double coef,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver,
        bool projectSPD = true,
        bool projectDBC = true) const;

    virtual void getEnergyValPerElem(const Mesh<dim>& data, Eigen::VectorXd& energyValPerElem, bool uniformWeight = false) const;

    virtual void compute_E(const Eigen::Matrix<double, dim, 1>& singularValues,
        double u, double lambda,
        double& E) const;
    virtual void compute_dE_div_dsigma(const Eigen::Matrix<double, dim, 1>& singularValues,
        double u, double lambda,
        Eigen::Matrix<double, dim, 1>& dE_div_dsigma) const;
    virtual void compute_d2E_div_dsigma2(const Eigen::Matrix<double, dim, 1>& singularValues,
        double u, double lambda,
        Eigen::Matrix<double, dim, dim>& d2E_div_dsigma2) const;
    virtual void compute_BLeftCoef(const Eigen::Matrix<double, dim, 1>& singularValues,
        double u, double lambda,
        Eigen::Matrix<double, dim*(dim - 1) / 2, 1>& BLeftCoef) const;
    virtual void compute_dE_div_dF(const Eigen::Matrix<double, dim, dim>& F,
        const AutoFlipSVD<Eigen::Matrix<double, dim, dim>>& svd,
        double u, double lambda,
        Eigen::Matrix<double, dim, dim>& dE_div_dF) const;
    virtual void checkEnergyVal(const Mesh<dim>& data) const; // check with isometric case

public:
    FixedCoRotEnergy(void);
};

} // namespace IPC

#endif /* FixedCoRotEnergy_hpp */
