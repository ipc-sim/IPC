//
//  Energy.hpp
//  IPC
//
//  Created by Minchen Li on 8/31/17.
//

#ifndef Energy_hpp
#define Energy_hpp

#include "Mesh.hpp"

#ifdef USE_CLOSEDFORMSVD2D
#include "ClosedFormSVD2d.hpp"
#else
#include "AutoFlipSVD.hpp"
#endif

#include "LinSysSolver.hpp"

#include <iostream>

namespace IPC {

// a class for energy terms in the objective of an optimization problem
template <int dim>
class Energy {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
    const bool needElemInvSafeGuard;

public:
    Energy(bool p_needElemInvSafeGuard);
    virtual ~Energy(void);

public:
    bool getNeedElemInvSafeGuard(void) const;

public:
    // interfaces to time integration methods
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

    virtual void checkEnergyVal(const Mesh<dim>& data) const;
    virtual void checkGradient(const Mesh<dim>& data) const; // check with finite difference method, according to energyVal
    virtual void checkHessian(const Mesh<dim>& data, bool triplet = false) const; // check with finite difference method, according to gradient

    // common computation framework implementation
    virtual void getEnergyValPerElemBySVD(const Mesh<dim>& data, int redoSVD,
        std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
        std::vector<Eigen::Matrix<double, dim, dim>>& F,
        Eigen::VectorXd& energyValPerElem,
        bool uniformWeight = false) const;
    virtual void computeEnergyValBySVD(const Mesh<dim>& data, int redoSVD,
        std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
        std::vector<Eigen::Matrix<double, dim, dim>>& F,
        double coef,
        double& energyVal) const;
    virtual void computeGradientByPK(const Mesh<dim>& data, bool redoSVD,
        std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
        std::vector<Eigen::Matrix<double, dim, dim>>& F,
        double coef,
        Eigen::VectorXd& gradient,
        bool projectDBC = true) const;
    virtual void computeHessianByPK(const Mesh<dim>& data, bool redoSVD,
        std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
        std::vector<Eigen::Matrix<double, dim, dim>>& F,
        double coef,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver,
        bool projectSPD = true,
        bool projectDBC = true) const;

    // common computation framework implementation -- per element computation
    virtual void computeGradientByPK(const Mesh<dim>& data,
        int elemI, bool redoSVD,
        AutoFlipSVD<Eigen::Matrix<double, dim, dim>>& svd,
        Eigen::Matrix<double, dim, dim>& F,
        double coef,
        Eigen::Matrix<double, dim*(dim + 1), 1>& gradient) const;
    virtual void computeHessianByPK(const Mesh<dim>& data,
        int elemI, bool redoSVD,
        AutoFlipSVD<Eigen::Matrix<double, dim, dim>>& svd,
        Eigen::Matrix<double, dim, dim>& F,
        double coef,
        Eigen::Matrix<double, dim*(dim + 1), dim*(dim + 1)>& hessian,
        Eigen::Matrix<int, 1, dim + 1>& vInd,
        bool projectSPD = true,
        bool projectDBC = true) const;

    // subclass interfaces -- components
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
    virtual void compute_dP_div_dF(const AutoFlipSVD<Eigen::Matrix<double, dim, dim>>& svd,
        double u, double lambda,
        Eigen::Matrix<double, dim * dim, dim * dim>& dP_div_dF,
        double w, bool projectSPD = true) const;

    virtual void filterStepSize(const Mesh<dim>& data,
        const Eigen::VectorXd& searchDir,
        double& stepSize) const;

    virtual void unitTest_dE_div_dsigma(std::ostream& os = std::cout) const;
    virtual void unitTest_d2E_div_dsigma2(std::ostream& os = std::cout) const;
    virtual void unitTest_BLeftCoef(std::ostream& os = std::cout) const;
    virtual void unitTest_dE_div_dF(std::ostream& os = std::cout) const;
    virtual void unitTest_dP_div_dF(std::ostream& os = std::cout) const;
};

} // namespace IPC

#endif /* Energy_hpp */
