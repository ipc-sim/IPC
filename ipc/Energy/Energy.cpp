//
//  Energy.cpp
//  IPC
//
//  Created by Minchen Li on 9/4/17.
//

#include "ipc/Energy/Energy.hpp"
#include "ipc/Utils/get_feasible_steps.hpp"
#include "ipc/Utils/IglUtils.hpp"
#include "ipc/Utils/Timer.hpp"

#include "ipc/LinSysSolver/LinSysSolver.hpp"

#include <igl/avg_edge_length.h>

#ifdef USE_TBB
#include <tbb/tbb.h>
#endif

#include <fstream>
#include <iostream>

extern std::string outputFolderPath;
extern std::ofstream logFile;
extern Timer timer_temp, timer_temp2;

namespace IPC {

template <int dim>
Energy<dim>::Energy(
    bool p_needElemInvSafeGuard)
    : needElemInvSafeGuard(p_needElemInvSafeGuard)
{
}

template <int dim>
Energy<dim>::~Energy(void)
{
}

template <int dim>
bool Energy<dim>::getNeedElemInvSafeGuard(void) const
{
    return needElemInvSafeGuard;
}

template <int dim>
void Energy<dim>::computeEnergyVal(const Mesh<dim>& data, int redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    double& energyVal) const
{
    assert(0 && "please implement this method in the subclass!");
}
template <int dim>
void Energy<dim>::computeGradient(const Mesh<dim>& data, bool redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    Eigen::VectorXd& gradient,
    bool projectDBC) const
{
    assert(0 && "please implement this method in the subclass!");
}
template <int dim>
void Energy<dim>::computeHessian(const Mesh<dim>& data, bool redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver,
    bool projectSPD,
    bool projectDBC) const
{
    assert(0 && "please implement this method in the subclass!");
}

template <int dim>
void Energy<dim>::checkEnergyVal(const Mesh<dim>& data) const
{
    assert(0 && "please implement this method in the subclass!");
}

template <int dim>
void Energy<dim>::checkGradient(const Mesh<dim>& data) const
{
    std::cout << "checking energy gradient computation..." << std::endl;

    std::vector<Eigen::Matrix<double, dim, dim>> F(data.F.rows());
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>> svd(data.F.rows());

    double energyVal0;
    computeEnergyVal(data, true, svd, F, 1.0, energyVal0);
    std::cout << "checkGrad E0 = " << energyVal0 << std::endl;
    const double h = 1.0e-6 * igl::avg_edge_length(data.V, data.F);
    Mesh<dim> perturbed = data;
    Eigen::VectorXd gradient_finiteDiff;
    gradient_finiteDiff.resize(data.V.rows() * dim);
    for (int vI = 0; vI < data.V.rows(); vI++) {
        for (int dimI = 0; dimI < dim; dimI++) {
            perturbed.V = data.V;
            perturbed.V(vI, dimI) += h;
            double energyVal_perturbed;
            computeEnergyVal(perturbed, true, svd, F, 1.0, energyVal_perturbed);
            gradient_finiteDiff[vI * dim + dimI] = (energyVal_perturbed - energyVal0) / h;
        }

        if (((vI + 1) % 100) == 0) {
            std::cout << vI + 1 << "/" << data.V.rows() << " vertices computed" << std::endl;
        }
    }
    for (const auto vI : data.DBCVertexIds) {
        gradient_finiteDiff.segment<dim>(dim * vI).setZero();
    }

    Eigen::VectorXd gradient_symbolic;
    computeGradient(data, true, svd, F, 1.0, gradient_symbolic);

    Eigen::VectorXd difVec = gradient_symbolic - gradient_finiteDiff;
    const double dif_L2 = difVec.norm();
    const double relErr = dif_L2 / gradient_finiteDiff.norm();

    std::cout << "L2 dist = " << dif_L2 << ", relErr = " << relErr << std::endl;

    logFile << "check gradient:" << std::endl;
    logFile << "g_symbolic =\n"
            << gradient_symbolic << std::endl;
    logFile << "g_finiteDiff = \n"
            << gradient_finiteDiff << std::endl;
}

template <int dim>
void Energy<dim>::checkHessian(const Mesh<dim>& data, bool triplet) const
{
    std::cout << "checking energy hessian computation..." << std::endl;

    std::vector<Eigen::Matrix<double, dim, dim>> F(data.F.rows());
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>> svd(data.F.rows());

    Eigen::VectorXd gradient0;
    computeGradientByPK(data, true, svd, F, 1.0, gradient0);
    const double h = 1.0e-6 * igl::avg_edge_length(data.V, data.F);
    Mesh<dim> perturbed = data;
    Eigen::SparseMatrix<double> hessian_finiteDiff;
    hessian_finiteDiff.resize(data.V.rows() * dim, data.V.rows() * dim);
    for (int vI = 0; vI < data.V.rows(); vI++) {
        if (data.DBCVertexIds.find(vI) != data.DBCVertexIds.end()) {
            hessian_finiteDiff.insert(vI * dim, vI * dim) = 1.0;
            hessian_finiteDiff.insert(vI * dim + 1, vI * dim + 1) = 1.0;
            { // Note: it was if constexpr (dim == 3) {
                hessian_finiteDiff.insert(vI * dim + 2, vI * dim + 2) = 1.0;
            }
            continue;
        }

        for (int dimI = 0; dimI < dim; dimI++) {
            perturbed.V = data.V;
            perturbed.V(vI, dimI) += h;
            Eigen::VectorXd gradient_perturbed;
            computeGradientByPK(perturbed, true, svd, F, 1.0, gradient_perturbed);
            Eigen::VectorXd hessian_colI = (gradient_perturbed - gradient0) / h;
            int colI = vI * dim + dimI;
            for (int rowI = 0; rowI < data.V.rows() * dim; rowI++) {
                if (data.DBCVertexIds.find(rowI / dim) == data.DBCVertexIds.end() && hessian_colI[rowI] != 0.0) {
                    hessian_finiteDiff.insert(rowI, colI) = hessian_colI[rowI];
                }
            }
        }

        if (((vI + 1) % 100) == 0) {
            std::cout << vI + 1 << "/" << data.V.rows() << " vertices computed" << std::endl;
        }
    }

    Eigen::SparseMatrix<double> hessian_symbolicPK;

    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver;
    linSysSolver = LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>::create(IPC_DEFAULT_LINSYSSOLVER);
    linSysSolver->set_pattern(data.vNeighbor, data.DBCVertexIds);
    linSysSolver->setZero();
    computeHessianByPK(data, true, svd, F, 1.0, linSysSolver, false);
    linSysSolver->getCoeffMtr(hessian_symbolicPK);

    Eigen::SparseMatrix<double> difMtrPK = hessian_symbolicPK - hessian_finiteDiff;
    const double difPK_L2 = difMtrPK.norm();
    const double relErrPK = difPK_L2 / hessian_finiteDiff.norm();
    std::cout << "PK L2 dist = " << difPK_L2 << ", relErr = " << relErrPK << std::endl;
    IglUtils::writeSparseMatrixToFile(outputFolderPath + "H_symbolicPK", hessian_symbolicPK, true);

    IglUtils::writeSparseMatrixToFile(outputFolderPath + "H_finiteDiff", hessian_finiteDiff, true);
}

template <int dim>
void Energy<dim>::getEnergyValPerElemBySVD(const Mesh<dim>& data, int redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    Eigen::VectorXd& energyValPerElem,
    bool uniformWeight) const
{
    energyValPerElem.resize(data.F.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)data.F.rows(), 1, [&](int triI)
#else
#pragma omp parallel for
    for (int triI = 0; triI < data.F.rows(); triI++)
#endif
        {
            if (redoSVD) {
                const Eigen::Matrix<int, 1, dim + 1>& triVInd = data.F.row(triI);
                Eigen::Matrix<double, dim, dim> Xt;
                Xt.col(0) = (data.V.row(triVInd[1]) - data.V.row(triVInd[0])).transpose();
                Xt.col(1) = (data.V.row(triVInd[2]) - data.V.row(triVInd[0])).transpose();
                { // Note: it was if constexpr (dim == 3) {
                    Xt.col(2) = (data.V.row(triVInd[3]) - data.V.row(triVInd[0])).transpose();
                }
                F[triI] = Xt * data.restTriInv[triI];

                svd[triI].compute(F[triI], Eigen::ComputeFullU | Eigen::ComputeFullV);
                //                Eigen::Matrix2d F = Xt * A;
                //                fprintf(out, "%le %le %le %le\n", F(0, 0), F(0, 1), F(1, 0), F(1, 1));
            }

            compute_E(svd[triI].singularValues(), data.matProps[triI], energyValPerElem[triI]);
            if (!uniformWeight) {
                energyValPerElem[triI] *= data.triArea[triI];
            }
        }
#ifdef USE_TBB
    );
#endif
}

template <int dim>
void Energy<dim>::computeEnergyValBySVD(const Mesh<dim>& data, int redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    double& energyVal) const
{
    Eigen::VectorXd energyValPerElem;
    getEnergyValPerElemBySVD(data, redoSVD, svd, F, energyValPerElem);
    energyVal = coef * energyValPerElem.sum();
}

template <int dim>
void Energy<dim>::computeGradientByPK(const Mesh<dim>& data, bool redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    Eigen::VectorXd& gradient,
    bool projectDBC) const
{
    std::vector<Eigen::Matrix<double, dim*(dim + 1), 1>> gradient_cont(data.F.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)data.F.rows(), 1, [&](int triI)
#else
#pragma omp parallel for
    for (int triI = 0; triI < data.F.rows(); triI++)
#endif
        {
            computeGradientByPK(data, triI, redoSVD,
                svd[triI], F[triI], coef,
                gradient_cont[triI]);
        }
#ifdef USE_TBB
    );
#endif

    gradient.conservativeResize(data.V.rows() * dim);
    gradient.setZero();
#ifdef USE_TBB
    tbb::parallel_for(0, (int)data.V.rows(), 1, [&](int vI)
#else
#pragma omp parallel for
    for (int vI = 0; vI < data.V.rows(); vI++)
#endif
        {
            int _dimVI = vI * dim;
            for (const auto FLocI : data.vFLoc[vI]) {
                gradient.segment<dim>(_dimVI) += gradient_cont[FLocI.first].segment(FLocI.second * dim, dim);
            }
        }
#ifdef USE_TBB
    );
#endif

    if (projectDBC) {
        for (const auto vI : data.DBCVertexIds) {
            gradient.segment<dim>(dim * vI).setZero();
        }
    }
}

template <int dim>
void Energy<dim>::computeHessianByPK(const Mesh<dim>& data, bool redoSVD,
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>>& svd,
    std::vector<Eigen::Matrix<double, dim, dim>>& F,
    double coef,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver,
    bool projectSPD,
    bool projectDBC) const
{
    std::vector<Eigen::Matrix<double, dim*(dim + 1), dim*(dim + 1)>> triHessians(data.F.rows());
    std::vector<Eigen::Matrix<int, 1, dim + 1>> vInds(data.F.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)data.F.rows(), 1, [&](int triI)
#else
#pragma omp parallel for
    for (int triI = 0; triI < data.F.rows(); triI++)
#endif
        {
            computeHessianByPK(data, triI, redoSVD,
                svd[triI], F[triI], coef,
                triHessians[triI], vInds[triI],
                projectSPD, projectDBC);
        }
#ifdef USE_TBB
    );
#endif
#ifdef USE_TBB
    tbb::parallel_for(0, (int)data.V.rows(), 1, [&](int vI)
#else
#pragma omp parallel for
    for (int vI = 0; vI < data.V.rows(); vI++)
#endif
        {
            for (const auto FLocI : data.vFLoc[vI]) {
                IglUtils::addBlockToMatrix<dim>(triHessians[FLocI.first].block(FLocI.second * dim, 0,
                                                    dim, dim * (dim + 1)),
                    vInds[FLocI.first], FLocI.second, linSysSolver);
            }
        }
#ifdef USE_TBB
    );
#endif
}

template <int dim>
void Energy<dim>::computeGradientByPK(const Mesh<dim>& data,
    int elemI, bool redoSVD,
    AutoFlipSVD<Eigen::Matrix<double, dim, dim>>& svd,
    Eigen::Matrix<double, dim, dim>& F,
    double coef,
    Eigen::Matrix<double, dim*(dim + 1), 1>& gradient) const
{
    const Eigen::Matrix<double, dim, dim>& A = data.restTriInv[elemI];

    if (redoSVD) {
        const Eigen::Matrix<int, 1, dim + 1>& triVInd = data.F.row(elemI);

        Eigen::Matrix<double, dim, dim> Xt;
        Xt.col(0) = (data.V.row(triVInd[1]) - data.V.row(triVInd[0])).transpose();
        Xt.col(1) = (data.V.row(triVInd[2]) - data.V.row(triVInd[0])).transpose();
        { // Note: it was if constexpr (dim == 3) {
            Xt.col(2) = (data.V.row(triVInd[3]) - data.V.row(triVInd[0])).transpose();
        }

        F = Xt * A;

        svd.compute(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    }

    Eigen::Matrix<double, dim, dim> P;
    compute_dE_div_dF(F, svd, data.matProps[elemI], P);

    const double w = coef * data.triArea[elemI];
    P *= w;

    IglUtils::dF_div_dx_mult(P, A, gradient);
}
template <int dim>
void Energy<dim>::computeHessianByPK(const Mesh<dim>& data,
    int elemI, bool redoSVD,
    AutoFlipSVD<Eigen::Matrix<double, dim, dim>>& svd,
    Eigen::Matrix<double, dim, dim>& F,
    double coef,
    Eigen::Matrix<double, dim*(dim + 1), dim*(dim + 1)>& hessian,
    Eigen::Matrix<int, 1, dim + 1>& vInd,
    bool projectSPD,
    bool projectDBC) const
{
    const Eigen::Matrix<int, 1, dim + 1>& triVInd = data.F.row(elemI);

    const Eigen::Matrix<double, dim, dim>& A = data.restTriInv[elemI];

    if (redoSVD) {
        Eigen::Matrix<double, dim, dim> Xt;
        Xt.col(0) = (data.V.row(triVInd[1]) - data.V.row(triVInd[0])).transpose();
        Xt.col(1) = (data.V.row(triVInd[2]) - data.V.row(triVInd[0])).transpose();
        { // Note: it was if constexpr (dim == 3) {
            Xt.col(2) = (data.V.row(triVInd[3]) - data.V.row(triVInd[0])).transpose();
        }
        F = Xt * A;
        svd.compute(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    }

    Eigen::Matrix<double, dim * dim, dim * dim> wdP_div_dF;
    const double w = coef * data.triArea[elemI];
    compute_dP_div_dF(svd, data.matProps[elemI], wdP_div_dF, w, projectSPD);

    Eigen::Matrix<double, dim*(dim + 1), dim * dim> wdP_div_dx;
    IglUtils::dF_div_dx_mult<dim * dim>(wdP_div_dF.transpose(), A, wdP_div_dx, false);
    IglUtils::dF_div_dx_mult<dim*(dim + 1)>(wdP_div_dx.transpose(), A, hessian, true);

    vInd[0] = data.isProjectDBCVertex(triVInd[0], projectDBC) ? (-triVInd[0] - 1) : triVInd[0];
    vInd[1] = data.isProjectDBCVertex(triVInd[1], projectDBC) ? (-triVInd[1] - 1) : triVInd[1];
    vInd[2] = data.isProjectDBCVertex(triVInd[2], projectDBC) ? (-triVInd[2] - 1) : triVInd[2];
    { // Note: it was if constexpr (dim == 3) {
        vInd[3] = data.isProjectDBCVertex(triVInd[3], projectDBC) ? (-triVInd[3] - 1) : triVInd[3];
    }
} // namespace IPC

template <int dim>
void Energy<dim>::compute_E(const Eigen::Matrix<double, dim, 1>& singularValues,
    const MaterialProps& mat,
    double& E) const
{
    assert(0 && "please implement this method in the subclass!");
}
template <int dim>
void Energy<dim>::compute_dE_div_dsigma(const Eigen::Matrix<double, dim, 1>& singularValues,
    const MaterialProps& mat,
    Eigen::Matrix<double, dim, 1>& dE_div_dsigma) const
{
    assert(0 && "please implement this method in the subclass!");
}
template <int dim>
void Energy<dim>::compute_d2E_div_dsigma2(const Eigen::Matrix<double, dim, 1>& singularValues,
    const MaterialProps& mat,
    Eigen::Matrix<double, dim, dim>& d2E_div_dsigma2) const
{
    assert(0 && "please implement this method in the subclass!");
}
template <int dim>
void Energy<dim>::compute_BLeftCoef(const Eigen::Matrix<double, dim, 1>& singularValues,
    const MaterialProps& mat,
    Eigen::Matrix<double, dim*(dim - 1) / 2, 1>& BLeftCoef) const
{
    assert(0 && "please implement this method in the subclass!");
}
template <int dim>
void Energy<dim>::compute_dE_div_dF(const Eigen::Matrix<double, dim, dim>& F,
    const AutoFlipSVD<Eigen::Matrix<double, dim, dim>>& svd,
    const MaterialProps& mat,
    Eigen::Matrix<double, dim, dim>& dE_div_dF) const
{
    assert(0 && "please implement this method in the subclass!");
}

template <int dim>
void Energy<dim>::compute_dP_div_dF(const AutoFlipSVD<Eigen::Matrix<double, dim, dim>>& svd,
    const MaterialProps& mat,
    Eigen::Matrix<double, dim * dim, dim * dim>& dP_div_dF,
    double w, bool projectSPD) const
{
    // compute A
    const Eigen::Matrix<double, dim, 1>& sigma = svd.singularValues();
    Eigen::Matrix<double, dim, 1> dE_div_dsigma;
    compute_dE_div_dsigma(sigma, mat, dE_div_dsigma);
    Eigen::Matrix<double, dim, dim> d2E_div_dsigma2;
    compute_d2E_div_dsigma2(sigma, mat, d2E_div_dsigma2);
    if (projectSPD) {
#if (DIM == 2)
        IglUtils::makePD2d(d2E_div_dsigma2);
#else
        IglUtils::makePD(d2E_div_dsigma2); //TODO: use implicit QR to accelerate
#endif
    }

    // compute B
    const int Cdim2 = dim * (dim - 1) / 2;
    Eigen::Matrix<double, Cdim2, 1> BLeftCoef;
    compute_BLeftCoef(sigma, mat, BLeftCoef);
    Eigen::Matrix2d B[Cdim2];
    for (int cI = 0; cI < Cdim2; cI++) {
        int cI_post = (cI + 1) % dim;

        double rightCoef = dE_div_dsigma[cI] + dE_div_dsigma[cI_post];
        double sum_sigma = sigma[cI] + sigma[cI_post];
        const double eps = 1.0e-6;
        if (sum_sigma < eps) {
            rightCoef /= 2.0 * eps;
        }
        else {
            rightCoef /= 2.0 * sum_sigma;
        }

        const double& leftCoef = BLeftCoef[cI];
        B[cI](0, 0) = B[cI](1, 1) = leftCoef + rightCoef;
        B[cI](0, 1) = B[cI](1, 0) = leftCoef - rightCoef;
        if (projectSPD) {
            IglUtils::makePD2d(B[cI]);
        }
    }

    // compute M using A(d2E_div_dsigma2) and B
    Eigen::Matrix<double, dim * dim, dim * dim> M;
    M.setZero();
    if constexpr (dim == 2) {
        M(0, 0) = w * d2E_div_dsigma2(0, 0);
        M(0, 3) = w * d2E_div_dsigma2(0, 1);
        M.block(1, 1, 2, 2) = w * B[0];
        M(3, 0) = w * d2E_div_dsigma2(1, 0);
        M(3, 3) = w * d2E_div_dsigma2(1, 1);
    }
    else {
        // A
        M(0, 0) = w * d2E_div_dsigma2(0, 0);
        M(0, 4) = w * d2E_div_dsigma2(0, 1);
        M(0, 8) = w * d2E_div_dsigma2(0, 2);
        M(4, 0) = w * d2E_div_dsigma2(1, 0);
        M(4, 4) = w * d2E_div_dsigma2(1, 1);
        M(4, 8) = w * d2E_div_dsigma2(1, 2);
        M(8, 0) = w * d2E_div_dsigma2(2, 0);
        M(8, 4) = w * d2E_div_dsigma2(2, 1);
        M(8, 8) = w * d2E_div_dsigma2(2, 2);
        // B01
        M(1, 1) = w * B[0](0, 0);
        M(1, 3) = w * B[0](0, 1);
        M(3, 1) = w * B[0](1, 0);
        M(3, 3) = w * B[0](1, 1);
        // B12
        M(5, 5) = w * B[1](0, 0);
        M(5, 7) = w * B[1](0, 1);
        M(7, 5) = w * B[1](1, 0);
        M(7, 7) = w * B[1](1, 1);
        // B20
        M(2, 2) = w * B[2](1, 1);
        M(2, 6) = w * B[2](1, 0);
        M(6, 2) = w * B[2](0, 1);
        M(6, 6) = w * B[2](0, 0);
    }

    // compute dP_div_dF
    Eigen::Matrix<double, dim * dim, dim* dim>& wdP_div_dF = dP_div_dF;
    const Eigen::Matrix<double, dim, dim>& U = svd.matrixU();
    const Eigen::Matrix<double, dim, dim>& V = svd.matrixV();
    for (int i = 0; i < dim; i++) {
        int _dim_i = i * dim;
        for (int j = 0; j < dim; j++) {
            int ij = _dim_i + j;
            for (int r = 0; r < dim; r++) {
                int _dim_r = r * dim;
                for (int s = 0; s < dim; s++) {
                    int rs = _dim_r + s;
                    if (ij > rs) {
                        // bottom left, same as upper right
                        continue;
                    }

                    if constexpr (dim == 2) {
                        wdP_div_dF(ij, rs) = M(0, 0) * U(i, 0) * V(j, 0) * U(r, 0) * V(s, 0) + M(0, 3) * U(i, 0) * V(j, 0) * U(r, 1) * V(s, 1) + M(1, 1) * U(i, 0) * V(j, 1) * U(r, 0) * V(s, 1) + M(1, 2) * U(i, 0) * V(j, 1) * U(r, 1) * V(s, 0) + M(2, 1) * U(i, 1) * V(j, 0) * U(r, 0) * V(s, 1) + M(2, 2) * U(i, 1) * V(j, 0) * U(r, 1) * V(s, 0) + M(3, 0) * U(i, 1) * V(j, 1) * U(r, 0) * V(s, 0) + M(3, 3) * U(i, 1) * V(j, 1) * U(r, 1) * V(s, 1);
                    }
                    else {
                        wdP_div_dF(ij, rs) = M(0, 0) * U(i, 0) * V(j, 0) * U(r, 0) * V(s, 0) + M(0, 4) * U(i, 0) * V(j, 0) * U(r, 1) * V(s, 1) + M(0, 8) * U(i, 0) * V(j, 0) * U(r, 2) * V(s, 2) + M(4, 0) * U(i, 1) * V(j, 1) * U(r, 0) * V(s, 0) + M(4, 4) * U(i, 1) * V(j, 1) * U(r, 1) * V(s, 1) + M(4, 8) * U(i, 1) * V(j, 1) * U(r, 2) * V(s, 2) + M(8, 0) * U(i, 2) * V(j, 2) * U(r, 0) * V(s, 0) + M(8, 4) * U(i, 2) * V(j, 2) * U(r, 1) * V(s, 1) + M(8, 8) * U(i, 2) * V(j, 2) * U(r, 2) * V(s, 2) + M(1, 1) * U(i, 0) * V(j, 1) * U(r, 0) * V(s, 1) + M(1, 3) * U(i, 0) * V(j, 1) * U(r, 1) * V(s, 0) + M(3, 1) * U(i, 1) * V(j, 0) * U(r, 0) * V(s, 1) + M(3, 3) * U(i, 1) * V(j, 0) * U(r, 1) * V(s, 0) + M(5, 5) * U(i, 1) * V(j, 2) * U(r, 1) * V(s, 2) + M(5, 7) * U(i, 1) * V(j, 2) * U(r, 2) * V(s, 1) + M(7, 5) * U(i, 2) * V(j, 1) * U(r, 1) * V(s, 2) + M(7, 7) * U(i, 2) * V(j, 1) * U(r, 2) * V(s, 1) + M(2, 2) * U(i, 0) * V(j, 2) * U(r, 0) * V(s, 2) + M(2, 6) * U(i, 0) * V(j, 2) * U(r, 2) * V(s, 0) + M(6, 2) * U(i, 2) * V(j, 0) * U(r, 0) * V(s, 2) + M(6, 6) * U(i, 2) * V(j, 0) * U(r, 2) * V(s, 0);
                    }

                    if (ij < rs) {
                        wdP_div_dF(rs, ij) = wdP_div_dF(ij, rs);
                    }
                }
            }
        }
    }
}

template <int dim>
void Energy<dim>::filterStepSize(const Mesh<dim>& data, const Eigen::VectorXd& searchDir, double& stepSize) const
{
    if (needElemInvSafeGuard) {
        Eigen::VectorXd output(data.F.rows());
        if constexpr (dim == 2) {
            computeInjectiveStepSize_2d(data.F, data.V, searchDir, 1.0e-6, 0.1, output.data());
        }
        else {
            computeInjectiveStepSize_3d(data.F, data.V, searchDir, 1.0e-6, 0.2, output.data());
        }

        double tentativeStepSize = output.minCoeff();
        if ((tentativeStepSize > 0.0) && (tentativeStepSize < stepSize)) {
            stepSize = tentativeStepSize;
        }
    }
}

template <int dim>
void Energy<dim>::unitTest_dE_div_dsigma(std::ostream& os) const
{
    std::vector<Eigen::Matrix<double, dim, 1>> testSigma;

    testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Ones());
    if (needElemInvSafeGuard) {
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
    }
    else {
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Zero());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
    }

    double YM = 100, PR = 0.4;
    MaterialProps mat = {YM, PR};

    const double h = 1.0e-6;
    int testI = 0;
    for (const auto& testSigmaI : testSigma) {
        os << "--- unitTest_dE_div_dsigma " << testI << " ---\n"
           << "sigma =\n"
           << testSigmaI << std::endl;

        double E0;
        compute_E(testSigmaI, mat, E0);

        Eigen::Matrix<double, dim, 1> dE_div_dsigma_FD;
        for (int dimI = 0; dimI < dim; dimI++) {
            Eigen::Matrix<double, dim, 1> sigma_perterb = testSigmaI;
            sigma_perterb[dimI] += h;
            double E;
            compute_E(sigma_perterb, mat, E);
            dE_div_dsigma_FD[dimI] = (E - E0) / h;
        }
        os << "dE_div_dsigma_FD =\n"
           << dE_div_dsigma_FD << std::endl;

        Eigen::Matrix<double, dim, 1> dE_div_dsigma_S;
        compute_dE_div_dsigma(testSigmaI, mat, dE_div_dsigma_S);
        os << "dE_div_dsigma_S =\n"
           << dE_div_dsigma_S << std::endl;

        double err = (dE_div_dsigma_FD - dE_div_dsigma_S).norm();
        os << "err = " << err << " (" << err / dE_div_dsigma_FD.norm() * 100 << "%)" << std::endl;

        ++testI;
    }
}

template <int dim>
void Energy<dim>::unitTest_d2E_div_dsigma2(std::ostream& os) const
{
    std::vector<Eigen::Matrix<double, dim, 1>> testSigma;

    testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Ones());
    if (needElemInvSafeGuard) {
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
    }
    else {
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Zero());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
    }

    double YM = 100, PR = 0.4;
    MaterialProps mat = {YM, PR};

    const double h = 1.0e-6;
    int testI = 0;
    for (const auto& testSigmaI : testSigma) {
        os << "--- unitTest_d2E_div_dsigma2 " << testI << " ---\n"
           << "sigma =\n"
           << testSigmaI << std::endl;

        Eigen::Matrix<double, dim, 1> dE_div_dsigma0;
        compute_dE_div_dsigma(testSigmaI, mat, dE_div_dsigma0);

        Eigen::Matrix<double, dim, dim> d2E_div_dsigma2_FD;
        for (int dimI = 0; dimI < dim; dimI++) {
            Eigen::Matrix<double, dim, 1> sigma_perterb = testSigmaI;
            sigma_perterb[dimI] += h;
            Eigen::Matrix<double, dim, 1> dE_div_dsigma;
            compute_dE_div_dsigma(sigma_perterb, mat, dE_div_dsigma);
            d2E_div_dsigma2_FD.row(dimI) = ((dE_div_dsigma - dE_div_dsigma0) / h).transpose();
        }
        os << "d2E_div_dsigma2_FD =\n"
           << d2E_div_dsigma2_FD << std::endl;

        Eigen::Matrix<double, dim, dim> d2E_div_dsigma2_S;
        compute_d2E_div_dsigma2(testSigmaI, mat, d2E_div_dsigma2_S);
        os << "d2E_div_dsigma2_S =\n"
           << d2E_div_dsigma2_S << std::endl;

        double err = (d2E_div_dsigma2_FD - d2E_div_dsigma2_S).norm();
        os << "err = " << err << " (" << err / d2E_div_dsigma2_FD.norm() * 100 << "%)" << std::endl;

        ++testI;
    }
}

template <int dim>
void Energy<dim>::unitTest_BLeftCoef(std::ostream& os) const
{
    std::vector<Eigen::Matrix<double, dim, 1>> testSigma;

    if (needElemInvSafeGuard) {
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() + Eigen::Matrix<double, dim, 1>::Ones());
    }
    else {
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Zero());
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
        testSigma.emplace_back(Eigen::Matrix<double, dim, 1>::Random() * 2);
    }

    double YM = 100, PR = 0.4;
    MaterialProps mat = {YM, PR};

    int testI = 0;
    for (const auto& testSigmaI : testSigma) {
        os << "--- unitTest_BLeftCoef " << testI << " ---\n"
           << "sigma =\n"
           << testSigmaI << std::endl;

        Eigen::Matrix<double, dim*(dim - 1) / 2, 1> BLeftCoef_div;
        Eigen::Matrix<double, dim, 1> dE_div_dsigma;
        compute_dE_div_dsigma(testSigmaI, mat, dE_div_dsigma);
        if constexpr (dim == 2) {
            BLeftCoef_div[0] = (dE_div_dsigma[0] - dE_div_dsigma[1]) / (testSigmaI[0] - testSigmaI[1]) / 2.0;
        }
        else {
            BLeftCoef_div[0] = (dE_div_dsigma[0] - dE_div_dsigma[1]) / (testSigmaI[0] - testSigmaI[1]) / 2.0;
            BLeftCoef_div[1] = (dE_div_dsigma[1] - dE_div_dsigma[2]) / (testSigmaI[1] - testSigmaI[2]) / 2.0;
            BLeftCoef_div[2] = (dE_div_dsigma[2] - dE_div_dsigma[0]) / (testSigmaI[2] - testSigmaI[0]) / 2.0;
        }
        os << "BLeftCoef_div =\n"
           << BLeftCoef_div << std::endl;

        Eigen::Matrix<double, dim*(dim - 1) / 2, 1> BLeftCoef_S;
        compute_BLeftCoef(testSigmaI, mat, BLeftCoef_S);
        os << "BLeftCoef_S =\n"
           << BLeftCoef_S << std::endl;

        double err = (BLeftCoef_div - BLeftCoef_S).norm();
        os << "err = " << err << " (" << err / BLeftCoef_div.norm() * 100 << "%)" << std::endl;

        ++testI;
    }
}

template <int dim>
void Energy<dim>::unitTest_dE_div_dF(std::ostream& os) const
{
    const double h = 1.0e-6;
    std::vector<Eigen::Matrix<double, dim, dim>> testF;

    testF.emplace_back(Eigen::Matrix<double, dim, dim>::Identity());
    if (needElemInvSafeGuard) {
        for (int testI = 0; testI < 6; testI++) {
            testF.emplace_back(Eigen::Matrix<double, dim, dim>::Random());
            IglUtils::flipDet_SVD(testF.back());
            testF.back() += 1.0e-2 * Eigen::Matrix<double, dim, dim>::Identity();
        }
    }
    else {
        testF.emplace_back(Eigen::Matrix<double, dim, dim>::Zero());
        for (int testI = 0; testI < 6; testI++) {
            testF.emplace_back(Eigen::Matrix<double, dim, dim>::Random());
        }
    }

    double YM = 100, PR = 0.4;
    MaterialProps mat = {YM, PR};

    int testI = 0;
    for (const auto& testFI : testF) {
        os << "--- unitTest_dE_div_dF " << testI << " ---\n"
           << "F =\n"
           << testFI << std::endl;

        AutoFlipSVD<Eigen::Matrix<double, dim, dim>> svd0(testFI, Eigen::ComputeFullU | Eigen::ComputeFullV);
        double E0;
        compute_E(svd0.singularValues(), mat, E0);

        Eigen::Matrix<double, dim, dim> P_FD;
        for (int dimI = 0; dimI < dim; dimI++) {
            for (int dimJ = 0; dimJ < dim; dimJ++) {
                Eigen::Matrix<double, dim, dim> F_perterb = testFI;
                F_perterb(dimI, dimJ) += h;

                AutoFlipSVD<Eigen::Matrix<double, dim, dim>> svd(F_perterb, Eigen::ComputeFullU | Eigen::ComputeFullV);
                double E;
                compute_E(svd.singularValues(), mat, E);

                P_FD(dimI, dimJ) = (E - E0) / h;
            }
        }
        os << "P_FD =\n"
           << P_FD << std::endl;

        Eigen::Matrix<double, dim, dim> P_S;
        compute_dE_div_dF(testFI, svd0, mat, P_S);
        os << "P_S =\n"
           << P_S << std::endl;

        double err = (P_FD - P_S).norm();
        os << "err = " << err << " (" << err / P_FD.norm() * 100 << "%)" << std::endl;

        ++testI;
    }
}

template <int dim>
void Energy<dim>::unitTest_dP_div_dF(std::ostream& os) const
{
    const double h = 1.0e-6;
    std::vector<Eigen::Matrix<double, dim, dim>> testF;

    testF.emplace_back(Eigen::Matrix<double, dim, dim>::Identity());
    if (needElemInvSafeGuard) {
        for (int testI = 0; testI < 6; testI++) {
            testF.emplace_back(Eigen::Matrix<double, dim, dim>::Random());
            IglUtils::flipDet_SVD(testF.back());
            testF.back() += 1.0e-2 * Eigen::Matrix<double, dim, dim>::Identity();
        }
    }
    else {
        testF.emplace_back(Eigen::Matrix<double, dim, dim>::Zero());
        for (int testI = 0; testI < 6; testI++) {
            testF.emplace_back(Eigen::Matrix<double, dim, dim>::Random());
        }
    }

    double YM = 100, PR = 0.4;
    MaterialProps mat = {YM, PR};

    int testI = 0;
    for (const auto& testFI : testF) {
        os << "--- unitTest_dP_div_dF " << testI << " ---\n"
           << "F =\n"
           << testFI << std::endl;

        AutoFlipSVD<Eigen::Matrix<double, dim, dim>> svd(testFI, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix<double, dim, dim> P0;
        compute_dE_div_dF(testFI, svd, mat, P0);

        Eigen::Matrix<double, dim * dim, dim * dim> dP_div_dF_FD;
        for (int dimI = 0; dimI < dim; dimI++) {
            for (int dimJ = 0; dimJ < dim; dimJ++) {
                Eigen::Matrix<double, dim, dim> F_perterb = testFI;
                F_perterb(dimI, dimJ) += h;
                Eigen::Matrix<double, dim, dim> P;
                AutoFlipSVD<Eigen::Matrix<double, dim, dim>> svd(F_perterb, Eigen::ComputeFullU | Eigen::ComputeFullV);
                compute_dE_div_dF(F_perterb, svd, mat, P);

                Eigen::Matrix<double, dim, dim> FD = (P - P0) / h;
                dP_div_dF_FD.block(0, dimI * dim + dimJ, dim, 1) = FD.row(0).transpose();
                dP_div_dF_FD.block(dim, dimI * dim + dimJ, dim, 1) = FD.row(1).transpose();
                { // Note: it was if constexpr (dim == 3) {
                    dP_div_dF_FD.block(dim * 2, dimI * dim + dimJ, dim, 1) = FD.row(2).transpose();
                }
            }
        }
        os << "dP_div_dF_FD =\n"
           << dP_div_dF_FD << std::endl;

        Eigen::Matrix<double, dim * dim, dim * dim> dP_div_dF_S;
        svd.compute(testFI, Eigen::ComputeFullU | Eigen::ComputeFullV);
        compute_dP_div_dF(svd, mat, dP_div_dF_S, 1.0, false);
        os << "dP_div_dF_S =\n"
           << dP_div_dF_S << std::endl;

        double err = (dP_div_dF_FD - dP_div_dF_S).norm();
        os << "err = " << err << " (" << err / dP_div_dF_FD.norm() * 100 << "%)" << std::endl;

        ++testI;
    }
}

template class Energy<DIM>;

} // namespace IPC
