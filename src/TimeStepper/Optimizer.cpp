//
//  Optimizer.cpp
//  IPC
//
//  Created by Minchen Li on 8/31/17.
//

#include "Optimizer.hpp"

#include "SelfCollisionHandler.hpp"

#ifdef LINSYSSOLVER_USE_CHOLMOD
#include "CHOLMODSolver.hpp"
#elif defined(LINSYSSOLVER_USE_AMGCL)
#include "AMGCLSolver.hpp"
#else
#include "EigenLibSolver.hpp"
#endif

#include "IglUtils.hpp"
#include "BarrierFunctions.hpp"
#include "Timer.hpp"

#include <igl/avg_edge_length.h>
#include <igl/face_areas.h>
#include <igl/edge_lengths.h>
#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>

#ifdef USE_TBB
#include <tbb/tbb.h>
#endif

#include <spdlog/spdlog.h>

#include <fstream>
#include <iostream>
#include <string>
#include <numeric>

#include <sys/stat.h> // for mkdir

extern const std::string outputFolderPath;

extern std::ofstream logFile;
extern Timer timer, timer_step, timer_temp, timer_temp2, timer_temp3, timer_mt;

extern Eigen::MatrixXi SF;
extern std::vector<int> sTri2Tet;
extern std::vector<bool> isSurfNode;
extern std::vector<int> tetIndToSurf;
extern std::vector<int> surfIndToTet;
extern Eigen::MatrixXd V_surf;
extern Eigen::MatrixXi F_surf;

extern std::vector<bool> isCENode;
extern std::vector<int> tetIndToCE;
extern std::vector<int> CEIndToTet;
extern Eigen::MatrixXd V_CE;
extern Eigen::MatrixXi F_CE;

extern std::vector<int> compVAccSize;
extern std::vector<int> compFAccSize;

extern int numOfCCDFail;

namespace IPC {

template <int dim>
Optimizer<dim>::Optimizer(const Mesh<dim>& p_data0,
    const std::vector<Energy<dim>*>& p_energyTerms, const std::vector<double>& p_energyParams,
    bool p_mute,
    const Eigen::MatrixXd& UV_bnds, const Eigen::MatrixXi& E, const Eigen::VectorXi& bnd,
    const Config& p_animConfig)
    : data0(p_data0), energyTerms(p_energyTerms), energyParams(p_energyParams), animConfig(p_animConfig), OSQPSolver(false)
{
    assert(energyTerms.size() == energyParams.size());

    gradient_ET.resize(energyTerms.size());
    energyVal_ET.resize(energyTerms.size());

    mute = p_mute;

    if (!mute) {
        file_iterStats.open(outputFolderPath + "iterStats.txt");
    }

    if (!data0.checkInversion()) {
        spdlog::error("element inverted in the initial mesh!");
        // exit(-1);
    }
    else {
        spdlog::info("no element inversion detected!");
    }

    globalIterNum = 0;
    relGL2Tol = 1.0e-8;
    innerIterAmt = 0;
    // TODO: Renable the Fischer-Burmeister check
    fbNormTol = std::numeric_limits<double>::infinity(); // 1.0e-6;
    bboxDiagSize2 = data0.matSpaceBBoxSize2(dim);
    dTolRel = 1.0e-9;
    if (animConfig.tuning.size() > 3) {
        dTolRel = animConfig.tuning[3];
    }
    dTol = bboxDiagSize2 * dTolRel * dTolRel;
    rho_DBC = 0.0;

    gravity.setZero();
    if (animConfig.withGravity) {
        gravity[1] = -9.80665;
    }
    setTime(10.0, 0.025);

    spdlog::info("dt and tol initialized");

#ifdef LINSYSSOLVER_USE_CHOLMOD
    linSysSolver = new CHOLMODSolver<Eigen::VectorXi, Eigen::VectorXd>();
    dampingMtr = new CHOLMODSolver<Eigen::VectorXi, Eigen::VectorXd>();
#elif defined(LINSYSSOLVER_USE_AMGCL)
    linSysSolver = new AMGCLSolver<Eigen::VectorXi, Eigen::VectorXd>();
    dampingMtr = new AMGCLSolver<Eigen::VectorXi, Eigen::VectorXd>();
#else
    linSysSolver = new EigenLibSolver<Eigen::VectorXi, Eigen::VectorXd>();
    dampingMtr = new EigenLibSolver<Eigen::VectorXi, Eigen::VectorXd>();
#endif

    setAnimScriptType(animConfig.animScriptType, animConfig.meshSeqFolderPath);
    solveWithQP = (animConfig.isConstrained && (animConfig.constraintSolverType == CST_QP));
    solveWithSQP = (animConfig.isConstrained && (animConfig.constraintSolverType == CST_SQP));
    solveIP = (animConfig.isConstrained && (animConfig.constraintSolverType == CST_IP));
    if (solveWithQP || solveWithSQP || solveIP) {
        activeSet.resize(animConfig.collisionObjects.size());
        activeSet_next.resize(animConfig.collisionObjects.size());
        activeSet_lastH.resize(animConfig.collisionObjects.size());
        lambda_lastH.resize(animConfig.collisionObjects.size());
        MMActiveSet.resize(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        MMActiveSet_next.resize(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        MMActiveSet_CCD.resize(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        paraEEMMCVIDSet.resize(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        paraEEeIeJSet.resize(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        buildConstraintStartIndsWithMM(activeSet, MMActiveSet, constraintStartInds);
        MMActiveSet_lastH.resize(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        MMLambda_lastH.resize(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        MMDistCoord.resize(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        MMTanBasis.resize(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        n_collPairs_sum.setZero(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        n_collPairs_max.setZero(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        n_convCollPairs_sum.setZero(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
        n_convCollPairs_max.setZero(animConfig.meshCollisionObjects.size() + animConfig.isSelfCollision);
    }
    n_collPairs_total_max = n_convCollPairs_total_max = 0;

    solveFric = false;
    for (const auto& coI : animConfig.collisionObjects) {
        if (coI->friction > 0.0) {
            solveFric = true;
            break;
        }
    }
    for (const auto& coI : animConfig.meshCollisionObjects) {
        if (coI->friction > 0.0) {
            solveFric = true;
            break;
        }
    }
    if (animConfig.selfFric > 0.0) {
        solveFric = true;
    }

    beta_NM = animConfig.beta;
    gamma_NM = animConfig.gamma;

    result = data0;
    animScripter.initAnimScript(result, animConfig.collisionObjects,
        animConfig.meshCollisionObjects, animConfig.DBCTimeRange, animConfig.NBCTimeRange);
    result.saveBCNodes(outputFolderPath + "/BC_0.txt");
    if (animConfig.restart) {
        std::ifstream in(animConfig.statusPath.c_str());
        assert(in.is_open());
        std::string line;
        while (std::getline(in, line)) {
            std::stringstream ss(line);
            std::string token;
            ss >> token;
            if (token == "timestep") {
                ss >> globalIterNum;
            }
            else if (token == "position") {
                spdlog::info("read restart position");
                int posRows, dim_in;
                ss >> posRows >> dim_in;
                assert(posRows <= result.V.rows());
                assert(dim_in == result.V.cols());
                for (int vI = 0; vI < posRows; ++vI) {
                    in >> result.V(vI, 0) >> result.V(vI, 1);
                    if constexpr (dim == 3) {
                        in >> result.V(vI, 2);
                    }
                }
            }
            else if (token == "velocity") {
                spdlog::info("read restart velocity");
                int velDim;
                ss >> velDim;
                assert(velDim <= result.V.rows() * dim);
                velocity.setZero(result.V.rows() * dim);
                for (int velI = 0; velI < velDim; ++velI) {
                    in >> velocity[velI];
                }
            }
            else if (token == "dx_Elastic") {
                spdlog::info("read restart dx_Elastic");
                int dxERows, dim_in;
                ss >> dxERows >> dim_in;
                assert(dxERows <= result.V.rows());
                assert(dim_in == dim);
                dx_Elastic.setZero(result.V.rows(), dim);
                for (int vI = 0; vI < dxERows; ++vI) {
                    in >> dx_Elastic(vI, 0) >> dx_Elastic(vI, 1);
                    if constexpr (dim == 3) {
                        in >> dx_Elastic(vI, 2);
                    }
                }
            }
        }
        in.close();
        //TODO: load acceleration, also save acceleration in the saveStatus() function
    }
    else {
        animScripter.initVelocity(result, animConfig.scriptParams, velocity);
        dx_Elastic = Eigen::MatrixXd::Zero(result.V.rows(), dim);
        acceleration = Eigen::MatrixXd::Zero(result.V.rows(), dim);
    }
    result.V_prev = result.V;
    computeXTilta();

    svd.resize(result.F.rows());
    F.resize(result.F.rows());
    if (energyTerms[0]->getNeedElemInvSafeGuard()) {
        if (!result.checkInversion()) {
            spdlog::error("scripted motion causes element inversion, end process");
            exit(-1);
        }
    }
    if (solveIP) {
        sh.build(result, result.avgEdgeLen / 3.0);
        if (isIntersected(result, sh, result.V, animConfig)) {
            spdlog::error("intersection detected in initial configuration!");
            exit(-1);
        }
    }
    spdlog::info("animScriptor set");

    updateTargetGRes();
    CN_MBC = std::sqrt(1.0e-4 * bboxDiagSize2 * dtSq);

    spdlog::info("Newton's solver for Backward Euler constructed");

    numOfLineSearch = 0;

    dHatEps = 1.0e-3;
    if (animConfig.tuning.size() > 1) {
        dHatEps = animConfig.tuning[1];
    }
    dHat = bboxDiagSize2 * dHatEps * dHatEps;

    dHatTarget = bboxDiagSize2 * 1.0e-6;
    if (animConfig.tuning.size() > 2) {
        dHatTarget = animConfig.tuning[2] * animConfig.tuning[2] * bboxDiagSize2;
    }

    fricDHatThres = dHat; // enable friction only when collision dHat is smaller than fricDHatThres if it's not enabled initially
    fricDHat0 = bboxDiagSize2 * 1.0e-6 * dtSq; // initial value of fricDHat
    if (animConfig.tuning.size() > 4) {
        fricDHat0 = bboxDiagSize2 * animConfig.tuning[4] * animConfig.tuning[4] * dtSq;
    }
    fricDHatTarget = bboxDiagSize2 * 1.0e-6 * dtSq;
    if (animConfig.tuning.size() > 5) {
        fricDHatTarget = bboxDiagSize2 * animConfig.tuning[5] * animConfig.tuning[5] * dtSq;
    }
    fricDHat = solveFric ? fricDHat0 : -1.0;

    mu_IP = 0.0;
    if (animConfig.tuning.size() > 0) {
        mu_IP = animConfig.tuning[0];
        upperBoundMu(mu_IP);
    }
    if (mu_IP == 0.0) {
        suggestMu(mu_IP);
    }

    // output wire.poly for rendering
    Eigen::MatrixXd V_surf(result.SVI.size(), 3);
    std::unordered_map<int, int> vI2SVI;
    for (int svI = 0; svI < result.SVI.size(); ++svI) {
        vI2SVI[result.SVI[svI]] = svI;
        V_surf.row(svI) = result.V.row(result.SVI[svI]);
    }

    Eigen::MatrixXi F_surf(result.SF.rows(), 3);
    for (int sfI = 0; sfI < result.SF.rows(); ++sfI) {
        F_surf(sfI, 0) = vI2SVI[result.SF(sfI, 0)];
        F_surf(sfI, 1) = vI2SVI[result.SF(sfI, 1)];
        F_surf(sfI, 2) = vI2SVI[result.SF(sfI, 2)];
    }

    FILE* out = fopen((outputFolderPath + "wire.poly").c_str(), "w");
    assert(out);
    fprintf(out, "POINTS\n");
    for (int vI = 0; vI < V_surf.rows(); ++vI) {
        fprintf(out, "%d: %le %le %le\n", vI + 1, V_surf(vI, 0),
            V_surf(vI, 1), V_surf(vI, 2));
    }
    fprintf(out, "POLYS\n");
    for (int fI = 0; fI < F_surf.rows(); ++fI) {
        int indStart = fI * 3;
        fprintf(out, "%d: %d %d\n", indStart + 1, F_surf(fI, 0) + 1, F_surf(fI, 1) + 1);
        fprintf(out, "%d: %d %d\n", indStart + 2, F_surf(fI, 1) + 1, F_surf(fI, 2) + 1);
        fprintf(out, "%d: %d %d\n", indStart + 3, F_surf(fI, 2) + 1, F_surf(fI, 0) + 1);
    }
    fprintf(out, "END\n");
    fclose(out);
}

template <int dim>
Optimizer<dim>::~Optimizer(void)
{
    if (file_iterStats.is_open()) {
        file_iterStats.close();
    }
    if (file_sysE.is_open()) {
        file_sysE.close();
    }
    if (file_sysM.is_open()) {
        file_sysM.close();
    }
    if (file_sysL.is_open()) {
        file_sysL.close();
    }
    delete linSysSolver;
    delete dampingMtr;
}

template <int dim>
Mesh<dim>& Optimizer<dim>::getResult(void)
{
    return result;
}

template <int dim>
int Optimizer<dim>::getIterNum(void) const
{
    return globalIterNum;
}
template <int dim>
int Optimizer<dim>::getInnerIterAmt(void) const
{
    return innerIterAmt;
}

template <int dim>
void Optimizer<dim>::setRelGL2Tol(double p_relTol)
{
    assert(p_relTol > 0.0);
    relGL2Tol = p_relTol * p_relTol;
    updateTargetGRes();
    logFile << globalIterNum << "th tol: " << targetGRes << std::endl;
}

template <int dim>
double Optimizer<dim>::getDt(void) const
{
    return dt;
}

template <int dim>
void Optimizer<dim>::setAnimScriptType(AnimScriptType animScriptType,
    const std::string& meshSeqFolderPath)
{
    animScripter.setAnimScriptType(animScriptType);
    if (animScriptType == AST_MESHSEQ_FROMFILE) {
        animScripter.setMeshSeqFolderPath(meshSeqFolderPath);
    }
}

template <int dim>
void Optimizer<dim>::setTime(double duration, double dt)
{
    this->dt = dt;
    dtSq = dt * dt;
    frameAmt = std::ceil(duration / dt);
    updateTargetGRes();
    gravityDtSq = dtSq * gravity;
    computeXTilta();
}

//    void Optimizer<dim>::fixDirection(void)
//    {
//        assert(result.V.rows() == result.V_rest.rows());
//
//        directionFix.clear();
//        const Eigen::RowVector2d& v0 = result.V.row(0);
//        int nbVI = *result.vNeighbor[0].begin();
//        const Eigen::RowVector2d& vi = result.V.row(nbVI);
//        Eigen::RowVector2d dif = vi - v0;
//        if(std::abs(dif[0]) > std::abs(dif[1])) {
//            double coef = dif[1] / dif[0];
//            directionFix[0] = coef;
//            directionFix[1] = -1.0;
//            directionFix[nbVI * 2] = -coef;
//            directionFix[nbVI * 2 + 1] = 1.0;
//        }
//        else {
//            double coef = dif[0] / dif[1];
//            directionFix[0] = -1.0;
//            directionFix[1] = coef;
//            directionFix[nbVI * 2] = 1.0;
//            directionFix[nbVI * 2 + 1] = -coef;
//        }
//    }

template <int dim>
void Optimizer<dim>::precompute(void)
{
    spdlog::info("precompute: start");
    if (!mute) { timer_step.start(1); }
    linSysSolver->set_pattern(result.vNeighbor, result.fixedVert);
    if (animConfig.dampingStiff) {
        dampingMtr->set_pattern(result.vNeighbor, result.fixedVert);
    }
    if (!mute) { timer_step.stop(); }
    spdlog::info("precompute: sparse matrix allocated");
    if (solveIP) {
        computeConstraintSets(result);
    }
    computePrecondMtr(result, true, linSysSolver, animConfig.dampingStiff);
    spdlog::info("precompute: sparse matrix entry computed");
    if (!mute) { timer_step.start(2); }
    linSysSolver->analyze_pattern();
    if (!mute) { timer_step.stop(); }
    spdlog::info("precompute: pattern analyzed");

    if (solveWithQP || solveWithSQP) {
        precomputeQPObjective(linSysSolver, P_QP, elemPtr_P_QP);
    }

    computeEnergyVal(result, false, lastEnergyVal);
    if (!mute) {
        spdlog::info("E_initial = {:g}", lastEnergyVal);
    }

    std::vector<double> sysE;
    std::vector<Eigen::Matrix<double, 1, dim>> sysM, sysL;
    computeSystemEnergy(sysE, sysM, sysL);

    file_sysE.open(outputFolderPath + "sysE.txt");
    for (const auto& i : sysE) {
        file_sysE << i << " ";
    }
    file_sysE << std::endl;

    file_sysM.open(outputFolderPath + "sysM.txt");
    for (const auto& i : sysM) {
        file_sysM << i << " ";
    }
    file_sysM << std::endl;

    file_sysL.open(outputFolderPath + "sysL.txt");
    for (const auto& i : sysL) {
        file_sysL << i << " ";
    }
    file_sysL << std::endl;
}

template <int dim>
int Optimizer<dim>::solve(int maxIter)
{
    static bool lastPropagate = false;
    int returnFlag = 0;
    for (int iterI = 0; iterI < maxIter; iterI++) {
        if (!mute) { timer.start(0); }
        timer_step.start(11);
        if (energyTerms[0]->getNeedElemInvSafeGuard()) {
            if (!result.checkInversion()) {
                spdlog::error("before scripted motion there's element inversion, end process");
                exit(-1);
            }
        }

        if (animScripter.stepAnimScript(result, sh,
                animConfig.collisionObjects, animConfig.meshCollisionObjects,
                animConfig.exactCCDMethod, dt, dHat, energyTerms, animConfig.isSelfCollision,
                /*forceIntersectionLineSearch=*/solveIP)) {
            result.saveBCNodes(outputFolderPath + "/BC_" + std::to_string(globalIterNum) + ".txt");
            updatePrecondMtrAndFactorize();
        }
#ifdef CHECK_ET_INTERSECTION
        sh.build(result, result.avgEdgeLen / 3.0);
        bool intersected = false;
        for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
            intersected |= !animConfig.meshCollisionObjects[coI]->checkEdgeTriIntersection(result, sh);
        }
        if (animConfig.isSelfCollision) {
            intersected |= !SelfCollisionHandler<dim>::checkEdgeTriIntersection(result, sh);
        }
        if (intersected) {
            getchar();
        }
#endif

        if (energyTerms[0]->getNeedElemInvSafeGuard()) {
            if (!result.checkInversion()) {
                spdlog::error("scripted motion causes element inversion, end process");
                exit(-1);
            }
        }
        timer_step.stop();

        if (globalIterNum >= frameAmt) {
            // converged
            if (!mute) { timer.stop(); }
            return 1;
        }
        else {
            if (solveIP) {
                if (fullyImplicit_IP()) {
                    returnFlag = 2;
                }
            }
            else {
                if (fullyImplicit()) {
                    returnFlag = 2;
                }
            }

            timer_step.start(11);
            switch (animConfig.timeIntegrationType) {
            case TIT_BE:
                dx_Elastic = result.V - xTilta;
                velocity = Eigen::Map<Eigen::MatrixXd>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>((result.V - result.V_prev).array() / dt).data(), velocity.rows(), 1);
                result.V_prev = result.V;
                computeXTilta();
                break;

            case TIT_NM:
                dx_Elastic = result.V - xTilta;
                velocity = velocity + dt * (1 - gamma_NM) * Eigen::Map<Eigen::MatrixXd>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(acceleration).data(), velocity.rows(), 1);
                acceleration = (result.V - xTilta) / (dtSq * beta_NM);
                acceleration.rowwise() += gravity.transpose();
                velocity += dt * gamma_NM * Eigen::Map<Eigen::MatrixXd>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(acceleration).data(), velocity.rows(), 1);
                result.V_prev = result.V;
                computeXTilta();
                break;
            }

            if (animConfig.dampingStiff) {
                computeDampingMtr(result, false, dampingMtr); //TODO: projectDBC?
            }
            timer_step.stop();
        }
        globalIterNum++;
        if (!mute) { timer.stop(); }
    }
    return returnFlag;
}

template <int dim>
void Optimizer<dim>::updatePrecondMtrAndFactorize(void)
{
    if (!mute) {
        spdlog::info("recompute proxy/Hessian matrix and factorize...");
    }

    if (!mute) { timer_step.start(1); }
    linSysSolver->set_pattern(solveIP ? vNeighbor_IP : result.vNeighbor, result.fixedVert);
    if (animConfig.dampingStiff) {
        dampingMtr->set_pattern(solveIP ? vNeighbor_IP : result.vNeighbor, result.fixedVert);
    }
    if (!mute) { timer_step.start(2); }
    linSysSolver->analyze_pattern();
    if (!mute) { timer_step.start(3); }

    computePrecondMtr(result, true, linSysSolver, animConfig.dampingStiff);

    if (!mute) { timer_step.start(3); }
    linSysSolver->factorize();
    if (!mute) { timer_step.stop(); }
}

template <int dim>
void Optimizer<dim>::precomputeQPObjective(
    const LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSys,
    Eigen::SparseMatrix<double>& P,
    std::vector<double*>& elemPtr_P) const
{
    // Construct matrix sparsity pattern using triplet
    std::vector<Eigen::Triplet<double>> coefs;
    coefs.reserve(linSys->getNumNonzeros() * 2 - linSys->getNumRows());
    for (int rowI = 0; rowI < linSys->getNumRows(); rowI++) {
        for (const auto& colIter : linSys->getIJ2aI()[rowI]) {
            coefs.emplace_back(rowI, colIter.first, 0.0);
            if (rowI != colIter.first) {
                coefs.emplace_back(colIter.first, rowI, 0.0);
            }
        }
    }
    P.resize(linSys->getNumRows(), linSys->getNumRows());
    P.setZero();
    P.reserve(coefs.size());
    P.setFromTriplets(coefs.begin(), coefs.end());

    // get matrix entry pointer for fast access
    elemPtr_P.resize(0);
    elemPtr_P.reserve(P.nonZeros());
    for (int rowI = 0; rowI < linSys->getNumRows(); rowI++) {
        for (const auto& colIter : linSys->getIJ2aI()[rowI]) {
            elemPtr_P.emplace_back(&P.coeffRef(rowI, colIter.first));
            if (rowI != colIter.first) {
                elemPtr_P.emplace_back(&P.coeffRef(colIter.first, rowI));
            }
        }
    }
}

template <int dim>
void Optimizer<dim>::updateQPObjective(
    const LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSys,
    const std::vector<double*>& elemPtr_P) const
{
    // Update Hessian using entry pointers
    int elemPtrI = 0;
    for (int rowI = 0; rowI < linSys->getNumRows(); rowI++) {
        for (const auto& colIter : linSys->getIJ2aI()[rowI]) {
            assert(elemPtrI < elemPtr_P.size());
            *elemPtr_P[elemPtrI++] = linSys->get_a()[colIter.second];
            if (rowI != colIter.first) {
                assert(elemPtrI < elemPtr_P.size());
                *elemPtr_P[elemPtrI++] = linSys->get_a()[colIter.second];
            }
        }
    }
}

template <int dim>
void Optimizer<dim>::computeQPInequalityConstraint(
    const Mesh<dim>& mesh,
    const std::vector<CollisionObject<dim>*>& collisionObjects,
    const std::vector<std::vector<int>>& activeSet,
    const int num_vars,
    std::vector<int>& constraintStartInds,
    Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) const
{
    // construct constraint matrix and vector
    std::vector<Eigen::Triplet<double>> A_triplet;

    constraintStartInds.resize(1);
    constraintStartInds[0] = 0;
    for (int coI = 0; coI < collisionObjects.size(); ++coI) {
        collisionObjects[coI]->updateConstraints_QP(
            mesh, activeSet[coI], A_triplet, b);
        constraintStartInds.emplace_back(b.size());
    }
    //TODO: use input argument
    for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
        animConfig.meshCollisionObjects[coI]->updateConstraints_QP(
            mesh, MMActiveSet[coI], animConfig.constraintType, A_triplet, b);
        constraintStartInds.emplace_back(b.size());
    }
    if (animConfig.isSelfCollision) {
        SelfCollisionHandler<dim>::updateConstraints_QP(
            mesh, MMActiveSet.back(), animConfig.constraintType,
            mesh_mmcvid_to_toi, A_triplet, b);
        constraintStartInds.emplace_back(b.size());
    }

    // Offset the constraints to allow for solver tolerances.
    b.array() += animConfig.constraintOffset;

    int num_constraints = b.size();
    A.resize(num_constraints, num_vars);
    if (num_constraints > 0) {
        A.setZero();
        A.reserve(A_triplet.size());
        A.setFromTriplets(A_triplet.begin(), A_triplet.end());
    }
}

template <int dim>
bool Optimizer<dim>::solveQP(
    const Mesh<dim>& mesh,
    const std::vector<CollisionObject<dim>*>& collisionObjects,
    const std::vector<std::vector<int>>& activeSet,
    const LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSys,
    Eigen::SparseMatrix<double>& P,
    const std::vector<double*>& elemPtr_P,
    Eigen::VectorXd& gradient,
    OSQP& OSQPSolver,
#ifdef USE_GUROBI
    Eigen::GurobiSparse& gurobiQPSolver,
#endif
    std::vector<int>& constraintStartInds,
    Eigen::VectorXd& searchDir,
    Eigen::VectorXd& dual,
    QPSolverType qpSolverType) const
{
    // minₓ ½xᵀQx + cᵀx
    if (!solveWithQP) {
        this->updateQPObjective(linSys, elemPtr_P);
    }
    // Ax ≥ b
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    this->computeQPInequalityConstraint(
        mesh, collisionObjects, activeSet, linSys->getNumRows(),
        constraintStartInds, A, b);

    switch (qpSolverType) {
    case QPSolverType::QP_OSQP:
        return solveQP_OSQP(
            P, gradient, A, b, OSQPSolver, searchDir, dual);
    case QPSolverType::QP_GUROBI:
#ifdef USE_GUROBI
        return solveQP_Gurobi(
            P, gradient, A, b, gurobiQPSolver, searchDir, dual);
#else
        spdlog::error("Gurobi QP solve is disabled. Set IPC_WITH_GUROBI=ON in cmake to enable.");
        throw "Gurobi QP solve is disabled. Set IPC_WITH_GUROBI=ON in cmake to enable.";
#endif
    }

    // std::string folderPath = outputFolderPath + std::to_string(globalIterNum) + "/";
    // mkdir(folderPath.c_str(), 0777);
    // IglUtils::writeSparseMatrixToFile(folderPath + "P" + std::to_string(innerIterAmt), Q, true);
    // IglUtils::writeVectorToFile(folderPath + "q" + std::to_string(innerIterAmt), c);
    // IglUtils::writeSparseMatrixToFile(folderPath + "A" + std::to_string(innerIterAmt), A, true);
    // IglUtils::writeVectorToFile(folderPath + "l" + std::to_string(innerIterAmt), b);
    // IglUtils::writeVectorToFile(folderPath + "primal" + std::to_string(innerIterAmt), x);
    // IglUtils::writeVectorToFile(folderPath + "dual" + std::to_string(innerIterAmt), dual);
}

template <int dim>
bool Optimizer<dim>::solveQP_OSQP(
    Eigen::SparseMatrix<double>& Q, Eigen::VectorXd& c,
    Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b,
    OSQP& QPSolver, Eigen::VectorXd& x, Eigen::VectorXd& dual) const
{
    c_int num_vars = c_int(c.rows()), num_constraints = c_int(b.size());
    Eigen::VectorXd u = Eigen::VectorXd::Constant(
        num_constraints, std::numeric_limits<double>::infinity());

    // setup QPSolver and solve
    QPSolver.setup(
        Q.valuePtr(), c_int(Q.nonZeros()), Q.innerIndexPtr(), Q.outerIndexPtr(),
        c.data(),
        A.valuePtr(), c_int(A.nonZeros()),
        A.innerIndexPtr(), A.outerIndexPtr(),
        b.data(), u.data(),
        num_vars, num_constraints);

    x.conservativeResize(num_vars);
    c_int status = QPSolver.solve();
    if ((status != -3) && (status != -4)) {
        memcpy(x.data(), QPSolver.getPrimal(),
            x.size() * sizeof(x[0]));
    }

    if (num_constraints > 0) {
        dual.conservativeResize(num_constraints);
        if ((status != -3) && (status != -4)) {
            memcpy(dual.data(), QPSolver.getDual(),
                dual.size() * sizeof(dual[0]));
        }
    }

    return (status != -3) && (status != -4);
}

#ifdef USE_GUROBI
template <int dim>
bool Optimizer<dim>::solveQP_Gurobi(
    const Eigen::SparseMatrix<double>& Q, const Eigen::VectorXd& c,
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b,
    Eigen::GurobiSparse& QPSolver,
    Eigen::VectorXd& x, Eigen::VectorXd& dual) const
{
    int num_vars = c.rows(), num_constraints = b.size();
    // setup QPSolver and solve
    QPSolver.displayOutput(false);
    QPSolver.problem(num_vars, /*nreq=*/0, num_constraints);

    Eigen::SparseVector<double> inf = Eigen::VectorXd::Constant(
        num_vars, std::numeric_limits<double>::infinity())
                                          .sparseView();

    bool success;
    try {
        success = QPSolver.solve(
            /*Q=*/Q, /*C=*/c.sparseView(),
            /*Aeq=*/Eigen::SparseMatrix<double>(), /*Beq=*/Eigen::SparseVector<double>(),
            // We want Ax ≥ b, but Gurobi expects Ax ≤ b
            /*Aineq=*/-A, /*Bineq=*/-b.sparseView(),
            /*XL=*/-inf, /*XU=*/inf);
    }
    catch (GRBException e) {
        spdlog::error("GRBException: {:s}", e.getMessage());
        success = false;
    }

    x.conservativeResize(num_vars);
    if (success) {
        memcpy(x.data(), QPSolver.result().data(), x.size() * sizeof(x[0]));
    }
    else {
        x.setZero();
    }

    if (num_constraints > 0) {
        dual.conservativeResize(num_constraints);
        if (success) {
            memcpy(dual.data(), QPSolver.dual_ineq().data(),
                dual.size() * sizeof(dual[0]));
        }
    }

    return success;
}
#endif

template <int dim>
void Optimizer<dim>::computeQPResidual(const Mesh<dim>& mesh,
    const std::vector<CollisionObject<dim>*>& collisionObjects,
    const std::vector<std::vector<int>>& activeSet,
    const std::vector<int>& constraintStartInds,
    const Eigen::VectorXd& gradient,
    const Eigen::VectorXd& dual,
    Eigen::VectorXd& grad_KKT,
    Eigen::VectorXd& constraintVal,
    Eigen::VectorXd& fb)
{
    grad_KKT = gradient;
    constraintVal.conservativeResize(0);
    for (int coI = 0; coI < collisionObjects.size(); ++coI) {
        int constraintAmtI = constraintStartInds[coI + 1] - constraintStartInds[coI];
        if (constraintAmtI) {
            collisionObjects[coI]->evaluateConstraintsQP(mesh, activeSet[coI], constraintVal);
            collisionObjects[coI]->leftMultiplyConstraintJacobianTQP(mesh, activeSet[coI],
                dual.segment(constraintStartInds[coI], constraintAmtI),
                grad_KKT);
        }
    }

    for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
        int constraintAmtI = constraintStartInds[coI + 1 + collisionObjects.size()] - constraintStartInds[coI + collisionObjects.size()];
        if (constraintAmtI) {
            animConfig.meshCollisionObjects[coI]->evaluateConstraintsQP(
                mesh, MMActiveSet[coI], animConfig.constraintType, constraintVal);
            animConfig.meshCollisionObjects[coI]->leftMultiplyConstraintJacobianTQP(
                mesh, MMActiveSet[coI],
                dual.segment(constraintStartInds[coI + collisionObjects.size()], constraintAmtI),
                animConfig.constraintType,
                grad_KKT);
        }
    }

    if (animConfig.isSelfCollision) {
        int constraintStartI = constraintStartInds[constraintStartInds.size() - 2];
        int constraintAmtI = constraintStartInds.back() - constraintStartI;
        if (constraintAmtI) {
            SelfCollisionHandler<dim>::evaluateConstraintsQP(
                mesh, MMActiveSet.back(), animConfig.constraintType,
                mesh_mmcvid_to_toi, constraintVal);
            SelfCollisionHandler<dim>::leftMultiplyConstraintJacobianTQP(
                mesh, MMActiveSet.back(),
                dual.segment(constraintStartI, constraintAmtI),
                animConfig.constraintType, mesh_mmcvid_to_toi, grad_KKT);
        }
    }

    fb.conservativeResize(constraintStartInds.back()); // Fischer-Burmeister measure for complementarity
    for (int constraintI = 0; constraintI < constraintStartInds.back(); ++constraintI) {
        const double dualI = -dual[constraintI]; // OSQP produce negative dual for >= constraints
        const double& constraint = constraintVal[constraintI];
        fb[constraintI] = (dualI + constraint - std::sqrt(dualI * dualI + constraint * constraint));
    }
}

template <int dim>
void Optimizer<dim>::initX(int option, std::vector<std::vector<int>>& p_activeSet_next)
{
    // global:
    searchDir.conservativeResize(result.V.rows() * dim);
    switch (option) {
    case 0:
        // already at last timestep config
        searchDir.setZero();
        break;

    case 1: // explicit Euler
#ifdef USE_TBB
        tbb::parallel_for(0, (int)result.V.rows(), 1, [&](int vI)
#else
        for (int vI = 0; vI < result.V.rows(); vI++)
#endif
            {
                if (result.isFixedVert[vI]) {
                    searchDir.segment<dim>(vI * dim).setZero();
                }
                else {
                    searchDir.segment<dim>(vI * dim) = dt * velocity.segment<dim>(vI * dim);
                }
            }
#ifdef USE_TBB
        );
#endif
        break;

    case 2: // xHat
        switch (animConfig.timeIntegrationType) {
        case TIT_BE:
#ifdef USE_TBB
            tbb::parallel_for(0, (int)result.V.rows(), 1, [&](int vI)
#else
            for (int vI = 0; vI < result.V.rows(); vI++)
#endif
                {
                    if (result.isFixedVert[vI]) {
                        searchDir.segment<dim>(vI * dim).setZero();
                    }
                    else {
                        searchDir.segment<dim>(vI * dim) = dt * velocity.segment<dim>(vI * dim) + gravityDtSq;
                    }
                }
#ifdef USE_TBB
            );
#endif
            break;

        case TIT_NM:
#ifdef USE_TBB
            tbb::parallel_for(0, (int)result.V.rows(), 1, [&](int vI)
#else
            for (int vI = 0; vI < result.V.rows(); vI++)
#endif
                {
                    if (result.isFixedVert[vI]) {
                        searchDir.segment<dim>(vI * dim).setZero();
                    }
                    else {
                        searchDir.segment<dim>(vI * dim) = dt * velocity.segment<dim>(vI * dim) + gravityDtSq / 2.0;
                    }
                }
#ifdef USE_TBB
            );
#endif
            break;
        }
        break;

    case 3: { // Symplectic Euler
        switch (animConfig.timeIntegrationType) {
        case TIT_BE:
#ifdef USE_TBB
            tbb::parallel_for(0, (int)result.V.rows(), 1, [&](int vI)
#else
            for (int vI = 0; vI < result.V.rows(); vI++)
#endif
                {
                    if (result.isFixedVert[vI]) {
                        searchDir.segment<dim>(vI * dim).setZero();
                    }
                    else {
                        searchDir.segment<dim>(vI * dim) = (dt * velocity.segment<dim>(vI * dim) + gravityDtSq + dx_Elastic.row(vI).transpose());
                    }
                }
#ifdef USE_TBB
            );
#endif
            break;

        case TIT_NM:
#ifdef USE_TBB
            tbb::parallel_for(0, (int)result.V.rows(), 1, [&](int vI)
#else
            for (int vI = 0; vI < result.V.rows(); vI++)
#endif
                {
                    if (result.isFixedVert[vI]) {
                        searchDir.segment<dim>(vI * dim).setZero();
                    }
                    else {
                        searchDir.segment<dim>(vI * dim) = dt * velocity.segment<dim>(vI * dim) + gravityDtSq / 2.0 + dx_Elastic.row(vI).transpose() * 2.0;
                    }
                }
#ifdef USE_TBB
            );
#endif
            break;
        }
        break;
    }

    case 4: { // uniformly accelerated motion approximation
        switch (animConfig.timeIntegrationType) {
        case TIT_BE:
#ifdef USE_TBB
            tbb::parallel_for(0, (int)result.V.rows(), 1, [&](int vI)
#else
            for (int vI = 0; vI < result.V.rows(); vI++)
#endif
                {
                    if (result.isFixedVert[vI]) {
                        searchDir.segment<dim>(vI * dim).setZero();
                    }
                    else {
                        searchDir.segment<dim>(vI * dim) = (dt * velocity.segment<dim>(vI * dim) + (gravityDtSq + 0.5 * dx_Elastic.row(vI).transpose()));
                    }
                }
#ifdef USE_TBB
            );
#endif
            break;

        case TIT_NM:
#ifdef USE_TBB
            tbb::parallel_for(0, (int)result.V.rows(), 1, [&](int vI)
#else
            for (int vI = 0; vI < result.V.rows(); vI++)
#endif
                {
                    if (result.isFixedVert[vI]) {
                        searchDir.segment<dim>(vI * dim).setZero();
                    }
                    else {
                        searchDir.segment<dim>(vI * dim) = dt * velocity.segment<dim>(vI * dim) + gravityDtSq / 2.0 + dx_Elastic.row(vI).transpose();
                    }
                }
#ifdef USE_TBB
            );
#endif
            break;
        }
        break;
    }

    case 5: { // Jacobi
        Eigen::VectorXd g;
        computeGradient(result, true, g);

        Eigen::VectorXi I, J;
        Eigen::VectorXd V;
        computePrecondMtr(result, false, linSysSolver);

#ifdef USE_TBB
        tbb::parallel_for(0, (int)result.V.rows(), 1, [&](int vI)
#else
        for (int vI = 0; vI < result.V.rows(); vI++)
#endif
            {
                if (result.isFixedVert[vI]) {
                    searchDir.segment<dim>(vI * dim).setZero();
                }
                else {
                    searchDir[vI * dim] = -g[vI * dim] / linSysSolver->coeffMtr(vI * dim, vI * dim);
                    searchDir[vI * dim + 1] = -g[vI * dim + 1] / linSysSolver->coeffMtr(vI * dim + 1, vI * dim + 1);
                    searchDir[vI * dim + 2] = -g[vI * dim + 2] / linSysSolver->coeffMtr(vI * dim + 2, vI * dim + 2);
                }
            }
#ifdef USE_TBB
        );
#endif

        break;
    }

    default:
        spdlog::warn("unkown primal initialization type, use last timestep instead");
        break;
    }

    if (option) {
        double stepSize = 1.0, slackness_a = 0.9, slackness_m = 0.8;
        energyTerms[0]->filterStepSize(result, searchDir, stepSize);

        if (solveWithQP || solveWithSQP) {
            // Always clear the active set at the start of time-step
            clearActiveSet();
            updateActiveSet_QP();
        }
        else if (solveIP) {
            for (int coI = 0; coI < animConfig.collisionObjects.size(); coI++) {
                p_activeSet_next[coI].resize(0);
                animConfig.collisionObjects[coI]->largestFeasibleStepSize(result, searchDir, slackness_a,
                    p_activeSet_next[coI], stepSize);
            }
            // always use full CCD for warm start
#ifdef USE_SH_LFSS
            timer_temp3.start(11);
            sh.build(result, searchDir, stepSize, result.avgEdgeLen / 3.0);
            timer_temp3.stop();
#endif
            for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); coI++) {
                MMActiveSet_next[coI].resize(0);
                if (animConfig.exactCCDMethod == ExactCCD::Method::NONE) {
                    animConfig.meshCollisionObjects[coI]->largestFeasibleStepSize_CCD(result, sh, searchDir, slackness_m, stepSize);
                }
                else {
                    animConfig.meshCollisionObjects[coI]->largestFeasibleStepSize_CCD_exact(result, sh, searchDir, animConfig.exactCCDMethod, stepSize);
                }
            }
            if (animConfig.isSelfCollision) {
                MMActiveSet_next.back().resize(0);
                if (animConfig.exactCCDMethod == ExactCCD::Method::NONE) {
                    std::vector<std::pair<int, int>> newCandidates;
                    SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD(result, sh, searchDir, slackness_m, newCandidates, stepSize);
                }
                else {
                    SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD_exact(result, sh, searchDir, animConfig.exactCCDMethod, stepSize);
                }
            }

            if (stepSize == 0.0) {
                spdlog::info("CCD gives 0 in initX()");
            }
        }

        if (!(solveWithQP || solveWithSQP)) {
            Eigen::MatrixXd V0 = result.V;
            stepForward(V0, result, stepSize);

            if (energyTerms[0]->getNeedElemInvSafeGuard()) {
                while (!result.checkInversion(true)) {
                    logFile << "element inversion detected during warmStart, backtrack!" << std::endl;
                    stepSize /= 2.0;
                    stepForward(V0, result, stepSize);
                }
            }

            if (solveIP) {
                sh.build(result, result.avgEdgeLen / 3.0);
                while (isIntersected(result, sh, V0, animConfig)) {
                    logFile << "intersection detected during warmStart, backtrack!" << std::endl;
                    stepSize /= 2.0;
                    stepForward(V0, result, stepSize);
                    sh.build(result, result.avgEdgeLen / 3.0);
                }
            }
        }
#ifdef CHECK_ET_INTERSECTION
        sh.build(result, result.avgEdgeLen / 3.0);
        bool intersected = false;
        for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
            intersected |= !animConfig.meshCollisionObjects[coI]->checkEdgeTriIntersection(result, sh);
        }
        if (animConfig.isSelfCollision) {
            intersected |= !SelfCollisionHandler<dim>::checkEdgeTriIntersection(result, sh);
        }
        if (intersected) {
            getchar();
        }
#endif

        spdlog::debug(
            "timestep={:d} init_stepSize={:g}", globalIterNum, stepSize);
    }
    else if (solveWithQP || solveWithSQP || solveIP) {
        clearActiveSet();
    }
}

template <int dim>
void Optimizer<dim>::computeXTilta(void)
{
    xTilta.conservativeResize(result.V.rows(), dim);
    switch (animConfig.timeIntegrationType) {
    case TIT_BE:
#ifdef USE_TBB
        tbb::parallel_for(0, (int)result.V.rows(), 1, [&](int vI)
#else
        for (int vI = 0; vI < result.V.rows(); vI++)
#endif
            {
                if (result.isFixedVert[vI]) {
                    xTilta.row(vI) = result.V_prev.row(vI);
                }
                else {
                    xTilta.row(vI) = (result.V_prev.row(vI) + (velocity.segment<dim>(vI * dim) * dt + gravityDtSq).transpose());
                }
            }
#ifdef USE_TBB
        );
#endif
        break;

    case TIT_NM:
#ifdef USE_TBB
        tbb::parallel_for(0, (int)result.V.rows(), 1, [&](int vI)
#else
        for (int vI = 0; vI < result.V.rows(); vI++)
#endif
            {
                if (result.isFixedVert[vI]) {
                    xTilta.row(vI) = result.V_prev.row(vI);
                }
                else {
                    xTilta.row(vI) = (result.V_prev.row(vI) + (velocity.segment<dim>(vI * dim) * dt + beta_NM * gravityDtSq + (0.5 - beta_NM) * (dtSq * acceleration.row(vI).transpose())).transpose());
                }
            }
#ifdef USE_TBB
        );
#endif
        break;
    }
}

template <int dim>
void Optimizer<dim>::clearActiveSet()
{
    for (int coI = 0; coI < animConfig.collisionObjects.size(); coI++) {
        activeSet_next[coI].clear();
    }
    for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); coI++) {
        MMActiveSet_next[coI].clear();
    }
    if (animConfig.isSelfCollision) {
        MMActiveSet_next.back().clear();
    }
}

template <int dim>
bool Optimizer<dim>::updateActiveSet_QP()
{
    bool newConstraintsAdded = false;
    // TODO: Explain this factor of two
    const double eta = 2 * animConfig.constraintOffset;
    for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
        // TODO: Modify this function to check if the set changed
        animConfig.collisionObjects[coI]->filterSearchDir_QP(
            result, searchDir, activeSet_next[coI]);
        activeSet[coI] = activeSet_next[coI];
    }
    for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); coI++) {
        newConstraintsAdded |= animConfig.meshCollisionObjects[coI]->updateActiveSet_QP(
            result, searchDir, animConfig.constraintType,
            MMActiveSet_next[coI], eta);
        MMActiveSet[coI] = MMActiveSet_next[coI];
    }
    if (animConfig.isSelfCollision) {
        newConstraintsAdded |= SelfCollisionHandler<dim>::updateActiveSet_QP(
            result, searchDir, animConfig.constraintType,
            MMActiveSet_next.back(), mesh_mmcvid_to_toi, eta);
        MMActiveSet.back() = MMActiveSet_next.back();
    }
    return newConstraintsAdded;
}

template <int dim>
int Optimizer<dim>::countConstraints()
{
    int constraintAmt = 0;
    for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
        constraintAmt += activeSet[coI].size();
    }
    for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); coI++) {
        constraintAmt += MMActiveSet[coI].size();
    }
    if (animConfig.isSelfCollision && MMActiveSet.size() > 0) {
        constraintAmt += MMActiveSet.back().size();
    }
    return constraintAmt;
}

template <int dim>
bool Optimizer<dim>::fullyImplicit(void)
{
    timer_step.start(10);

    // initial warm start with inital guess search direction p₀
    // and compute initial active set A_0 checking for collisions from x_t to x_t + p_0
    initX(animConfig.warmStart, activeSet_next);

    computeEnergyVal(result, true, lastEnergyVal);
    computeGradient(result, false, gradient);
    if (solveWithQP) {
        // Compute the hessian once at the beginning of the optimization
        computePrecondMtr(result, false, linSysSolver);
        updateQPObjective(linSysSolver, elemPtr_P_QP);
    }
    timer_step.stop();

    int constraintAmt = countConstraints();
    spdlog::debug("after initX: E={:g} ||g||²={:g} targetGRes={:g} constraintAmt={:d}",
        lastEnergyVal, gradient.squaredNorm(), targetGRes, constraintAmt);
    file_iterStats << globalIterNum << " 0 " << lastEnergyVal << " "
                   << gradient.squaredNorm() << " 0" << std::endl;

    // Termination criteria variables
    Eigen::VectorXd fb;
    double sqn_g = __DBL_MAX__, sqn_g0;
    int iterCap = 10000, curIterI = 0;
    bool newConstraintsAdded = false;

    double meshScale = (result.V.colwise().maxCoeff() - result.V.colwise().minCoeff()).maxCoeff();

    do {
        file_iterStats << globalIterNum << " ";

        if (solve_oneStep()) {
            spdlog::error(
                "timestep={:d} iteration={:d} "
                "msg=\"line search failed!!!\"",
                globalIterNum, curIterI);
            logFile << "\tline search with Armijo's rule failed!!!"
                    << std::endl;
            return true;
        }
        innerIterAmt++;

        if (searchDir.norm() > 100 * meshScale) {
            spdlog::warn(
                "timestep={:d} iteration={:05d} msg=\"possible displacement blow-up\" meshScale={:g} searchDir.norm()={:g}",
                globalIterNum, curIterI, meshScale, searchDir.norm());
            result.saveSurfaceMesh("output/blowup_t" + std::to_string(globalIterNum) + "_i" + std::to_string(curIterI) + ".obj");
#ifdef EXIT_UPON_QP_FAIL
            exit(2);
#endif
        }

        timer_step.start(10);
        Eigen::VectorXd grad_KKT, constraintVal;
        if (solveWithQP || solveWithSQP) {
            computeQPResidual(result, animConfig.collisionObjects, activeSet,
                constraintStartInds, gradient, dual_QP,
                grad_KKT, constraintVal, fb);
            // || H(xᵗ)dx + ∇E(xᵗ) + Cᵀλ ||²
            if (solveWithQP) {
                sqn_g = (P_QP * searchDir + grad_KKT).squaredNorm();
            }
            else {
                sqn_g = grad_KKT.squaredNorm();
            }
            newConstraintsAdded = updateActiveSet_QP();
        }
        else {
            sqn_g = gradient.squaredNorm();
        }

        if (curIterI == 0) {
            sqn_g0 = sqn_g;
        }

        file_iterStats << lastEnergyVal << " " << sqn_g << " " << (constraintStartInds.empty() ? 0 : constraintStartInds.back()) << std::endl;
        timer_step.stop();

        constraintAmt = countConstraints();
        double currentEnergyVal;
        computeEnergyVal(result, true, currentEnergyVal);
        spdlog::debug(
            "timestep={:d} iteration={:05d} currentEnergyVal={:.6e} "
            "gradient.squaredNorm()={:g} sqn_g={:.6e} fb.norm()={:.6e} "
            "constraintAmt={:03d} newConstraintsAdded={}",
            globalIterNum, curIterI, currentEnergyVal, gradient.squaredNorm(),
            sqn_g, fb.norm(), constraintAmt, newConstraintsAdded);

        if (++curIterI >= iterCap) {
            spdlog::warn(
                "timestep={:d} msg=\"optimization did not converge in at most {:d} iterations\" "
                "sqn_g={:g} targetGRes={:g} fb.norm()={:g} fbNormTol={:g}",
                globalIterNum, iterCap, sqn_g, targetGRes, fb.norm(), fbNormTol);
            break;
        }

        double curEnergyVal;
        computeEnergyVal(result, true, curEnergyVal);
        if ((curEnergyVal - lastEnergyVal) / lastEnergyVal > 1e9 && curEnergyVal > 1e9) {
            spdlog::warn(
                "timestep={:d} iteration={:05d} msg=\"possible energy blow-up\" lastEnergyVal={:g} curEnergyVal={:g}",
                globalIterNum, curIterI, lastEnergyVal, curEnergyVal);
            result.saveSurfaceMesh("output/blowup_t" + std::to_string(globalIterNum) + "_i" + std::to_string(curIterI) + ".obj");
#ifdef EXIT_UPON_QP_FAIL
            exit(2);
#endif
        }

        if ((sqn_g - sqn_g0) / sqn_g0 > 1e9 && sqn_g > 1e9) {
            spdlog::warn(
                "timestep={:d} iteration={:05d} msg=\"possible gradient blow-up\" sqn_g0={:g} sqn_g={:g}",
                globalIterNum, curIterI, sqn_g0, sqn_g);
            result.saveSurfaceMesh("output/blowup_t" + std::to_string(globalIterNum) + "_i" + std::to_string(curIterI) + ".obj");
#ifdef EXIT_UPON_QP_FAIL
            exit(2);
#endif
        }
    } while ((animConfig.useActiveSetConvergence && newConstraintsAdded)
        || (sqn_g > targetGRes) || (solveWithSQP && (fb.norm() > fbNormTol)));

    spdlog::info(
        "timestep={:d} total_iterations={:d} sqn_g={:.6e} fb.norm()={:.6e} "
        "dual_QP.maxCoeff()={:s} dual_QP.minCoeff()={:s}",
        globalIterNum, curIterI, sqn_g, fb.norm(),
        dual_QP.size() ? std::to_string(dual_QP.maxCoeff()) : "N/a",
        dual_QP.size() ? std::to_string(dual_QP.minCoeff()) : "N/a");

    sh.build(result, result.avgEdgeLen / 3.0);
    bool isIntersecting = false;
    for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
        if (!animConfig.meshCollisionObjects[coI]->checkEdgeTriIntersection(result, sh)) {
            spdlog::error(
                "timestep={:d} msg=\"intersecting state of the mesh and collision object {:d}\"",
                globalIterNum, coI);
            isIntersecting = true;
        }
    }
    if (animConfig.isSelfCollision && !SelfCollisionHandler<dim>::checkEdgeTriIntersection(result, sh)) {
        spdlog::error(
            "timestep={:d} msg=\"self-intersecting state\"", globalIterNum);
        isIntersecting = true;
    }
    if (isIntersecting) {
        result.saveSurfaceMesh("output/intersecting_t" + std::to_string(globalIterNum) + ".obj");
#ifdef EXIT_UPON_QP_FAIL
        exit(1);
#endif
    }

    logFile << "Timestep" << globalIterNum << " innerIterAmt = " << innerIterAmt << ", accumulated line search steps " << numOfLineSearch << std::endl;

    std::vector<double> sysE;
    std::vector<Eigen::Matrix<double, 1, dim>> sysM, sysL;
    computeSystemEnergy(sysE, sysM, sysL);

    for (const auto& i : sysE) {
        file_sysE << i << " ";
    }
    file_sysE << std::endl;

    for (const auto& i : sysM) {
        file_sysM << i << " ";
    }
    file_sysM << std::endl;

    for (const auto& i : sysL) {
        file_sysL << i << " ";
    }
    file_sysL << std::endl;

    return (curIterI >= iterCap);
}

template <int dim>
bool Optimizer<dim>::fullyImplicit_IP(void)
{
    spdlog::info("Timestep {:d}:", globalIterNum);

    for (auto& i : activeSet_lastH) {
        i.resize(0);
    }
    for (auto& i : MMActiveSet_lastH) {
        i.resize(0);
    }

    timer_step.start(10);
    initX(animConfig.warmStart, activeSet_next);
    timer_step.stop();

    fricDHat = solveFric ? fricDHat0 : -1.0;
    dHat = bboxDiagSize2 * dHatEps * dHatEps;
    computeConstraintSets(result);

    mu_IP = 0.0;
    if (animConfig.tuning.size() > 0) {
        mu_IP = animConfig.tuning[0];
        upperBoundMu(mu_IP);
    }
    if (mu_IP == 0.0) {
        suggestMu(mu_IP);
    }
#ifdef ADAPTIVE_MU
    initMu_IP(mu_IP);
#endif

    timer_step.start(10);
    Eigen::VectorXd constraintVal;
    bool activeSetChanged = false;
    for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
        int startCI = constraintVal.size();
        animConfig.collisionObjects[coI]->evaluateConstraints(result,
            activeSet[coI], constraintVal);

        if (animConfig.collisionObjects[coI]->friction > 0.0) {
            lambda_lastH[coI].resize(constraintVal.size() - startCI);
            //TODO: parallelize
            for (int i = 0; i < lambda_lastH[coI].size(); ++i) {
                compute_g_b(constraintVal[startCI + i], dHat, lambda_lastH[coI][i]);
                lambda_lastH[coI][i] *= -mu_IP * 2.0 * std::sqrt(constraintVal[startCI + i]);
            }

            if (activeSet_lastH[coI] != activeSet[coI]) {
                activeSetChanged = true;
                activeSet_lastH[coI] = activeSet[coI];
            }
        }
    }
    for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
        animConfig.meshCollisionObjects[coI]->evaluateConstraints(result,
            MMActiveSet[coI], constraintVal);
    }
    if (animConfig.isSelfCollision) {
        int startCI = constraintVal.size();
        SelfCollisionHandler<dim>::evaluateConstraints(result, MMActiveSet.back(), constraintVal);

        if (animConfig.selfFric > 0.0) {
            MMLambda_lastH.back().resize(constraintVal.size() - startCI);
            //TODO: parallelize
            for (int i = 0; i < MMLambda_lastH.back().size(); ++i) {
                compute_g_b(constraintVal[startCI + i], dHat, MMLambda_lastH.back()[i]);
                MMLambda_lastH.back()[i] *= -mu_IP * 2.0 * std::sqrt(constraintVal[startCI + i]);
                if (MMActiveSet.back()[i][3] < -1) {
                    // PP or PE duplication
                    MMLambda_lastH.back()[i] *= -MMActiveSet.back()[i][3];
                }
            }

            SelfCollisionHandler<dim>::computeDistCoordAndTanBasis(
                result, MMActiveSet.back(), MMDistCoord.back(), MMTanBasis.back());
            if (MMActiveSet_lastH.back() != MMActiveSet.back()) {
                activeSetChanged = true;
                MMActiveSet_lastH.back() = MMActiveSet.back();
            }
        }
    }
    timer_step.stop();

    timer_step.start(10);
    computeEnergyVal(result, true, lastEnergyVal); // for possible warmstart and pre-scripting Dirichlet boundary conditions
    file_iterStats << globalIterNum << " 0.0 0" << std::endl;
    // globaIterNum lineSearchStepSize constraintNum
    timer_step.stop();
    int fricIterI = 0;
    while (true) {
        initSubProb_IP();
        solveSub_IP(mu_IP, activeSet_next, MMActiveSet_next);
        ++fricIterI;

        timer_step.start(10);
        Eigen::VectorXd constraintVal;
        bool activeSetChanged = false;
        for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
            int startCI = constraintVal.size();
            animConfig.collisionObjects[coI]->evaluateConstraints(result,
                activeSet[coI], constraintVal);

            if (animConfig.collisionObjects[coI]->friction > 0.0) {
                lambda_lastH[coI].resize(constraintVal.size() - startCI);
                //TODO: parallelize
                for (int i = 0; i < lambda_lastH[coI].size(); ++i) {
                    compute_g_b(constraintVal[startCI + i], dHat, lambda_lastH[coI][i]);
                    lambda_lastH[coI][i] *= -mu_IP * 2.0 * std::sqrt(constraintVal[startCI + i]);
                }

                if (activeSet_lastH[coI] != activeSet[coI]) {
                    activeSetChanged = true;
                    activeSet_lastH[coI] = activeSet[coI];
                }
            }
        }
        for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
            animConfig.meshCollisionObjects[coI]->evaluateConstraints(result,
                MMActiveSet[coI], constraintVal);
        }
        if (animConfig.isSelfCollision) {
            int startCI = constraintVal.size();
            SelfCollisionHandler<dim>::evaluateConstraints(result, MMActiveSet.back(), constraintVal);

            if (animConfig.selfFric > 0.0) {
                MMLambda_lastH.back().resize(constraintVal.size() - startCI);
                //TODO: parallelize
                for (int i = 0; i < MMLambda_lastH.back().size(); ++i) {
                    compute_g_b(constraintVal[startCI + i], dHat, MMLambda_lastH.back()[i]);
                    MMLambda_lastH.back()[i] *= -mu_IP * 2.0 * std::sqrt(constraintVal[startCI + i]);
                    if (MMActiveSet.back()[i][3] < -1) {
                        // PP or PE duplication
                        MMLambda_lastH.back()[i] *= -MMActiveSet.back()[i][3];
                    }
                }

                SelfCollisionHandler<dim>::computeDistCoordAndTanBasis(
                    result, MMActiveSet.back(), MMDistCoord.back(), MMTanBasis.back());
                if (MMActiveSet_lastH.back() != MMActiveSet.back()) {
                    activeSetChanged = true;
                    MMActiveSet_lastH.back() = MMActiveSet.back();
                }
            }
        }
        timer_step.stop();

        bool updateDHat = true, updateFricDHat = true;
        if (constraintVal.size()) {
            logFile << "d range: [" << constraintVal.minCoeff() << ", " << constraintVal.maxCoeff() << "]" << std::endl;

            timer_step.start(10);
            Eigen::VectorXd fb(constraintVal.size());
#ifdef USE_TBB
            tbb::parallel_for(0, (int)fb.size(), 1, [&](int i)
#else
            for (int i = 0; i < fb.size(); ++i)
#endif
                {
                    double dualI;
                    compute_g_b(constraintVal[i], dHat, dualI);
                    dualI *= -mu_IP;
                    const double& constraint = constraintVal[i];
                    fb[i] = (dualI + constraint - std::sqrt(dualI * dualI + constraint * constraint));
                }
#ifdef USE_TBB
            );
#endif
            double fbNorm = fb.norm();
            spdlog::info("||fb|| = {:g}", fbNorm);
            timer_step.stop();

            if (constraintVal.maxCoeff() < fricDHatThres && fricDHat < 0.0) {
                // setup friction homotopy
                fricDHat = fricDHat0 * 2.0;
            }

#ifdef USE_DISCRETE_CMS
            if (constraintVal.maxCoeff() < dHatTarget) {
                logFile << "Discrete complimentarity slackness condition satisfied" << std::endl;
                updateDHat = false;
            }
            else if (constraintVal.minCoeff() < dTol) {
                logFile << "Tiny distance fail-safe triggered" << std::endl;
                break;
            }

            if (solveFric) {
                // test on tangent space convergence
                if (fricDHat <= fricDHatTarget) {
                    // compute gradient with updated tangent,
                    // if gradient still below CNTol,
                    // then friction tangent space has converged
                    computeGradient(result, false, gradient, true);
                    computePrecondMtr(result, false, linSysSolver, false, true);
                    if (!linSysSolver->factorize()) {
                        linSysSolver->precondition_diag(gradient, searchDir);
                    }
                    else {
                        Eigen::VectorXd minusG = -gradient;
                        linSysSolver->solve(minusG, searchDir);
                    }

                    if (searchDir.cwiseAbs().maxCoeff() < targetGRes) {
                        updateFricDHat = false;
                    }
                    if (animConfig.fricIterAmt > 0 && fricIterI >= animConfig.fricIterAmt) {
                        updateFricDHat = false;
                    }
                }
            }
            else {
                updateFricDHat = false;
            }

#else // USE_DISCRETE_CMS
            if (fbNorm < fbNormTol || constraintVal.minCoeff() < dTol) {
                break;
            }
#endif // USE_DISCRETE_CMS
        }
        else {
            spdlog::info("no collision in this time step");
            break;
        }

        // if (!updateDHat && !updateFricDHat && !activeSetChanged) {
        if (!updateDHat && !updateFricDHat) {
            break;
        }

        if constexpr (HOMOTOPY_VAR == 0) {
            mu_IP *= 0.5;
            spdlog::info("mu decreased to {:g}", mu_IP);
        }
        else if (HOMOTOPY_VAR == 1) {
            if (updateDHat) {
                dHat *= 0.5;
                if (dHat < dHatTarget) {
                    dHat = dHatTarget;
                }
                computeConstraintSets(result);
                spdlog::info("dHat decreased to {:g}", dHat);
#ifdef ADAPTIVE_MU
                initMu_IP(mu_IP);
#endif
            }

            if (updateFricDHat && fricDHat > 0.0) {
                fricDHat *= 0.5;
                if (fricDHat < fricDHatTarget) {
                    fricDHat = fricDHatTarget;
                }
            }
        }
        else {
            spdlog::error("needs to define HOMOTOPY_VAR to either 0 (mu) or 1 (dHat)");
            exit(-1);
        }
    }

    for (int i = 0; i < MMActiveSet.size(); ++i) {
        n_convCollPairs_sum[i] += MMActiveSet[i].size();
        if (n_convCollPairs_max[i] < MMActiveSet[i].size()) {
            n_convCollPairs_max[i] = MMActiveSet[i].size();
        }
    }
    if (n_convCollPairs_total_max < constraintStartInds.back()) {
        n_convCollPairs_total_max = constraintStartInds.back();
    }

    std::vector<double> sysE;
    std::vector<Eigen::Matrix<double, 1, dim>> sysM, sysL;
    computeSystemEnergy(sysE, sysM, sysL);

    for (const auto& i : sysE) {
        file_sysE << i << " ";
    }
    file_sysE << std::endl;

    for (const auto& i : sysM) {
        file_sysM << i << " ";
    }
    file_sysM << std::endl;

    for (const auto& i : sysL) {
        file_sysL << i << " ";
    }
    file_sysL << std::endl;

    return false;
}

template <int dim>
bool Optimizer<dim>::solveSub_IP(double mu, std::vector<std::vector<int>>& AHat,
    std::vector<std::vector<MMCVID>>& MMAHat)
{
    int iterCap = 10000, k = 0; // 10000 for running comparative schemes that can potentially not converge
    m_projectDBC = true;
    rho_DBC = 0.0;
    double lastMove = animScripter.getCompletedStepSize();
    for (; k < iterCap; ++k) {
        // output constraints for debug
        // Eigen::VectorXd constraintVal;
        // for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
        //     animConfig.meshCollisionObjects[coI]->evaluateConstraints(result, MMActiveSet[coI], constraintVal);
        // }
        // if (animConfig.isSelfCollision) {
        //     SelfCollisionHandler<dim>::evaluateConstraints(result, MMActiveSet.back(), constraintVal);
        // }
        // std::cout << constraintStartInds.back() << " constraints:" << std::endl;
        // int i = 0;
        // for (const auto& coI : MMActiveSet) {
        //     for (const auto& cI : coI) {
        //         std::cout << cI[0] << " " << cI[1] << " " << cI[2] << " " << cI[3] << " " << constraintVal[i++] << std::endl;
        //     }
        //     std::cout << "-----------" << std::endl;
        // }
        // std::cout << "size = " << MMActiveSet.back().size() << std::endl;
        // std::cout << "dHat = " << dHat << std::endl;
        buildConstraintStartIndsWithMM(activeSet, MMActiveSet, constraintStartInds);
        for (int i = 0; i < MMActiveSet.size(); ++i) {
            n_collPairs_sum[i] += MMActiveSet[i].size();
            if (n_collPairs_max[i] < MMActiveSet[i].size()) {
                n_collPairs_max[i] = MMActiveSet[i].size();
            }
        }
        if (n_collPairs_total_max < constraintStartInds.back()) {
            n_collPairs_total_max = constraintStartInds.back();
        }

        // compute gradient
        timer_step.start(12);
        computeGradient(result, false, gradient, m_projectDBC);
        timer_step.stop();
        //TODO: constraint val can be computed and reused for E, g, H

        spdlog::info("# constraint = {:d}", MMActiveSet.back().size());
        spdlog::info("# paraEE = {:d}", paraEEMMCVIDSet.back().size());

        // check convergence
        double gradSqNorm = gradient.squaredNorm();
        double distToOpt_PN = searchDir.cwiseAbs().maxCoeff();
        spdlog::info("mu = {:g}, dHat = {:g}, {:g}, subproblem ||g||^2 = {:g}", mu_IP, dHat, fricDHat, gradSqNorm);
        spdlog::info("distToOpt_PN = {:g}, targetGRes = {:g}", distToOpt_PN, targetGRes);
        bool gradVanish = (distToOpt_PN < targetGRes);
        if (!useGD && k && gradVanish && (animScripter.getCompletedStepSize() > 1.0 - 1.0e-3)
            && (animScripter.getCOCompletedStepSize() > 1.0 - 1.0e-3)) {
            // subproblem converged
            logFile << k << " Newton-type iterations for mu = " << mu_IP << ", dHat = " << dHat << std::endl;
            return false;
        }
        innerIterAmt++;

        computeSearchDir(k, m_projectDBC);

        // largest feasible step size
        timer_step.start(13);
        double alpha = 1.0, slackness_a = 0.9, slackness_m = 0.8;
        energyTerms[0]->filterStepSize(result, searchDir, alpha);
        for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
            animConfig.collisionObjects[coI]->largestFeasibleStepSize(result, searchDir, slackness_a, AHat[coI], alpha);
        }
#if (CFL_FOR_CCD == 0) // use full CCD
#ifdef USE_SH_LFSS
        timer_temp3.start(11);
        sh.build(result, searchDir, alpha, result.avgEdgeLen / 3.0);
        timer_temp3.stop();
#endif
#endif
        // full CCD or partial CCD
        for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
            if (animConfig.exactCCDMethod == ExactCCD::Method::NONE) {
                animConfig.meshCollisionObjects[coI]->largestFeasibleStepSize(result, sh, searchDir, slackness_m, MMActiveSet_CCD[coI], alpha);
            }
            else {
                animConfig.meshCollisionObjects[coI]->largestFeasibleStepSize_exact(result, sh, searchDir, animConfig.exactCCDMethod, MMActiveSet_CCD[coI], alpha);
            }
        }
        std::vector<std::pair<int, int>> newCandidates;
        if (animConfig.isSelfCollision) {
            if (animConfig.exactCCDMethod == ExactCCD::Method::NONE) {
                SelfCollisionHandler<dim>::largestFeasibleStepSize(result, sh, searchDir, slackness_m, MMActiveSet_CCD.back(), newCandidates, alpha);
            }
            else {
                SelfCollisionHandler<dim>::largestFeasibleStepSize_exact(result, sh, searchDir, animConfig.exactCCDMethod, MMActiveSet_CCD.back(), alpha);
            }
        }

#if (CFL_FOR_CCD != 0) // CFL for CCD acceleration
        if (MMActiveSet_CCD.size()) { // if there is mesh-mesh collision (CO or selfCollision)
            Eigen::VectorXd pMag(result.SVI.size());
            for (int i = 0; i < pMag.size(); ++i) {
                pMag[i] = searchDir.template segment<dim>(result.SVI[i] * dim).norm();
            }
            double alpha_CFL = std::sqrt(dHat) / (pMag.maxCoeff() * 2.0);
            spdlog::info("partial CCD {:g}, CFL {:g}", alpha, alpha_CFL);

            if constexpr (CFL_FOR_CCD == 1) {
                alpha = std::min(alpha, alpha_CFL);
                spdlog::info("take the smaller one {:g}", alpha);
            }
            else if (CFL_FOR_CCD == 2) {
                if ((!k && alpha > alpha_CFL) || alpha > 2.0 * alpha_CFL) {
                    // full CCD
#ifdef USE_SH_LFSS
                    timer_temp3.start(11);
                    sh.build(result, searchDir, alpha, result.avgEdgeLen / 3.0);
                    timer_temp3.stop();
#endif
                    for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
                        if (animConfig.exactCCDMethod == ExactCCD::Method::NONE) {
                            animConfig.meshCollisionObjects[coI]->largestFeasibleStepSize_CCD(result, sh, searchDir, slackness_m, alpha);
                        }
                        else {
                            animConfig.meshCollisionObjects[coI]->largestFeasibleStepSize_CCD_exact(result, sh, searchDir, animConfig.exactCCDMethod, alpha);
                        }
                    }
                    if (animConfig.isSelfCollision) {
                        if (animConfig.exactCCDMethod == ExactCCD::Method::NONE) {
                            SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD(result, sh, searchDir, slackness_m, newCandidates, alpha);
                        }
                        else {
                            SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD_exact(result, sh, searchDir, animConfig.exactCCDMethod, alpha);
                        }
                    }

                    if (alpha < alpha_CFL) {
                        alpha = alpha_CFL;
                        spdlog::info("take a full CCD instead but failed, take CFL {:g}", alpha);
                    }
                    else {
                        spdlog::info("take a full CCD instead {:g}", alpha);
                    }
                }
                else {
                    alpha = std::min(alpha, alpha_CFL);
                    spdlog::info("take the smaller one {:g}", alpha);
                }
            }
        }
#endif
        // double largestFeasibleStepSize = alpha;

        if (alpha == 0.0) {
            spdlog::info("CCD gives 0 step size");
            exit(-1);

            alpha = 1.0; // fail-safe, let safe-guard in line search find the stepsize
        }
        timer_step.stop();

        file_iterStats << globalIterNum << " ";

        // Eigen::MatrixXd Vlast = result.V;
        double alpha_feasible = alpha;
        lineSearch(alpha);
        if (alpha_feasible < 1.0e-8) {
            logFile << "tiny feasible step size " << alpha_feasible << " " << alpha << std::endl;
        }
        else if (alpha < 1.0e-8) {
            logFile << "tiny step size after armijo " << alpha_feasible << " " << alpha << std::endl;
        }

#ifdef LINSYSSOLVER_USE_AMGCL
        if (alpha < 1.0e-8) { // only for debugging tiny step size issue
            // output system in AMGCL format
            // Eigen::VectorXd minusG = -gradient;
            // dynamic_cast<AMGCLSolver<Eigen::VectorXi, Eigen::VectorXd>*>(linSysSolver)->write_AMGCL("Stiff", minusG);
            // exit(-1);
#else
        if (false) {
#endif
            // tiny step size issue
            // if (useGD) {
            // MMCVID i = paraEEMMCVIDSet.back()[0];
            // for (int j = 0; j < 4; ++j) {
            //     printf("%.20le %.20le %.20le\n",
            //         Vlast(i[j], 0), Vlast(i[j], 1), Vlast(i[j], 2));
            // }
            // for (int j = 0; j < 4; ++j) {
            //     printf("%.20le %.20le %.20le\n",
            //         searchDir[i[j] * dim], searchDir[i[j] * dim + 1], searchDir[i[j] * dim + 2]);
            // }
            // FILE* out = fopen((outputFolderPath + "E.txt").c_str(), "w");
            // for (int i = -1000; i < 1000; ++i) {
            //     stepForward(Vlast, result, i / 1000.0 * largestFeasibleStepSize);
            //     double E;
            //     computeConstraintSets(result);
            //     computeEnergyVal(result, true, E);
            //     fprintf(out, "%.20le %.20le %.20le %.20le %.20le %.20le\n", i / 1000.0 * largestFeasibleStepSize,
            //         energies[0], energies[1], energies[2], energies[3], E);
            //     if (i == -1) {
            //         for (int j = -999; j < 0; ++j) {
            //             stepForward(Vlast, result, j / 1000000.0 * largestFeasibleStepSize);
            //             double E;
            //             computeConstraintSets(result);
            //             computeEnergyVal(result, true, E);
            //             fprintf(out, "%.20le %.20le %.20le %.20le %.20le %.20le\n", j / 1000000.0 * largestFeasibleStepSize,
            //                 energies[0], energies[1], energies[2], energies[3], E);
            //             if (j > -15) {
            //                 std::cout << "j=" << j << std::endl;
            //                 for (const auto& i : MMActiveSet.back()) {
            //                     std::cout << i[0] << " " << i[1] << " " << i[2] << " " << i[3] << std::endl;
            //                 }
            //             }
            //         }
            //     }
            //     else if (i == 0) {
            //         for (int j = 1; j < 1000; ++j) {
            //             stepForward(Vlast, result, j / 1000000.0 * largestFeasibleStepSize);
            //             double E;
            //             computeConstraintSets(result);
            //             computeEnergyVal(result, true, E);
            //             fprintf(out, "%.20le %.20le %.20le %.20le %.20le %.20le\n", j / 1000000.0 * largestFeasibleStepSize,
            //                 energies[0], energies[1], energies[2], energies[3], E);
            //             if (j < 15) {
            //                 std::cout << "j=" << j << std::endl;
            //                 for (const auto& i : MMActiveSet.back()) {
            //                     std::cout << i[0] << " " << i[1] << " " << i[2] << " " << i[3] << std::endl;
            //                 }
            //             }
            //         }
            //     }
            // }
            // fclose(out);
            //     exit(0);
            // }

            saveStatus("tinyStepSize");
            useGD = true;
            // // mu_IP *= 1.1;
            // if constexpr (dim == 3) {
            //     spdlog::info("check EE pairs:");
            //     for (const auto cI : MMActiveSet.back()) {
            //         if (cI[0] >= 0) { // EE
            //             const Eigen::Matrix<double, 1, dim> e01 = result.V.row(cI[1]) - result.V.row(cI[0]);
            //             const Eigen::Matrix<double, 1, dim> e23 = result.V.row(cI[3]) - result.V.row(cI[2]);
            //             double sine = std::sqrt(e01.cross(e23).squaredNorm() / (e01.squaredNorm() * e23.squaredNorm()));
            //             if (sine < 0.02) {
            //                 spdlog::info("{:d}-{:d} {:d}-{:d} sine = {:g}", cI[0], cI[1], cI[2], cI[3], sine);
            //             }
            //         }
            //     }
            //     spdlog::info("cosine p g = {:g}", searchDir.dot(gradient) / std::sqrt(searchDir.squaredNorm() * gradient.squaredNorm()));
            // }
        }
        else {
            useGD = false;
        }

#ifdef CHECK_ET_INTERSECTION
        bool intersected = false;
        for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
            intersected |= !animConfig.meshCollisionObjects[coI]->checkEdgeTriIntersection(result, sh);
        }
        if (animConfig.isSelfCollision) {
            intersected |= !SelfCollisionHandler<dim>::checkEdgeTriIntersection(result, sh);
        }
        if (intersected) {
            getchar();
        }
#endif

        postLineSearch(alpha);

        file_iterStats << constraintStartInds.back() << std::endl;

        if (m_projectDBC) {
            if (animScripter.getCompletedStepSize() < 1.0 - 1.0e-3) {
                // setup KKT solving option
                m_projectDBC = false;
                rho_DBC = 1.0e6;
                spdlog::info("setup penalty solve");
            }
        }
        else {
            double completedStepSize = animScripter.computeCompletedStepSize(result);
            spdlog::info("completed step size of MBC = {:g}", completedStepSize);
            if (completedStepSize > 1.0 - 1.0e-3) {
                m_projectDBC = true;
                spdlog::info("penalty solve finished");
            }
            else {
                if (completedStepSize < lastMove && rho_DBC < 1.0e8) {
                    rho_DBC *= 2.0;
                    spdlog::info("set rho_DBC to {:g}", rho_DBC);
                }
                else {
                    bool safeToPull = (searchDir.cwiseAbs().maxCoeff() < CN_MBC);
                    if (safeToPull) {
                        if (completedStepSize < 0.99 && rho_DBC < 1.0e8) {
                            rho_DBC *= 2.0;
                            spdlog::info("set rho_DBC to {:g}", rho_DBC);
                        }
                        else {
                            animScripter.updateLambda(result, rho_DBC);
                            //TODO: use second order dual update
                            spdlog::info("MBC lambda updated with rho_DBC = {:g}", rho_DBC);
                        }
                    }
                }
            }
        }
    }

    if (k >= iterCap) {
        logFile << "iteration cap reached for IP subproblem solve!!!" << std::endl;
        exit(0); // for stopping comparative schemes that can potentially not converge
    }

    return true;
} // namespace IPC

template <int dim>
void Optimizer<dim>::upperBoundMu(double& mu)
{
    double H_b;
    compute_H_b(1.0e-16 * bboxDiagSize2, dHat, H_b);
    double muMax = 1.0e13 * result.avgNodeMass(dim) / (4.0e-16 * bboxDiagSize2 * H_b);
    if (mu > muMax) {
        mu = muMax;
        logFile << "upper bounded mu to " << muMax << " at dHat = " << dHat << std::endl;
    }
}

template <int dim>
void Optimizer<dim>::suggestMu(double& mu)
{
    double H_b;
    compute_H_b(1.0e-16 * bboxDiagSize2, dHat, H_b);
    mu = 1.0e11 * result.avgNodeMass(dim) / (4.0e-16 * bboxDiagSize2 * H_b);
}

template <int dim>
void Optimizer<dim>::initMu_IP(double& mu)
{
    //TODO: optimize implementation and pipeline
    buildConstraintStartIndsWithMM(activeSet, MMActiveSet, constraintStartInds);
    if (constraintStartInds.back()) {
        Eigen::VectorXd g_E;
        solveIP = false;
        computeGradient(result, true, g_E);
        solveIP = true;

        Eigen::VectorXd constraintVal, g_c;
        g_c.setZero(g_E.size());
        for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
            animConfig.collisionObjects[coI]->evaluateConstraints(result, activeSet[coI], constraintVal);
            int startCI = constraintStartInds[coI];
            for (int cI = startCI; cI < constraintStartInds[coI + 1]; ++cI) {
                compute_g_b(constraintVal[cI], dHat, constraintVal[cI]);
            }
            animConfig.collisionObjects[coI]->leftMultiplyConstraintJacobianT(result, activeSet[coI],
                constraintVal.segment(startCI, activeSet[coI].size()), g_c);
        }
        for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
            animConfig.meshCollisionObjects[coI]->evaluateConstraints(result, MMActiveSet[coI], constraintVal);
            int startCI = constraintStartInds[coI + animConfig.collisionObjects.size()];
            for (int cI = startCI; cI < constraintStartInds[coI + 1 + animConfig.collisionObjects.size()]; ++cI) {
                compute_g_b(constraintVal[cI], dHat, constraintVal[cI]);
            }
            animConfig.meshCollisionObjects[coI]->leftMultiplyConstraintJacobianT(result, MMActiveSet[coI],
                constraintVal.segment(startCI, MMActiveSet[coI].size()), g_c);
        }
        if (animConfig.isSelfCollision) {
            SelfCollisionHandler<dim>::evaluateConstraints(result, MMActiveSet.back(), constraintVal);
            int startCI = constraintStartInds[animConfig.meshCollisionObjects.size() + animConfig.collisionObjects.size()];
            for (int cI = startCI; cI < constraintVal.size(); ++cI) {
                compute_g_b(constraintVal[cI], dHat, constraintVal[cI]);
            }
            SelfCollisionHandler<dim>::leftMultiplyConstraintJacobianT(result, MMActiveSet.back(),
                constraintVal.segment(startCI, MMActiveSet.back().size()), g_c);
        }
        for (const auto& fixedVI : result.fixedVert) {
            g_c.segment<dim>(fixedVI * dim).setZero();
        }

        // balance current gradient at constrained DOF
        double mu_trial = -g_c.dot(g_E) / g_c.squaredNorm();
        if (mu_trial > 0.0) {
            mu = mu_trial;
        }
        suggestMu(mu_trial);
        if (mu < mu_trial) {
            mu = mu_trial;
        }
        upperBoundMu(mu);

        // minimize approximated next search direction
        //            Eigen::VectorXd HInvGc, HInvGE;
        //            linSysSolver->solve(g_c, HInvGc);
        //            linSysSolver->solve(g_E, HInvGE);
        //            mu = HInvGE.dot(HInvGc) / HInvGc.squaredNorm();
        //            if(mu < 0.0) {
        //                mu = g_c.dot(g_E) / g_c.squaredNorm();
        //            }

        // minimize approximated next search direction at constrained DOF
        //            Eigen::VectorXd HInvGc;
        //            linSysSolver->solve(g_c, HInvGc);
        //            mu = g_E.dot(HInvGc) / g_c.dot(HInvGc);
        //            if(mu < 0.0) {
        //                mu = g_c.dot(g_E) / g_c.squaredNorm();
        //            }

        spdlog::info("mu initialized to {:g}", mu);
    }
    else {
        spdlog::info("no constraint detected, start with default mu");
    }
}

template <int dim>
void Optimizer<dim>::initSubProb_IP(void)
{
    closeConstraintID.resize(0);
    closeMConstraintID.resize(0);
    closeConstraintVal.resize(0);
    closeMConstraintVal.resize(0);
}
template <int dim>
void Optimizer<dim>::computeSearchDir(int k, bool projectDBC)
{
    // compute matrix
    computePrecondMtr(result, false, linSysSolver, false, projectDBC);

    timer_step.start(3);
    Eigen::VectorXd minusG = -gradient;
    if (useGD || !linSysSolver->factorize()) {
        if (!useGD) {
            spdlog::warn("matrix not positive definite");
            logFile << "matrix not positive definite" << std::endl;
        }

        // std::cout << "writing it to disk for debug..." << std::endl;
        // IglUtils::writeSparseMatrixToFile(outputFolderPath + "H_IP", linSysSolver, true);
        // std::cout << "writing complete, terminate program" << std::endl;
        // exit(0);

        spdlog::warn("use diagonal preconditioned -gradient direction");
        logFile << "use diagonal preconditioned -gradient direction" << std::endl;
        timer_step.start(4);
        linSysSolver->precondition_diag(minusG, searchDir);
        // searchDir = minusG;
        timer_step.stop();
    }
    else {
        // solve for searchDir
        timer_step.start(4);
        linSysSolver->solve(minusG, searchDir);
        timer_step.stop();
    }
}
template <int dim>
void Optimizer<dim>::postLineSearch(double alpha)
{
#ifdef ADAPTIVE_MU
    if (mu_IP == 0.0) {
        initMu_IP(mu_IP);
    }
    else {
        //TODO: avoid recomputation of constraint functions
        bool updateMu = false;
        for (int i = 0; i < closeConstraintID.size(); ++i) {
            double d;
            animConfig.collisionObjects[closeConstraintID[i].first]->evaluateConstraint(result,
                closeConstraintID[i].second, d);
            if (d <= closeConstraintVal[i]) {
                updateMu = true;
                break;
            }
        }
        if (!updateMu) {
            for (int i = 0; i < closeMConstraintID.size(); ++i) {
                double d;
                if (closeMConstraintID[i].first < animConfig.meshCollisionObjects.size()) {
                    animConfig.meshCollisionObjects[closeMConstraintID[i].first]->evaluateConstraint(result,
                        closeMConstraintID[i].second, d);
                }
                else {
                    SelfCollisionHandler<dim>::evaluateConstraint(result, closeMConstraintID[i].second, d);
                }
                if (d <= closeMConstraintVal[i]) {
                    updateMu = true;
                    break;
                }
            }
        }
        if (updateMu) {
            mu_IP *= 2.0;
            upperBoundMu(mu_IP);
        }
#endif

#ifdef ADAPTIVE_MU
        closeConstraintID.resize(0);
        closeMConstraintID.resize(0);
        closeConstraintVal.resize(0);
        closeMConstraintVal.resize(0);
        Eigen::VectorXd constraintVal;
        for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
            int constraintValIndStart = constraintVal.size();
            animConfig.collisionObjects[coI]->evaluateConstraints(result,
                activeSet[coI], constraintVal);

            for (int i = 0; i < activeSet[coI].size(); ++i) {
                if (constraintVal[constraintValIndStart + i] < dTol) {
                    closeConstraintID.emplace_back(coI, activeSet[coI][i]);
                    closeConstraintVal.emplace_back(constraintVal[constraintValIndStart + i]);
                }
            }
        }
        for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
            int constraintValIndStart = constraintVal.size();
            animConfig.meshCollisionObjects[coI]->evaluateConstraints(result,
                MMActiveSet[coI], constraintVal);

            for (int i = 0; i < MMActiveSet[coI].size(); ++i) {
                // std::cout << MMActiveSet[coI][i][0] << " " << MMActiveSet[coI][i][1] << " " << MMActiveSet[coI][i][2] << " " << MMActiveSet[coI][i][3] << ", d=" << constraintVal[constraintValIndStart + i] << std::endl;
                if (constraintVal[constraintValIndStart + i] < dTol) {
                    closeMConstraintID.emplace_back(coI, MMActiveSet[coI][i]);
                    closeMConstraintVal.emplace_back(constraintVal[constraintValIndStart + i]);
                }
            }
        }
        if (animConfig.isSelfCollision) {
            int constraintValIndStart = constraintVal.size();
            SelfCollisionHandler<dim>::evaluateConstraints(result, MMActiveSet.back(), constraintVal);

            for (int i = 0; i < MMActiveSet.back().size(); ++i) {
                // std::cout << MMActiveSet.back()[i][0] << " " << MMActiveSet.back()[i][1] << " " << MMActiveSet.back()[i][2] << " " << MMActiveSet.back()[i][3] << ", d=" << constraintVal[constraintValIndStart + i] << std::endl;
                if (constraintVal[constraintValIndStart + i] < dTol) {
                    closeMConstraintID.emplace_back(MMActiveSet.size() - 1, MMActiveSet.back()[i]);
                    closeMConstraintVal.emplace_back(constraintVal[constraintValIndStart + i]);
                }
            }
        }
        if (constraintVal.size()) {
            spdlog::info("min distance^2 = {:g}", constraintVal.minCoeff());
        }
    }
#endif
}

template <int dim>
void Optimizer<dim>::computeConstraintSets(const Mesh<dim>& data, bool rehash)
{
    timer_step.start(14);
#ifndef CCD_FILTERED_CS
#ifdef USE_SH_CCS
    if (rehash) {
        timer_temp3.start(7);
        sh.build(data, data.avgEdgeLen / 3.0);
        timer_temp3.stop();
    }
#endif
#endif
    for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
        animConfig.collisionObjects[coI]->computeConstraintSet(data, dHat, activeSet[coI]);
    }
    for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
        animConfig.meshCollisionObjects[coI]->computeConstraintSet(data, sh, dHat, MMActiveSet[coI], paraEEMMCVIDSet[coI], paraEEeIeJSet[coI], CFL_FOR_CCD, MMActiveSet_CCD[coI]);
    }
    if (animConfig.isSelfCollision) {
        SelfCollisionHandler<dim>::computeConstraintSet(data, sh, dHat, MMActiveSet.back(), paraEEMMCVIDSet.back(), paraEEeIeJSet.back(), CFL_FOR_CCD, MMActiveSet_CCD.back());
    }
    timer_step.stop();
}

template <int dim>
void Optimizer<dim>::buildConstraintStartInds(const std::vector<std::vector<int>>& activeSet,
    std::vector<int>& constraintStartInds)
{
    constraintStartInds.resize(1);
    constraintStartInds[0] = 0;
    for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
        constraintStartInds.emplace_back(constraintStartInds.back() + activeSet[coI].size());
    }
}
template <int dim>
void Optimizer<dim>::buildConstraintStartIndsWithMM(const std::vector<std::vector<int>>& activeSet,
    const std::vector<std::vector<MMCVID>>& MMActiveSet,
    std::vector<int>& constraintStartInds)
{
    constraintStartInds.resize(1);
    constraintStartInds[0] = 0;
    for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
        constraintStartInds.emplace_back(constraintStartInds.back() + activeSet[coI].size());
    }
    for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
        constraintStartInds.emplace_back(constraintStartInds.back() + MMActiveSet[coI].size());
    }
    if (animConfig.isSelfCollision) {
        constraintStartInds.emplace_back(constraintStartInds.back() + MMActiveSet.back().size());
    }
}

template <int dim>
bool Optimizer<dim>::solve_oneStep(void)
{
    if (!solveWithQP) {
        // for the changing hessian
        // Update matrix entries
        computePrecondMtr(result, false, linSysSolver);

        try {
            // Numerically factorizing Hessian/Proxy matrix
            if (!mute) { timer_step.start(3); }
            if (!solveWithSQP) {
                linSysSolver->factorize();
                // output factorization and exit
                // linSysSolver->outputFactorization("/Users/mincli/Desktop/test/IPC/output/L");
                // std::cout << "matrix written" << std::endl;
                // exit(0);
            }
            if (!mute) { timer_step.stop(); }
        }
        catch (std::exception e) {
            IglUtils::writeSparseMatrixToFile(outputFolderPath + "mtr_numFacFail",
                linSysSolver, true);
            spdlog::error("numerical factorization failed, matrix written into {}mtr_numFacFail", outputFolderPath);
            exit(-1);
        }
    }

    // Back solve
    if (!mute) { timer_step.start(4); }
    if (solveWithQP || solveWithSQP) {
#ifdef EXIT_UPON_QP_FAIL
        bool prevQPsuccess = solveQPSuccess;
#endif
        solveQPSuccess = solveQP(
            result, animConfig.collisionObjects, activeSet,
            linSysSolver, P_QP, elemPtr_P_QP, gradient,
            OSQPSolver,
#ifdef USE_GUROBI
            gurobiQPSolver,
#endif
            constraintStartInds, searchDir, dual_QP,
            animConfig.qpSolverType);

#ifdef USE_GUROBI
        if (animConfig.qpSolverType == QP_GUROBI) {
            if (solveQPSuccess) {
                spdlog::debug(
                    "timestep={:d} Gurobi_status={:d} description=\"{:s}\"",
                    globalIterNum, gurobiQPSolver.status(), gurobiQPSolver.statusDescription());
            }
            else {
                spdlog::error(
                    "timestep={:d} Gurobi_status={:d} description=\"{:s}\"",
                    globalIterNum, gurobiQPSolver.status(), gurobiQPSolver.statusDescription());
            }
        }
#endif

        if (!solveQPSuccess) {
#ifdef EXIT_UPON_QP_FAIL
            if (!prevQPsuccess) {
                spdlog::error(
                    "timestep={:d} msg=\"unable to solve QP for two successive iterations; exiting\"",
                    globalIterNum);
                exit(3);
            }
#endif
            spdlog::error("timestep={:d} msg=\"unable to solve QP; resetting active set\"", globalIterNum);

            searchDir.setZero(gradient.size());
            for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
                activeSet_next[coI].clear();
            }
            for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
                MMActiveSet_next[coI].clear();
            }
            if (animConfig.isSelfCollision) {
                MMActiveSet_next.back().clear();
            }
        }
    }
    else {
        Eigen::VectorXd minusG = -gradient;
        linSysSolver->solve(minusG, searchDir);
    }
    if (!mute) { timer_step.stop(); }

    double stepSize;
    bool stopped;
    if (solveWithQP || solveWithSQP) {
        // stopped = lineSearch(stepSize);
        stepSize = 1.0;
        Eigen::MatrixXd resultV0 = result.V;
        if (animConfig.energyType == ET_NH) {
            // Line search to prevent element inversions
            stopped = true; // Fail until a valid stepSize is found
            energyTerms[0]->filterStepSize(result, searchDir, stepSize);
            while (stepSize > 0) {
                stepSize /= 2.0;
                stepForward(resultV0, result, stepSize);
                if (result.checkInversion(true)) {
                    // No inversions found
                    stopped = false;
                    break;
                }
            }
        }
        else {
            stopped = false;
            stepForward(resultV0, result, stepSize);
        }
    }
    else {
        stopped = lineSearch(stepSize); // Also takes the step
    }
    if (stopped) {
        // IglUtils::writeSparseMatrixToFile(outputFolderPath + "precondMtr_stopped_" + std::to_string(globalIterNum), precondMtr);
        // logFile << "descent step stopped at overallIter" << globalIterNum << " for no prominent energy decrease." << std::endl;
    }
    if (!solveWithQP) { // Do not update the gradient for a single QP
        computeGradient(result, false, gradient);
    }
    return stopped;
}

template <int dim>
bool isIntersected(
    const Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    const Eigen::MatrixXd& V0,
    const Config& p_animConfig)
{
#ifdef NO_CCD_FAILSAFE
    return false;
#endif

    for (const auto& obj : p_animConfig.collisionObjects) {
        if (obj->isIntersected(mesh, V0)) {
            logFile << "mesh intersects with analytical collision object" << std::endl;
            return true;
        }
    }

    for (int coI = 0; coI < p_animConfig.meshCollisionObjects.size(); ++coI) {
        if (!p_animConfig.meshCollisionObjects[coI]->checkEdgeTriIntersectionIfAny(mesh, sh)) {
            logFile << "mesh intersects with mesh collision object" << std::endl;
            return true;
        }
    }

    if (p_animConfig.isSelfCollision) {
        if (!SelfCollisionHandler<dim>::checkEdgeTriIntersectionIfAny(mesh, sh)) {
            logFile << "mesh intersects with itself" << std::endl;
            return true;
        }
    }

    return false;
}

template <int dim>
bool Optimizer<dim>::lineSearch(double& stepSize,
    double armijoParam,
    double lowerBound)
{
    timer_step.start(5);
    std::stringstream msg;

    bool outputLineSearch = false;
    if (outputLineSearch) {
        result.saveAsMesh(outputFolderPath + ((dim == 2) ? "lineSearchBase.obj" : "lineSearchBase.msh"));
    }

    bool stopped = false;

    if (!solveIP) {
        initStepSize(result, stepSize);
    }
    else {
        timer_step.start(9);
        computeEnergyVal(result, false, lastEnergyVal); //TODO: only for updated constraints set and mu
        timer_step.start(5);
    }
    msg << "E_last = " << lastEnergyVal << " stepSize: " << stepSize << " -> ";

    // const double m = searchDir.dot(gradient);
    // const double c1m = 1.0e-4 * m;
    double c1m = 0.0;
    if (armijoParam > 0.0) {
        c1m = armijoParam * searchDir.dot(gradient);
    }
    Eigen::MatrixXd resultV0 = result.V;
    // Mesh<dim> temp = result; // TEST

    if (outputLineSearch) {
        Eigen::VectorXd energyValPerElem;
        energyTerms[0]->getEnergyValPerElemBySVD(result, false, svd, F,
            energyValPerElem, false);
        IglUtils::writeVectorToFile(outputFolderPath + "E0.txt", energyValPerElem);
    }

#ifdef OUTPUT_CCD_FAIL
    bool output = false;
    if (solveIP && animConfig.exactCCDMethod != ExactCCD::Method::NONE) {
        result.saveSurfaceMesh(outputFolderPath + "before" + std::to_string(numOfCCDFail) + ".obj");
        output = true;
    }
#endif
    stepForward(resultV0, result, stepSize);
    if (energyTerms[0]->getNeedElemInvSafeGuard()) {
        while (!result.checkInversion(true)) {
            logFile << "element inversion detected during line search, backtrack!" << std::endl;
            stepSize /= 2.0;
            stepForward(resultV0, result, stepSize);
        }
        msg << stepSize << "(safeGuard) -> ";
    }
    bool rehash = true;
    if (solveIP) {
        sh.build(result, result.avgEdgeLen / 3.0);
        while (isIntersected(result, sh, resultV0, animConfig)) {
            logFile << "intersection detected during line search, backtrack!" << std::endl;
#ifdef OUTPUT_CCD_FAIL
            if (output) {
                result.saveSurfaceMesh(outputFolderPath + "after" + std::to_string(numOfCCDFail) + ".obj");
                ++numOfCCDFail;
                output = false;
            }
#endif
            stepSize /= 2.0;
            stepForward(resultV0, result, stepSize);
            sh.build(result, result.avgEdgeLen / 3.0);
        }
        msg << stepSize << "(safeGuard_IP) -> ";
        rehash = false;
    }
    double testingE;
    //        Eigen::VectorXd testingG;
    timer_step.stop();
    if (solveIP) {
        computeConstraintSets(result, rehash);
    }
    timer_step.start(9);
    computeEnergyVal(result, 2, testingE);
    timer_step.start(5);
    //        computeGradient(testingData, testingG);

    if (outputLineSearch) {
        Eigen::VectorXd energyValPerElem;
        energyTerms[0]->getEnergyValPerElemBySVD(result, false, svd, F,
            energyValPerElem, false);
        IglUtils::writeVectorToFile(outputFolderPath + "E.txt", energyValPerElem);
    }

#define ARMIJO_RULE
#ifdef ARMIJO_RULE
    //        while((testingE > lastEnergyVal + stepSize * c1m) ||
    //              (searchDir.dot(testingG) < c2m)) // Wolfe condition
    // FILE* out = fopen((outputFolderPath + "Els.txt").c_str(), "w");
    double LFStepSize = stepSize;
    while ((testingE > lastEnergyVal + stepSize * c1m) && // Armijo condition
        (stepSize > lowerBound)) {
        // fprintf(out, "%.9le %.9le\n", stepSize, testingE);
        if (stepSize == 1.0) {
            // can try cubic interpolation here
            stepSize /= 2.0;
        }
        else {
            stepSize /= 2.0;
        }

        ++numOfLineSearch;
        if (stepSize == 0.0) {
            stopped = true;
            if (!mute) {
                logFile << "testingE" << globalIterNum << " " << testingE << " > " << lastEnergyVal << " " << stepSize * c1m << std::endl;
                //                    logFile << "testingG" << globalIterNum << " " << searchDir.dot(testingG) << " < " << c2m << std::endl;
            }
            break;
        }

        stepForward(resultV0, result, stepSize);
        // while (isIntersected(result, resultV0, animConfig)) {
        //     stepSize /= 2.0;
        //     stepForward(resultV0, result, stepSize);
        // }
        msg << stepSize << " -> ";

        timer_step.stop();
        if (solveIP) {
            computeConstraintSets(result);
        }
        timer_step.start(9);
        computeEnergyVal(result, 2, testingE);
        timer_step.start(5);
        //            computeGradient(testingData, testingG);
    }

    if (stepSize < LFStepSize && solveIP) {
        bool needRecomputeCS = false;
        while (isIntersected(result, sh, resultV0, animConfig)) {
            logFile << "intersection detected after line search, backtrack!" << std::endl;
            stepSize /= 2.0;
            stepForward(resultV0, result, stepSize);
            sh.build(result, result.avgEdgeLen / 3.0);
            needRecomputeCS = true;
        }
        if (needRecomputeCS) {
            computeConstraintSets(result, false);
        }
    }

#ifdef CHECK_RATIONAL_CCD_GLOBAL
    Eigen::MatrixXd V_bk = result.V;

    result.V = resultV0;
    sh.build(result, searchDir, stepSize, result.avgEdgeLen / 3.0);
    double stepSize_rational = stepSize;
    SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD_exact(result, sh, searchDir, ExactCCD::Method::RATIONAL_ROOT_PARITY, stepSize_rational);
    if (stepSize_rational != stepSize) {
        std::cout << "rational CCD detects interpenetration but inexact didn't" << std::endl;
    }

    result.V = V_bk;
#endif

    // fclose(out);
#endif
    // // further search
    // double stepSize_fs = stepSize / 2.0;
    // stepForward(resultV0, result, stepSize_fs);
    // double E_fs;
    // timer_step.stop();
    // computeConstraintSets(result);
    // timer_step.start(9);
    // computeEnergyVal(result, 2, E_fs);
    // timer_step.start(5);
    // while (E_fs < testingE) {
    //     testingE = E_fs;
    //     stepSize_fs /= 2.0;
    //     stepForward(resultV0, result, stepSize_fs);

    //     timer_step.stop();
    //     computeConstraintSets(result);
    //     timer_step.start(9);
    //     computeEnergyVal(result, 2, E_fs);
    //     timer_step.start(5);
    // }
    // stepForward(resultV0, result, stepSize_fs * 2.0);
    // timer_step.stop();
    // computeConstraintSets(result);
    // timer_step.start(9);
    // computeEnergyVal(result, 2, testingE);
    // timer_step.start(5);

    msg << stepSize << "(armijo) ";

    // if (energyTerms[0]->getNeedElemInvSafeGuard()) {
    //     if (!result.checkInversion(true)) {
    //         logFile << "element inversion detected after line search, backtrack!" << std::endl;
    //         stepSize /= 2.0;
    //         stepForward(resultV0, result, stepSize);
    //     }
    // }
    // if (solveIP) {
    //     while (isIntersected(result, sh, resultV0, animConfig)) {
    //         logFile << "intersection detected after line search, backtrack!" << std::endl;
    //         stepSize /= 2.0;
    //         stepForward(resultV0, result, stepSize);
    //         sh.build(result, result.avgEdgeLen / 3.0);
    //     }
    // }

    // while((!result.checkInversion()) ||
    //       ((scaffolding) && (!scaffold.airMesh.checkInversion())))
    // {
    //     assert(0 && "element inversion after armijo shouldn't happen!");
    //
    //     stepSize /= 2.0;
    //     if(stepSize == 0.0) {
    //         assert(0 && "line search failed!");
    //         stopped = true;
    //         break;
    //     }
    //
    //     stepForward(resultV0, scaffoldV0, result, scaffold, stepSize);
    //     computeEnergyVal(result, scaffold, testingE);
    // }
    //
    // // projection method for collision handling
    // for(int vI = 0; vI < result.V.rows(); vI++) {
    //     if(result.V(vI, 1) < 0.0) {
    //         result.V(vI, 1) = 0.0;
    //     }
    // }
    // computeEnergyVal(result, scaffold, true, testingE);

    lastEnergyVal = testingE;

    if (!mute) {
        msg << stepSize << " (E = " << testingE << ")";
        spdlog::info("{:s}", msg.str());
        spdlog::info("stepLen = {:g}", (stepSize * searchDir).squaredNorm());
    }
    file_iterStats << stepSize << " ";

    if (outputLineSearch) {
        IglUtils::writeVectorToFile(outputFolderPath + "searchDir.txt",
            searchDir);
        exit(0);
    }

    timer_step.stop();

    return stopped;
}

template <int dim>
void Optimizer<dim>::stepForward(const Eigen::MatrixXd& dataV0,
    Mesh<dim>& data,
    double stepSize) const
{
    assert(dataV0.rows() == data.V.rows());
    assert(data.V.rows() * dim == searchDir.size());
    assert(data.V.rows() == result.V.rows());

#ifdef USE_TBB
    tbb::parallel_for(0, (int)data.V.rows(), 1, [&](int vI)
#else
    for (int vI = 0; vI < data.V.rows(); vI++)
#endif
        {
            data.V.row(vI) = dataV0.row(vI) + stepSize * searchDir.segment<dim>(vI * dim).transpose();
        }
#ifdef USE_TBB
    );
#endif
}

template <int dim>
void Optimizer<dim>::updateTargetGRes(void)
{
    targetGRes = std::sqrt(relGL2Tol * bboxDiagSize2 * dtSq);
}

template <int dim>
void Optimizer<dim>::getFaceFieldForVis(Eigen::VectorXd& field)
{
    //        field = Eigen::VectorXd::Zero(result.F.rows());
    field = result.u;
}

template <int dim>
void Optimizer<dim>::initStepSize(const Mesh<dim>& data, double& stepSize)
{
    stepSize = 1.0;

    for (int eI = 0; eI < energyTerms.size(); eI++) {
        energyTerms[eI]->filterStepSize(data, searchDir, stepSize);
    }
}

template <int dim>
void Optimizer<dim>::saveStatus(const std::string& appendStr)
{
    FILE* out = fopen((outputFolderPath + "status" + std::to_string(globalIterNum) + appendStr).c_str(), "w");
    assert(out);

    fprintf(out, "timestep %d\n", globalIterNum);

    fprintf(out, "\nposition %ld %ld\n", result.V.rows(), result.V.cols());
    for (int vI = 0; vI < result.V.rows(); ++vI) {
        fprintf(out, "%le %le", result.V(vI, 0),
            result.V(vI, 1));
        if constexpr (dim == 3) {
            fprintf(out, " %le\n", result.V(vI, 2));
        }
        else {
            fprintf(out, "\n");
        }
    }

    fprintf(out, "\nvelocity %ld\n", velocity.size());
    for (int velI = 0; velI < velocity.size(); ++velI) {
        fprintf(out, "%le\n", velocity[velI]);
    }

    fprintf(out, "\ndx_Elastic %ld %d\n", dx_Elastic.rows(), dim);
    for (int velI = 0; velI < dx_Elastic.rows(); ++velI) {
        fprintf(out, "%le %le", dx_Elastic(velI, 0),
            dx_Elastic(velI, 1));
        if constexpr (dim == 3) {
            fprintf(out, " %le\n", dx_Elastic(velI, 2));
        }
        else {
            fprintf(out, "\n");
        }
    }

    fclose(out);

    // surface mesh, including codimensional surface CO
#ifdef USE_TBB
    tbb::parallel_for(0, (int)V_surf.rows(), 1, [&](int vI)
#else
    for (int vI = 0; vI < V_surf.rows(); ++vI)
#endif
        {
            V_surf.row(vI) = result.V.row(surfIndToTet[vI]);
        }
#ifdef USE_TBB
    );
#endif

    igl::writeOBJ(outputFolderPath + std::to_string(globalIterNum) + ".obj", V_surf, F_surf);

    // codimensional segment CO
    if (result.CE.rows()) {
#ifdef USE_TBB
        tbb::parallel_for(0, (int)V_CE.rows(), 1, [&](int vI)
#else
        for (int vI = 0; vI < V_CE.rows(); ++vI)
#endif
            {
                V_CE.row(vI) = result.V.row(CEIndToTet[vI]);
            }
#ifdef USE_TBB
        );
#endif
        IglUtils::writeSEG(outputFolderPath + std::to_string(globalIterNum) + ".seg",
            V_CE, F_CE);
    }

    // codimensional point CO
    bool hasPointCO = false;
    for (const auto& coDimI : result.componentCoDim) {
        if (coDimI == 0) {
            hasPointCO = true;
            break;
        }
    }
    if (hasPointCO) {
        FILE* out = fopen((outputFolderPath + "pt" + std::to_string(globalIterNum) + ".obj").c_str(), "w");
        for (int compI = 0; compI < result.componentCoDim.size(); ++compI) {
            if (result.componentCoDim[compI] == 0) {
                for (int vI = result.componentNodeRange[compI];
                     vI < result.componentNodeRange[compI + 1]; ++vI) {
                    fprintf(out, "v %le %le %le\n", result.V(vI, 0), result.V(vI, 1), result.V(vI, 2));
                }
            }
        }
        fclose(out);
    }
}

template <int dim>
void Optimizer<dim>::outputCollStats(std::ostream& out)
{
    for (int i = 0; i < n_collPairs_sum.size(); ++i) {
        out << "avg # collision pairs / iter = " << n_collPairs_sum[i] / (innerIterAmt ? innerIterAmt : 1) << std::endl;
        out << "max # collision pairs / iter = " << n_collPairs_max[i] << std::endl;
        out << "avg # collision pairs / end of ts = " << n_convCollPairs_sum[i] / globalIterNum << std::endl;
        out << "max # collision pairs / end of ts = " << n_convCollPairs_max[i] << std::endl;
        out << std::endl;
    }

    out << "Total:" << std::endl;
    out << "avg # collision pairs / iter = " << n_collPairs_sum.sum() / (innerIterAmt ? innerIterAmt : 1) << std::endl;
    out << "max # collision pairs / iter = " << n_collPairs_total_max << std::endl;
    out << "avg # collision pairs / end of ts = " << n_convCollPairs_sum.sum() / globalIterNum << std::endl;
    out << "max # collision pairs / end of ts = " << n_convCollPairs_total_max << std::endl;
    out << std::endl;
}

template <int dim>
void Optimizer<dim>::computeEnergyVal(const Mesh<dim>& data, int redoSVD, double& energyVal)
{
    //TODO: write inertia and augmented Lagrangian term into energyTerms
    //        if(!mute) { timer_step.start(0); }
    switch (animConfig.timeIntegrationType) {
    case TIT_BE: {
        energyTerms[0]->computeEnergyVal(data, redoSVD, svd, F,
            dtSq * energyParams[0], energyVal_ET[0]);
        energyVal = energyVal_ET[0];
        for (int eI = 1; eI < energyTerms.size(); eI++) {
            energyTerms[eI]->computeEnergyVal(data, redoSVD, svd, F,
                dtSq * energyParams[eI], energyVal_ET[eI]);
            energyVal += energyVal_ET[eI];
        }
        break;
    }

    case TIT_NM: {
        energyTerms[0]->computeEnergyVal(data, redoSVD, svd, F, dtSq * beta_NM * energyParams[0], energyVal_ET[0]);
        energyVal = energyVal_ET[0];
        for (int eI = 1; eI < energyTerms.size(); eI++) {
            energyTerms[eI]->computeEnergyVal(data, redoSVD, svd, F, dtSq * beta_NM * energyParams[eI], energyVal_ET[eI]);
            energyVal += energyVal_ET[eI];
        }
        break;
    }
    }

    Eigen::VectorXd energyVals(data.V.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)data.V.rows(), 1, [&](int vI)
#else
    for (int vI = 0; vI < data.V.rows(); vI++)
#endif
        {
            energyVals[vI] = ((data.V.row(vI) - xTilta.row(vI)).squaredNorm() * data.massMatrix.coeff(vI, vI) / 2.0);
        }
#ifdef USE_TBB
    );
#endif
    energyVal += energyVals.sum();

    if (animScripter.isNBCActive()) {
        for (const auto& NMI : data.NeumannBC) {
            if (!data.isFixedVert[NMI.first]) {
                energyVal -= dtSq * data.massMatrix.coeff(NMI.first, NMI.first) * data.V.row(NMI.first).dot(NMI.second);
            }
        }
    }

    if (solveIP) {
        Eigen::VectorXd constraintVals, bVals;
        for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
            int startCI = constraintVals.size();
            animConfig.collisionObjects[coI]->evaluateConstraints(data, activeSet[coI], constraintVals);
            bVals.conservativeResize(constraintVals.size());
            for (int cI = startCI; cI < constraintVals.size(); ++cI) {
                if (constraintVals[cI] <= 0.0) {
                    spdlog::error("constraintVal = {:g} when evaluating constraints", constraintVals[cI]);
                    exit(0);
                }
                else {
                    compute_b(constraintVals[cI], dHat, bVals[cI]);
                }
            }
        }
        for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
            int startCI = constraintVals.size();
            animConfig.meshCollisionObjects[coI]->evaluateConstraints(data, MMActiveSet[coI], constraintVals);
            bVals.conservativeResize(constraintVals.size());
            for (int cI = startCI; cI < constraintVals.size(); ++cI) {
                if (constraintVals[cI] <= 0.0) {
                    spdlog::error("constraintVal = {:g} when evaluating constraints", constraintVals[cI]);
                    exit(0);
                }
                else {
                    compute_b(constraintVals[cI], dHat, bVals[cI]);
                    int duplication = MMActiveSet[coI][cI - startCI][3];
                    if (duplication < -1) {
                        // PP or PE or EP, handle duplication
                        bVals[cI] *= -duplication;
                    }
                }
            }

            animConfig.meshCollisionObjects[coI]->augmentParaEEEnergy(data, paraEEMMCVIDSet[coI], paraEEeIeJSet[coI],
                constraintVals, bVals, dHat, 1.0);
        }
        if (animConfig.isSelfCollision) {
            int startCI = constraintVals.size();
            SelfCollisionHandler<dim>::evaluateConstraints(data, MMActiveSet.back(), constraintVals);
            bVals.conservativeResize(constraintVals.size());
            //TODO: parallelize
            for (int cI = startCI; cI < constraintVals.size(); ++cI) {
                if (constraintVals[cI] <= 0.0) {
                    spdlog::error("constraintVal = {:g} when evaluating constraints {:d} {:d} {:d} {:d}",
                        constraintVals[cI], MMActiveSet.back()[cI][0], MMActiveSet.back()[cI][1],
                        MMActiveSet.back()[cI][2], MMActiveSet.back()[cI][3]);
                    std::cout << data.V.row(MMActiveSet.back()[cI][0] >= 0 ? MMActiveSet.back()[cI][0] : (-MMActiveSet.back()[cI][0] - 1)) << std::endl;
                    std::cout << data.V.row(MMActiveSet.back()[cI][1] >= 0 ? MMActiveSet.back()[cI][1] : (-MMActiveSet.back()[cI][1] - 1)) << std::endl;
                    std::cout << data.V.row(MMActiveSet.back()[cI][2] >= 0 ? MMActiveSet.back()[cI][2] : (-MMActiveSet.back()[cI][2] - 1)) << std::endl;
                    std::cout << data.V.row(MMActiveSet.back()[cI][3] >= 0 ? MMActiveSet.back()[cI][3] : (-MMActiveSet.back()[cI][3] - 1)) << std::endl;
                    data.saveSurfaceMesh("debug.obj");
                    exit(0);
                }
                else {
                    compute_b(constraintVals[cI], dHat, bVals[cI]);
                    int duplication = MMActiveSet.back()[cI - startCI][3];
                    if (duplication < -1) {
                        // PP or PE, handle duplication
                        bVals[cI] *= -duplication;
                    }
                }
            }

            startCI = constraintVals.size();
            SelfCollisionHandler<dim>::evaluateConstraints(data, paraEEMMCVIDSet.back(), constraintVals);
            bVals.conservativeResize(constraintVals.size());
            for (int cI = startCI; cI < constraintVals.size(); ++cI) {
                if (constraintVals[cI] <= 0.0) {
                    spdlog::error("constraintVal = {:g} when evaluating constraints {:d} {:d} {:d} {:d}",
                        constraintVals[cI], MMActiveSet.back()[cI][0], MMActiveSet.back()[cI][1],
                        MMActiveSet.back()[cI][2], MMActiveSet.back()[cI][3]);
                    std::cout << data.V.row(MMActiveSet.back()[cI][0] >= 0 ? MMActiveSet.back()[cI][0] : (-MMActiveSet.back()[cI][0] - 1)) << std::endl;
                    std::cout << data.V.row(MMActiveSet.back()[cI][1] >= 0 ? MMActiveSet.back()[cI][1] : (-MMActiveSet.back()[cI][1] - 1)) << std::endl;
                    std::cout << data.V.row(MMActiveSet.back()[cI][2] >= 0 ? MMActiveSet.back()[cI][2] : (-MMActiveSet.back()[cI][2] - 1)) << std::endl;
                    std::cout << data.V.row(MMActiveSet.back()[cI][3] >= 0 ? MMActiveSet.back()[cI][3] : (-MMActiveSet.back()[cI][3] - 1)) << std::endl;
                    data.saveSurfaceMesh("debug.obj");
                    exit(0);
                }
                else {
                    const MMCVID& MMCVIDI = paraEEMMCVIDSet.back()[cI - startCI];
                    double eps_x, e;
                    if (MMCVIDI[3] >= 0) {
                        // EE
                        compute_eps_x(data, MMCVIDI[0], MMCVIDI[1], MMCVIDI[2], MMCVIDI[3], eps_x);
                        compute_e(data.V.row(MMCVIDI[0]), data.V.row(MMCVIDI[1]), data.V.row(MMCVIDI[2]), data.V.row(MMCVIDI[3]), eps_x, e);
                    }
                    else {
                        // PP or PE
                        const std::pair<int, int>& eIeJ = paraEEeIeJSet.back()[cI - startCI];
                        const std::pair<int, int>& eI = data.SFEdges[eIeJ.first];
                        const std::pair<int, int>& eJ = data.SFEdges[eIeJ.second];
                        compute_eps_x(data, eI.first, eI.second, eJ.first, eJ.second, eps_x);
                        compute_e(data.V.row(eI.first), data.V.row(eI.second), data.V.row(eJ.first), data.V.row(eJ.second), eps_x, e);
                    }
                    compute_b(constraintVals[cI], dHat, bVals[cI]);
                    bVals[cI] *= e;
                }
            }
        }
        energyVal += mu_IP * bVals.sum();
        // std::cout << ", E_c=" << mu_IP * bVals.sum();
        // energies[1] = mu_IP * bVals.sum();

        for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
            // friction
            if (activeSet_lastH[coI].size() && fricDHat > 0.0 && animConfig.collisionObjects[coI]->friction > 0.0) {
                double Ef;
                animConfig.collisionObjects[coI]->computeFrictionEnergy(data.V, result.V_prev, activeSet_lastH[coI],
                    lambda_lastH[coI], Ef, fricDHat, 1.0);
                energyVal += Ef;
                // std::cout << ", E_f" << coI << "=" << Ef;
                // energies[2] = Ef;
            }
        }
        if (animConfig.isSelfCollision) {
            if (MMActiveSet_lastH.back().size() && fricDHat > 0.0 && animConfig.selfFric > 0.0) {
                double Ef;
                SelfCollisionHandler<dim>::computeFrictionEnergy(data.V, result.V_prev, MMActiveSet_lastH.back(),
                    MMLambda_lastH.back(), MMDistCoord.back(), MMTanBasis.back(), Ef, fricDHat, animConfig.selfFric);
                energyVal += Ef;
                // std::cout << ", E_fs=" << Ef;
                // energies[3] = Ef;
            }
        }
    }
    // std::cout << std::endl;

    if (animConfig.dampingStiff > 0.0) {
        Eigen::VectorXd displacement(data.V.rows() * dim), Adx;
#ifdef USE_TBB
        tbb::parallel_for(0, (int)data.V.rows(), 1, [&](int vI)
#else
        for (int vI = 0; vI < data.V.rows(); vI++)
#endif
            {
                displacement.segment<dim>(vI * dim) = (data.V.row(vI) - result.V_prev.row(vI)).transpose();
            }
#ifdef USE_TBB
        );
#endif
        for (const auto& fixedVI : data.fixedVert) {
            displacement.segment<dim>(fixedVI * dim).setZero();
        }

        dampingMtr->multiply(displacement, Adx);
        energyVal += 0.5 * displacement.dot(Adx);
    }

    if (rho_DBC) {
        animScripter.augmentMDBCEnergy(data, energyVal, rho_DBC);
    }

    //        if(!mute) { timer_step.stop(); }
}
template <int dim>
void Optimizer<dim>::computeGradient(const Mesh<dim>& data,
    bool redoSVD, Eigen::VectorXd& gradient, bool projectDBC)
{
    //        if(!mute) { timer_step.start(0); }

    switch (animConfig.timeIntegrationType) {
    case TIT_BE: {
        energyTerms[0]->computeGradient(data, redoSVD, svd, F,
            dtSq * energyParams[0], gradient_ET[0], projectDBC);
        gradient = gradient_ET[0];
        for (int eI = 1; eI < energyTerms.size(); eI++) {
            energyTerms[eI]->computeGradient(data, redoSVD, svd, F,
                dtSq * energyParams[eI], gradient_ET[eI], projectDBC);
            gradient += gradient_ET[eI];
        }
        break;
    }

    case TIT_NM: {
        energyTerms[0]->computeGradient(data, redoSVD, svd, F, dtSq * beta_NM * energyParams[0], gradient_ET[0], projectDBC);
        gradient = gradient_ET[0];
        for (int eI = 1; eI < energyTerms.size(); eI++) {
            energyTerms[eI]->computeGradient(data, redoSVD, svd, F, dtSq * beta_NM * energyParams[eI], gradient_ET[eI], projectDBC);
            gradient += gradient_ET[eI];
        }
        break;
    }
    }

#ifdef USE_TBB
    tbb::parallel_for(0, (int)data.V.rows(), 1, [&](int vI)
#else
    for (int vI = 0; vI < data.V.rows(); vI++)
#endif
        {
            if (!data.isFixedVert[vI] || !projectDBC) {
                gradient.segment<dim>(vI * dim) += (data.massMatrix.coeff(vI, vI) * (data.V.row(vI) - xTilta.row(vI)).transpose());
            }
        }
#ifdef USE_TBB
    );
#endif

    if (animScripter.isNBCActive()) {
        for (const auto& NMI : data.NeumannBC) {
            if (!data.isFixedVert[NMI.first]) {
                gradient.template segment<dim>(NMI.first * dim) -= dtSq * data.massMatrix.coeff(NMI.first, NMI.first) * NMI.second.transpose();
            }
        }
    }

    if (solveIP) {
        Eigen::VectorXd constraintVal;
        for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
            int startCI = constraintVal.size();
            animConfig.collisionObjects[coI]->evaluateConstraints(data, activeSet[coI], constraintVal);
            for (int cI = startCI; cI < constraintVal.size(); ++cI) {
                compute_g_b(constraintVal[cI], dHat, constraintVal[cI]);
            }
            animConfig.collisionObjects[coI]->leftMultiplyConstraintJacobianT(data, activeSet[coI],
                constraintVal.segment(startCI, activeSet[coI].size()), gradient, mu_IP);

            // friction
            if (activeSet_lastH[coI].size() && fricDHat > 0.0 && animConfig.collisionObjects[coI]->friction > 0.0) {
                animConfig.collisionObjects[coI]->augmentFrictionGradient(data.V, result.V_prev, activeSet_lastH[coI],
                    lambda_lastH[coI], gradient, fricDHat, 1.0);
            }
        }
        for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
            int startCI = constraintVal.size();
            animConfig.meshCollisionObjects[coI]->evaluateConstraints(data, MMActiveSet[coI], constraintVal);
            for (int cI = startCI; cI < constraintVal.size(); ++cI) {
                compute_g_b(constraintVal[cI], dHat, constraintVal[cI]);
            }
            animConfig.meshCollisionObjects[coI]->leftMultiplyConstraintJacobianT(data, MMActiveSet[coI],
                constraintVal.segment(startCI, MMActiveSet[coI].size()), gradient, mu_IP);

            animConfig.meshCollisionObjects[coI]->augmentParaEEGradient(data,
                paraEEMMCVIDSet[coI], paraEEeIeJSet[coI], gradient, dHat, mu_IP);
        }
        if (animConfig.isSelfCollision) {
            int startCI = constraintVal.size();
            SelfCollisionHandler<dim>::evaluateConstraints(data, MMActiveSet.back(), constraintVal);
            for (int cI = startCI; cI < constraintVal.size(); ++cI) {
                compute_g_b(constraintVal[cI], dHat, constraintVal[cI]);
            }
            SelfCollisionHandler<dim>::leftMultiplyConstraintJacobianT(data, MMActiveSet.back(),
                constraintVal.segment(startCI, MMActiveSet.back().size()), gradient, mu_IP);

            SelfCollisionHandler<dim>::augmentParaEEGradient(data,
                paraEEMMCVIDSet.back(), paraEEeIeJSet.back(), gradient, dHat, mu_IP);

            if (MMActiveSet_lastH.back().size() && fricDHat > 0.0 && animConfig.selfFric > 0.0) {
                SelfCollisionHandler<dim>::augmentFrictionGradient(data.V, result.V_prev, MMActiveSet_lastH.back(),
                    MMLambda_lastH.back(), MMDistCoord.back(), MMTanBasis.back(), gradient, fricDHat, animConfig.selfFric);
            }
        }
        if (projectDBC) {
            for (const auto& fixedVI : data.fixedVert) {
                gradient.segment<dim>(fixedVI * dim).setZero();
            }
        }
    }

    if (animConfig.dampingStiff > 0.0) {
        Eigen::VectorXd displacement(data.V.rows() * dim), Adx;
#ifdef USE_TBB
        tbb::parallel_for(0, (int)data.V.rows(), 1, [&](int vI)
#else
        for (int vI = 0; vI < data.V.rows(); vI++)
#endif
            {
                displacement.segment<dim>(vI * dim) = (data.V.row(vI) - result.V_prev.row(vI)).transpose();
            }
#ifdef USE_TBB
        );
#endif
        if (projectDBC) {
            for (const auto& fixedVI : data.fixedVert) {
                displacement.segment<dim>(fixedVI * dim).setZero();
            }
        }

        dampingMtr->multiply(displacement, Adx);
        gradient += Adx;
    }

    if (!projectDBC && rho_DBC) {
        animScripter.augmentMDBCGradient(data, gradient, rho_DBC);
    }

    //        if(!mute) { timer_step.stop(); }
}
template <int dim>
void Optimizer<dim>::computePrecondMtr(const Mesh<dim>& data,
    bool redoSVD,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* p_linSysSolver,
    bool updateDamping, bool projectDBC)
{
    if (!mute) { timer_step.start(0); }

    if (solveIP) {
        // figure out vNeighbor_IP in this iteration
        timer_step.start(1);
        if (animConfig.isSelfCollision && (MMActiveSet.back().size() + paraEEeIeJSet.back().size() + MMActiveSet_lastH.back().size())) {
            // there is extra connectivity from self-contact
            timer_mt.start(8);
            std::vector<std::set<int>> vNeighbor_IP_new = data.vNeighbor;
            SelfCollisionHandler<dim>::augmentConnectivity(data, MMActiveSet.back(), vNeighbor_IP_new);
            SelfCollisionHandler<dim>::augmentConnectivity(data, paraEEMMCVIDSet.back(), paraEEeIeJSet.back(), vNeighbor_IP_new);
            if (MMActiveSet_lastH.back().size() && fricDHat > 0.0 && animConfig.selfFric > 0.0) {
                SelfCollisionHandler<dim>::augmentConnectivity(data, MMActiveSet_lastH.back(), vNeighbor_IP_new);
            }
            timer_mt.stop();

            if (vNeighbor_IP_new != vNeighbor_IP) {
                timer_mt.start(8);
                vNeighbor_IP = vNeighbor_IP_new;
                timer_mt.start(9);
                p_linSysSolver->set_pattern(vNeighbor_IP, data.fixedVert);
                timer_mt.stop();

                timer_step.start(2);
                p_linSysSolver->analyze_pattern();
                timer_step.start(1);
            }
        }
        else if (vNeighbor_IP != data.vNeighbor) {
            // no extra connectivity in this iteration but there is in the last iteration
            timer_mt.start(8);
            vNeighbor_IP = data.vNeighbor;
            timer_mt.start(9);
            p_linSysSolver->set_pattern(vNeighbor_IP, data.fixedVert);
            timer_mt.stop();

            timer_step.start(2);
            p_linSysSolver->analyze_pattern();
            timer_step.start(1);
        }
        timer_step.start(0);
    }

    timer_mt.start(10);
    if (updateDamping && animConfig.dampingStiff) {
        computeDampingMtr(data, redoSVD, dampingMtr, projectDBC);

        switch (animConfig.timeIntegrationType) {
        case TIT_BE: {
            p_linSysSolver->setCoeff(dampingMtr,
                1.0 + dtSq * dt / animConfig.dampingStiff);
            break;
        }

        case TIT_NM: {
            p_linSysSolver->setCoeff(dampingMtr,
                1.0 + dtSq * beta_NM * dt / animConfig.dampingStiff);
            break;
        }
        }
    }
    else {
        p_linSysSolver->setZero();
        switch (animConfig.timeIntegrationType) {
        case TIT_BE: {
            for (int eI = 0; eI < energyTerms.size(); eI++) {
                energyTerms[eI]->computeHessian(data, redoSVD, svd, F,
                    energyParams[eI] * dtSq,
                    p_linSysSolver, true, projectDBC);
            }
            break;
        }

        case TIT_NM: {
            for (int eI = 0; eI < energyTerms.size(); eI++) {
                energyTerms[eI]->computeHessian(data, redoSVD, svd, F, energyParams[eI] * dtSq * beta_NM, p_linSysSolver, true, projectDBC);
            }
            break;
        }
        }
    }

    timer_mt.start(11);
#ifdef USE_TBB
    tbb::parallel_for(0, (int)data.V.rows(), 1, [&](int vI)
#else
    for (int vI = 0; vI < data.V.rows(); vI++)
#endif
        {
            if (!data.isFixedVert[vI] || !projectDBC) {
                double massI = data.massMatrix.coeff(vI, vI);
                int ind0 = vI * dim;
                int ind1 = ind0 + 1;
                p_linSysSolver->addCoeff(ind0, ind0, massI);
                p_linSysSolver->addCoeff(ind1, ind1, massI);
                if constexpr (dim == 3) {
                    int ind2 = ind0 + 2;
                    p_linSysSolver->addCoeff(ind2, ind2, massI);
                }
            }
            else {
                // for Dirichlet boundary condition
                int ind0 = vI * dim;
                int ind1 = ind0 + 1;
                p_linSysSolver->setCoeff(ind0, ind0, 1.0);
                p_linSysSolver->setCoeff(ind1, ind1, 1.0);
                if constexpr (dim == 3) {
                    int ind2 = ind0 + 2;
                    p_linSysSolver->setCoeff(ind2, ind2, 1.0);
                }
            }
        }
#ifdef USE_TBB
    );
#endif
    timer_mt.stop();

    if (solveIP) {
        timer_mt.start(12);
        for (int coI = 0; coI < animConfig.collisionObjects.size(); ++coI) {
            animConfig.collisionObjects[coI]->augmentIPHessian(data,
                activeSet[coI], p_linSysSolver, dHat, mu_IP, projectDBC);

            // friction
            if (activeSet_lastH[coI].size() && fricDHat > 0.0 && animConfig.collisionObjects[coI]->friction > 0.0) {
                animConfig.collisionObjects[coI]->augmentFrictionHessian(data,
                    result.V_prev, activeSet_lastH[coI], lambda_lastH[coI],
                    p_linSysSolver, fricDHat, 1.0, projectDBC);
            }
        }
        timer_mt.start(13);
        for (int coI = 0; coI < animConfig.meshCollisionObjects.size(); ++coI) {
            animConfig.meshCollisionObjects[coI]->augmentIPHessian(data, MMActiveSet[coI], p_linSysSolver, dHat, mu_IP, projectDBC);

            animConfig.meshCollisionObjects[coI]->augmentParaEEHessian(data, paraEEMMCVIDSet[coI], paraEEeIeJSet[coI],
                p_linSysSolver, dHat, mu_IP, projectDBC);
        }
        timer_mt.start(14);
        if (animConfig.isSelfCollision) {
            SelfCollisionHandler<dim>::augmentIPHessian(data, MMActiveSet.back(), p_linSysSolver, dHat, mu_IP, projectDBC);

            SelfCollisionHandler<dim>::augmentParaEEHessian(data, paraEEMMCVIDSet.back(), paraEEeIeJSet.back(),
                p_linSysSolver, dHat, mu_IP, projectDBC);

            if (MMActiveSet_lastH.back().size() && fricDHat > 0.0 && animConfig.selfFric > 0.0) {
                SelfCollisionHandler<dim>::augmentFrictionHessian(data, result.V_prev, MMActiveSet_lastH.back(),
                    MMLambda_lastH.back(), MMDistCoord.back(), MMTanBasis.back(),
                    p_linSysSolver, fricDHat, animConfig.selfFric, projectDBC);
            }
        }
        timer_mt.stop();
    }

    if (animConfig.dampingStiff && (!updateDamping)) {
        p_linSysSolver->addCoeff(dampingMtr, 1.0);
    }

    if (!projectDBC && rho_DBC) {
        animScripter.augmentMDBCHessian(data, p_linSysSolver, rho_DBC);
    }

    if (!mute) { timer_step.stop(); }
    // output matrix and exit
    //        IglUtils::writeSparseMatrixToFile("/Users/mincli/Desktop/IPC/output/A", p_linSysSolver, true);
    //        std::cout << "matrix written" << std::endl;
    //        exit(0);
}

template <int dim>
void Optimizer<dim>::computeDampingMtr(const Mesh<dim>& data,
    bool redoSVD,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* p_dampingMtr,
    bool projectDBC)
{
    p_dampingMtr->setZero();
    for (int eI = 0; eI < energyTerms.size(); eI++) {
        energyTerms[eI]->computeHessian(data, redoSVD, svd, F,
            energyParams[eI] * animConfig.dampingStiff / dt,
            p_dampingMtr, true, projectDBC);
    }
}

template <int dim>
void Optimizer<dim>::setupDampingMtr(const Mesh<dim>& data,
    bool redoSVD,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* p_dampingMtr)
{
    dampingMtr->set_pattern(data.vNeighbor, data.fixedVert);
    computeDampingMtr(data, redoSVD, p_dampingMtr);
}

template <int dim>
void Optimizer<dim>::computeSystemEnergy(std::vector<double>& sysE,
    std::vector<Eigen::Matrix<double, 1, dim>>& sysM,
    std::vector<Eigen::Matrix<double, 1, dim>>& sysL)
{
    Eigen::VectorXd energyValPerElem;
    energyTerms[0]->getEnergyValPerElemBySVD(result, true, svd, F, energyValPerElem, false);

    sysE.resize(compVAccSize.size());
    sysM.resize(compVAccSize.size());
    sysL.resize(compVAccSize.size());
    for (int compI = 0; compI < compVAccSize.size(); ++compI) {
        sysE[compI] = 0.0;
        sysM[compI].setZero();
        sysL[compI].setZero();

        for (int fI = (compI ? compFAccSize[compI - 1] : 0); fI < compFAccSize[compI]; ++fI) {
            sysE[compI] += energyValPerElem[fI];
        }
        for (int vI = (compI ? compVAccSize[compI - 1] : 0); vI < compVAccSize[compI]; ++vI) {
            sysE[compI] += result.massMatrix.coeff(vI, vI) * ((result.V.row(vI) - result.V_prev.row(vI)).squaredNorm() / dtSq / 2.0 - gravity.dot(result.V.row(vI).transpose()));

            Eigen::Matrix<double, 1, dim> p = result.massMatrix.coeff(vI, vI) / dt * (result.V.row(vI) - result.V_prev.row(vI));
            sysM[compI] += p;

            if constexpr (dim == 3) {
                sysL[compI] += Eigen::Matrix<double, 1, dim>(result.V.row(vI)).cross(p);
            }
            else {
                sysL[compI][0] += result.V(vI, 0) * p[1] - result.V(vI, 1) * p[0];
            }
        }
    }
}

template <int dim>
void Optimizer<dim>::checkGradient(void)
{
    spdlog::info("checking energy gradient computation...");

    double energyVal0;
    computeConstraintSets(result);
    computeEnergyVal(result, true, energyVal0);

    const double h = 1.0e-6 * igl::avg_edge_length(result.V, result.F);
    Mesh<dim> perturbed = result;
    Eigen::VectorXd gradient_finiteDiff;
    gradient_finiteDiff.resize(result.V.rows() * dim);
    for (int vI = 0; vI < result.V.rows(); vI++) {
        for (int dimI = 0; dimI < dim; dimI++) {
            perturbed.V = result.V;
            perturbed.V(vI, dimI) += h;
            double energyVal_perturbed;
            computeConstraintSets(perturbed);
            computeEnergyVal(perturbed, true, energyVal_perturbed);
            gradient_finiteDiff[vI * dim + dimI] = (energyVal_perturbed - energyVal0) / h;
        }

        if (((vI + 1) % 100) == 0) {
            spdlog::info("{:d}/{:d} vertices computed", vI + 1, result.V.rows());
        }
    }
    for (const auto fixedVI : result.fixedVert) {
        gradient_finiteDiff.segment<dim>(dim * fixedVI).setZero();
    }

    Eigen::VectorXd gradient_symbolic;
    computeConstraintSets(result);
    computeGradient(result, true, gradient_symbolic);

    Eigen::VectorXd difVec = gradient_symbolic - gradient_finiteDiff;
    const double dif_L2 = difVec.norm();
    const double relErr = dif_L2 / gradient_finiteDiff.norm();

    spdlog::info("L2 dist = {:g}, relErr = {:g}", dif_L2, relErr);

    logFile << "check gradient:" << std::endl;
    logFile << "g_symbolic =\n"
            << gradient_symbolic << std::endl;
    logFile << "g_finiteDiff = \n"
            << gradient_finiteDiff << std::endl;
}

template <int dim>
void Optimizer<dim>::checkHessian(void)
{
    //TODO: needs to turn off SPD projection
    spdlog::info("checking hessian computation...");

    Eigen::VectorXd gradient0;
    computeConstraintSets(result);
    computeGradient(result, true, gradient0);
    const double h = 1.0e-6 * igl::avg_edge_length(result.V, result.F);
    Mesh<dim> perturbed = result;
    Eigen::SparseMatrix<double> hessian_finiteDiff;
    hessian_finiteDiff.resize(result.V.rows() * dim, result.V.rows() * dim);
    for (int vI = 0; vI < result.V.rows(); vI++) {
        if (result.fixedVert.find(vI) != result.fixedVert.end()) {
            hessian_finiteDiff.insert(vI * dim, vI * dim) = 1.0;
            hessian_finiteDiff.insert(vI * dim + 1, vI * dim + 1) = 1.0;
            if constexpr (dim == 3) {
                hessian_finiteDiff.insert(vI * dim + 2, vI * dim + 2) = 1.0;
            }
            continue;
        }

        for (int dimI = 0; dimI < dim; dimI++) {
            perturbed.V = result.V;
            perturbed.V(vI, dimI) += h;
            Eigen::VectorXd gradient_perturbed;
            // computeConstraintSets(perturbed); // the transition is only C1 continuous
            computeGradient(perturbed, true, gradient_perturbed);
            Eigen::VectorXd hessian_colI = (gradient_perturbed - gradient0) / h;
            int colI = vI * dim + dimI;
            for (int rowI = 0; rowI < result.V.rows() * dim; rowI++) {
                if ((result.fixedVert.find(rowI / dim) == result.fixedVert.end()) && (hessian_colI[rowI] != 0.0)) {
                    hessian_finiteDiff.insert(rowI, colI) = hessian_colI[rowI];
                }
            }
        }

        if (((vI + 1) % 100) == 0) {
            spdlog::info("{:d}/{:d} vertices computed", vI + 1, result.V.rows());
        }
    }

    Eigen::SparseMatrix<double> hessian_symbolicPK;
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver;
#ifdef LINSYSSOLVER_USE_CHOLMOD
    linSysSolver = new CHOLMODSolver<Eigen::VectorXi, Eigen::VectorXd>();
#elif defined(LINSYSSOLVER_USE_AMGCL)
    linSysSolver = new AMGCLSolver<Eigen::VectorXi, Eigen::VectorXd>();
#else
    linSysSolver = new EigenLibSolver<Eigen::VectorXi, Eigen::VectorXd>();
#endif
    linSysSolver->set_pattern(result.vNeighbor, result.fixedVert);
    computeConstraintSets(result);
    computePrecondMtr(result, true, linSysSolver);
    linSysSolver->getCoeffMtr(hessian_symbolicPK);

    Eigen::SparseMatrix<double> difMtrPK = hessian_symbolicPK - hessian_finiteDiff;
    const double difPK_L2 = difMtrPK.norm();
    const double relErrPK = difPK_L2 / hessian_finiteDiff.norm();
    spdlog::info("PK L2 dist = {:g}, relErr = {:g}", difPK_L2, relErrPK);
    IglUtils::writeSparseMatrixToFile(outputFolderPath + "H_symbolicPK", hessian_symbolicPK, true);

    IglUtils::writeSparseMatrixToFile(outputFolderPath + "H_finiteDiff", hessian_finiteDiff, true);
}

template class Optimizer<DIM>;
} // namespace IPC
