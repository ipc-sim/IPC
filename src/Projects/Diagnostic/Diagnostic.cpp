//
//  Diagnostic.hpp
//  IPC
//
//  Created by Minchen Li on 1/31/18.
//

#ifndef Diagnostic_hpp
#define Diagnostic_hpp

#include "Utils/Types.hpp"
#include <cmath>
#include "Utils/BarrierFunctions.hpp"
#include "CollisionObject/MeshCollisionUtils.hpp"
#include "NeoHookeanEnergy.hpp"
#include "FixedCoRotEnergy.hpp"

#include "MeshCollisionUtils.hpp"
#include "FrictionUtils.hpp"

#include "IglUtils.hpp"
#include "Mesh.hpp"
#include "GIF.hpp"
#include "OSQPWrapper.h"
#include "CTCD.h"

#include "Optimizer.hpp"

#include "LinSysSolver.hpp"
#include "AMGCLSolver.hpp"

// amgcl
#ifdef USE_AMGCL
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#endif

#include "Timer.hpp"

std::ofstream logFile;
std::string outputFolderPath = "output/";

#include <igl/readOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>

#include <cstdio>

namespace IPC {
class Diagnostic {
public:
    static void run(int argc, char* argv[])
    {
        if (argc > 2) {
            int diagMode = 0;
            diagMode = std::stoi(argv[2]);
            switch (diagMode) {

            case 11: { // a unit test for computing dE/dsigma
                std::vector<Energy<DIM>*> e;
                e.emplace_back(new FixedCoRotEnergy<DIM>);
                e.emplace_back(new NeoHookeanEnergy<DIM>);
                for (const auto eI : e) {
                    eI->unitTest_dE_div_dsigma();
                    eI->unitTest_d2E_div_dsigma2();
                    eI->unitTest_BLeftCoef();
                    eI->unitTest_dE_div_dF();
                    eI->unitTest_dP_div_dF();
                    delete eI;
                }
                break;
            }

            case 14: {
                std::cout << "PP:" << std::endl;
                derivTest_PP();
                std::cout << "PE:" << std::endl;
                derivTest_PE();
                std::cout << "PT:" << std::endl;
                derivTest_PT();
                std::cout << "EE:" << std::endl;
                derivTest_EE();
                break;
            }

            case 15: {
                checkDType();
                break;
            }

            case 16: {
                checkEdgeEdgeD();
                break;
            }

            case 17: { // EE CCD unit test
                std::cout << "v = ";
                Eigen::Vector3d v[8];
                for (int i = 0; i < 8; ++i) {
                    std::cin >> v[i][0] >> v[i][1] >> v[i][2];
                }

                std::cout << "etaRatio = ";
                double etaRatio;
                std::cin >> etaRatio;

                double d_sqrt;
                computeEdgeEdgeD(v[0], v[1], v[2], v[3], d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlpha = 1.0;
                if (CTCD::edgeEdgeCTCD(v[0], v[1], v[2], v[3],
                        v[4], v[5], v[6], v[7],
                        etaRatio * d_sqrt,
                        largestAlpha)) {
                    std::cout << "hit" << std::endl;
                }
                std::cout << "largestAlpha = " << largestAlpha << std::endl;

                for (int i = 0; i < 4; ++i) {
                    std::cout << v[i].transpose() << std::endl;
                }
                for (int i = 4; i < 8; ++i) {
                    std::cout << (v[i - 4] + (v[i] - v[i - 4]) * largestAlpha).transpose() << std::endl;
                }

                break;
            }

            case 18: { // PT CCD unit test
                std::cout << "v = ";
                Eigen::Vector3d v[8];
                for (int i = 0; i < 8; ++i) {
                    std::cin >> v[i][0] >> v[i][1] >> v[i][2];
                }

                std::cout << "etaRatio = ";
                double etaRatio;
                std::cin >> etaRatio;

                double d_sqrt;
                computePointTriD(v[0], v[1], v[2], v[3], d_sqrt);
                d_sqrt = std::sqrt(d_sqrt);

                double largestAlpha = 1.0;
                if (CTCD::vertexFaceCTCD(v[0], v[1], v[2], v[3],
                        v[4], v[5], v[6], v[7],
                        etaRatio * d_sqrt,
                        largestAlpha)) {
                    std::cout << "hit" << std::endl;
                }
                std::cout << "largestAlpha = " << largestAlpha << std::endl;

                for (int i = 0; i < 4; ++i) {
                    std::cout << v[i].transpose() << std::endl;
                }
                for (int i = 4; i < 8; ++i) {
                    std::cout << (v[i - 4] + (v[i] - v[i - 4]) * largestAlpha).transpose() << std::endl;
                }

                break;
            }

            case 19: { // check elements from status file
                if (argc <= 4) {
                    std::cout << "please input rest mesh file path and status file path " << std::endl;
                    exit(0);
                }

                std::ifstream in(argv[3]);
                if (!in.is_open()) {
                    std::cout << "rest mesh file path not found" << std::endl;
                    exit(0);
                }
                in.close();

                Eigen::MatrixXd V0;
                Eigen::MatrixXi F, SF;
                IPC::IglUtils::readTetMesh(argv[3], V0, F, SF);

                in.open(argv[4]);
                if (!in.is_open()) {
                    std::cout << "status file path not found" << std::endl;
                    exit(0);
                }

                std::string line;
                int timeStepNum = 0;
                Eigen::MatrixXd V;
                Eigen::VectorXd velocity;
                while (std::getline(in, line)) {
                    std::stringstream ss(line);
                    std::string token;
                    ss >> token;
                    if (token == "timestep") {
                        ss >> timeStepNum;
                    }
                    else if (token == "position") {
                        std::cout << "read restart position" << std::endl;
                        int posRows, dim_in;
                        ss >> posRows >> dim_in;
                        V.resize(posRows, dim_in);
                        for (int vI = 0; vI < posRows; ++vI) {
                            in >> V(vI, 0) >> V(vI, 1);
                            if constexpr (DIM == 3) {
                                in >> V(vI, 2);
                            }
                        }
                    }
                    else if (token == "velocity") {
                        std::cout << "read restart velocity" << std::endl;
                        int velDim;
                        ss >> velDim;
                        assert(velDim == V.rows() * DIM);
                        velocity.conservativeResize(velDim);
                        for (int velI = 0; velI < velDim; ++velI) {
                            in >> velocity[velI];
                        }
                    }
                }
                in.close();

                // check element inversion:
                IPC::Mesh<DIM> mesh(V0, F, SF, Eigen::MatrixXi(0, 2), V,
                    std::vector<int>(), std::vector<int>(), std::vector<int>(), std::vector<int>(),
                    std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>>(),
                    std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>>(),
                    std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>>(),
                    std::vector<std::pair<Eigen::Vector3i, std::array<Eigen::Vector3d, 2>>>(),
                    std::vector<DirichletBC>(),
                    std::vector<NeumannBC>(),
                    std::vector<std::pair<int, std::string>>(),
                    1e5, 0.4, 1e3);
                int numOfInv = 0, numOfBInv = 0;
                for (int i = 0; i < mesh.F.rows(); ++i) {
                    if (!mesh.checkInversion(i, true)) {
                        ++numOfInv;
                        if (mesh.isBoundaryVert(mesh.F(i, 0)) || mesh.isBoundaryVert(mesh.F(i, 1)) || mesh.isBoundaryVert(mesh.F(i, 2)) || mesh.isBoundaryVert(mesh.F(i, 3))) {
                            ++numOfBInv;
                        }
                    }
                }
                std::cout << numOfInv << " inversions out of " << mesh.F.rows() << " tetrahedra" << std::endl;
                std::cout << "where " << numOfBInv << " of the elements are on boundary" << std::endl;

                mesh.saveAsMesh("test.msh", false, SF);

                FILE* out = fopen("test.obj", "w");
                for (int i = 0; i < V.rows(); ++i) {
                    fprintf(out, "v %le %le %le\n", V(i, 0), V(i, 1), V(i, 2));
                }
                fclose(out);

                break;
            }

            case 20: { // check friction stencil
                std::cout << "check friction stencil" << std::endl;

                double eps = 1.0e-6;
                Eigen::RowVector3d v0PE(0.5, 1.0, -eps);
                Eigen::RowVector3d v0PT(0.5, 1.0, eps);
                Eigen::RowVector3d v0PT2(0.5, 1.0, 0.25);
                Eigen::RowVector3d v1(0.0, 0.0, 0.0);
                Eigen::RowVector3d v2(0.0, 0.0, 1.0);
                Eigen::RowVector3d v3(1.0, 0.0, 0.0);

                Eigen::Matrix<double, 12, 1> displacement_PT;
                displacement_PT.setZero();
                displacement_PT[2] = -1.0;
                displacement_PT[5] = 1.0;
                displacement_PT[8] = 1.0;
                displacement_PT[11] = 1.0;
                Eigen::Matrix<double, 9, 1> displacement_PE;
                displacement_PE.setZero();
                displacement_PE.head(6) = displacement_PT.head(6);
                displacement_PE.tail(3) = displacement_PT.tail(3);

                // PT1
                Eigen::Vector2d weight_PT;
                computeClosestPoint_PT(v0PT, v1, v2, v3, weight_PT);
                Eigen::Matrix<double, 3, 2> basis_PT;
                computeTangentBasis_PT(v0PT, v1, v2, v3, basis_PT);
                Eigen::Matrix<double, 12, 12> TTT_PT;
                computeTTT_PT(basis_PT, weight_PT[0], weight_PT[1], TTT_PT);
                Eigen::Matrix<double, 12, 1> grad_PT;
                g_PT(v0PT, v1, v2, v3, grad_PT);

                std::cout << "grad_PT1 = " << grad_PT << std::endl;
                std::cout << "basis_PT1 = " << basis_PT << std::endl;
                std::cout << "weight_PT1 = " << weight_PT << std::endl;
                std::cout << "fric_PT1 = " << grad_PT.norm() * TTT_PT * displacement_PT << std::endl;

                // PT2
                computeClosestPoint_PT(v0PT2, v1, v2, v3, weight_PT);
                computeTangentBasis_PT(v0PT2, v1, v2, v3, basis_PT);
                computeTTT_PT(basis_PT, weight_PT[0], weight_PT[1], TTT_PT);
                g_PT(v0PT2, v1, v2, v3, grad_PT);

                std::cout << "grad_PT2 = " << grad_PT << std::endl;
                std::cout << "basis_PT2 = " << basis_PT << std::endl;
                std::cout << "weight_PT2 = " << weight_PT << std::endl;
                std::cout << "fric_PT2 = " << grad_PT.norm() * TTT_PT * displacement_PT << std::endl;

                // PE
                double weight_PE;
                computeClosestPoint_PE(v0PE, v1, v3, weight_PE);
                Eigen::Matrix<double, 3, 2> basis_PE;
                computeTangentBasis_PE(v0PE, v1, v3, basis_PE);
                Eigen::Matrix<double, 9, 9> TTT_PE;
                computeTTT_PE(basis_PE, weight_PE, TTT_PE);
                Eigen::Matrix<double, 9, 1> grad_PE;
                g_PE(v0PE, v1, v3, grad_PE);

                std::cout << "grad_PE = " << grad_PE << std::endl;
                std::cout << "basis_PE = " << basis_PE << std::endl;
                std::cout << "weight_PE = " << weight_PE << std::endl;
                std::cout << "fric_PE = " << grad_PE.norm() * TTT_PE * displacement_PE << std::endl;

                break;
            }

            case 21: { // check computing EE cross energy
                derivTest_EECross();
                derivTest_e(1.0e+1);
                derivTest_e(1.0);
                break;
            }

#ifdef USE_AMGCL
            case 22: { // test AMGCL
                std::vector<int> ptr, col;
                std::vector<double> val, rhs;
                int n = poisson(1000, ptr, col, val, rhs);

                typedef amgcl::backend::builtin<double> Backend;

                // Use AMG as preconditioner:
                typedef amgcl::make_solver<
                    // Use AMG as preconditioner:
                    amgcl::amg<
                        Backend,
                        amgcl::coarsening::smoothed_aggregation,
                        amgcl::relaxation::spai0>,
                    // And BiCGStab as iterative solver:
                    amgcl::solver::bicgstab<Backend>>
                    Solver;

                Solver solve(std::tie(n, ptr, col, val));
                std::cout << solve.precond() << std::endl;

                std::vector<double> x(n, 0.0);
                int iters;
                double error;
                std::tie(iters, error) = solve(rhs, x);
                std::cout << "Iterations: " << iters << std::endl
                          << "Error:      " << error << std::endl;

                break;
            }
#endif

            case 23: { // test linear solver
                LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver
                    = LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>::create(IPC_DEFAULT_LINSYSSOLVER);

                std::vector<std::set<int>> vNeighbor(10);
                std::set<int> fixedVert;

                linSysSolver->set_pattern(vNeighbor, fixedVert);

                linSysSolver->setZero();
                for (int i = 0; i < vNeighbor.size() * DIM; ++i) {
                    linSysSolver->addCoeff(i, i, 10);
                }

                linSysSolver->analyze_pattern();

                linSysSolver->factorize();

                Eigen::VectorXd rhs(vNeighbor.size() * DIM), x;
                rhs.setOnes();
                linSysSolver->solve(rhs, x);

                std::cout << x << std::endl;

                break;
            }

            case 24: { // nearly parallel EE test
                Eigen::RowVector3d v0, v1, v2, v3, p0, p1, p2, p3;
                v0 << 5.80045067825167737219e-01, 1.29804572829927900024e+00, 4.11720521740031375479e-01;
                v1 << 2.64648850028170068427e-01, 1.03339687827110893181e+00, 4.11720521740029710145e-01;
                v2 << 3.15396217796997779814e-01, 1.61344194609627678005e+00, 4.11720521740033651437e-01;
                v3 << 1.95835346779988131006e-16, 1.34879309606810671163e+00, 4.11720521740031986102e-01;
                p0 << -1.92857301711554124974e-14, 1.46066226496805817429e-15, 4.25214990292590645001e-04;
                p1 << -6.08688292835163873411e-16, 4.82431649207783937812e-14, -4.25214990338334164650e-04;
                p2 << 5.44092882834986712421e-15, -3.69845570938928918346e-14, -4.25214990330826389867e-04;
                p3 << 1.99012418663158535008e-14, -9.01951096008029847640e-15, 4.25214990314199173770e-04;

                FILE* out = fopen("output/diag24.txt", "w");
                assert(out);
                for (int i = -1000; i <= 1000; ++i) {
                    double stepSize = i * 1.0e-9;
                    double dist;
                    Eigen::RowVector3d va = v0 + stepSize * p0;
                    Eigen::RowVector3d vb = v1 + stepSize * p1;
                    Eigen::RowVector3d vc = v2 + stepSize * p2;
                    Eigen::RowVector3d vd = v3 + stepSize * p3;
                    computeEdgeEdgeD(va, vb, vc, vd, dist);
                    fprintf(out, "%.20le %.20le %.20le %.20le\n", stepSize, dist,
                        (vb - va).cross(vd - vc).squaredNorm(), (vb - va).squaredNorm() * (vd - vc).squaredNorm());
                }
                fclose(out);

                break;
            }

            case 25: { // test exact CCD
                // 0.01 step size has been verified to not have any intersection
                Eigen::RowVector3d v0, v1, v2, v3, p0, p1, p2, p3;
                v0 << 1.27272800072121783366e+00, 8.04974006455849711017e-01, 6.70184760591664985263e-01;
                v1 << 1.09091003529400110139e+00, 8.76000651912412720002e-01, 5.94896702132349797232e-01;
                v2 << 1.17436223258123839663e+00, 8.34969460592872803062e-01, 6.38397253693771982164e-01;
                v3 << 1.14174988443290525986e+00, 8.83018816218291591191e-01, 5.87459106340128545121e-01;
                p0 << -1.08873390871259774953e-10, -1.35602114355536183116e-03, 8.80732558556843311490e-04;
                p1 << -3.17599195294374979933e-09, -7.34018788780557924004e-04, 1.44070755109990665012e-03;
                p2 << 1.00750930498716656340e-03, -1.54735808182153216847e-03, 1.59672278950023175531e-03;
                p3 << 1.53219607460276839930e-03, -2.16017817697321734202e-03, 3.07026181598997471855e-03;

                double stepSize = 1;
                while (ccd::edgeEdgeCCD(
                    v0, v1, v2, v3,
                    v0 + stepSize * p0, v1 + stepSize * p1, v2 + stepSize * p2, v3 + stepSize * p3,
                    ccd::CCDMethod::ROOT_PARITY)) {
                    stepSize /= 2.0;
                }

                std::cout << "output largest feasible stepSize = " << stepSize << std::endl;
                std::cout << "%% matlab code to visualize EE CCD: " << std::endl;
                std::cout << "v = [" << std::endl;
                printf("%.20le %.20le %.20le\n", v0[0], v0[1], v0[2]);
                printf("%.20le %.20le %.20le\n", v1[0], v1[1], v1[2]);
                printf("%.20le %.20le %.20le\n", v2[0], v2[1], v2[2]);
                printf("%.20le %.20le %.20le\n", v3[0], v3[1], v3[2]);
                printf("%.20le %.20le %.20le\n", v0[0] + stepSize * p0[0], v0[1] + stepSize * p0[1], v0[2] + stepSize * p0[2]);
                printf("%.20le %.20le %.20le\n", v1[0] + stepSize * p1[0], v1[1] + stepSize * p1[1], v1[2] + stepSize * p1[2]);
                printf("%.20le %.20le %.20le\n", v2[0] + stepSize * p2[0], v2[1] + stepSize * p2[1], v2[2] + stepSize * p2[2]);
                printf("%.20le %.20le %.20le\n", v3[0] + stepSize * p3[0], v3[1] + stepSize * p3[1], v3[2] + stepSize * p3[2]);
                std::cout << "];" << std::endl;
                std::cout << "hold on" << std::endl;
                std::cout << "view(3)" << std::endl;
                std::cout << "trisurf([1 2 5; 2 5 6], v(:,1), v(:,2), v(:,3), [0 0]);" << std::endl;
                std::cout << "trisurf([3 4 7; 4 7 8], v(:,1), v(:,2), v(:,3), [1 1]);" << std::endl;
                std::cout << "scatter3(v(1,1), v(1,2), v(1,3));" << std::endl;
                std::cout << "scatter3(v(2,1), v(2,2), v(2,3));" << std::endl;
                std::cout << "scatter3(v(3,1), v(3,2), v(3,3));" << std::endl;
                std::cout << "scatter3(v(4,1), v(4,2), v(4,3));" << std::endl;
                std::cout << "hold off" << std::endl;
                break;
            }

            case 26: { // ./src/Projects/Diagnostic/diagnostic 0 26
                double dHat = 0.8;

                Eigen::RowVector3d v0, v1, v2, v3;
                v1 << 0.0, 0.0, 1.0;
                v2 << -1.0, 0.0, -1.0;
                v3 << 1.0, 0.0, -1.0;
                v0 << 0.0, 0.5, 0.0;

                double d;
                std::cout << "distance function value and derivatives:" << std::endl;
                switch (dType_PT(v0, v1, v2, v3)) {
                case 0: {
                    d_PP(v0, v1, d);
                    Eigen::Matrix<double, 6, 1> g;
                    Eigen::Matrix<double, 6, 6> h;
                    g_PP(v0, v1, g);
                    H_PP(h);
                    std::cout << d << std::endl;
                    std::cout << g << std::endl
                              << h << std::endl;
                    break;
                }

                case 1: {
                    d_PP(v0, v2, d);
                    Eigen::Matrix<double, 6, 1> g;
                    Eigen::Matrix<double, 6, 6> h;
                    g_PP(v0, v2, g);
                    H_PP(h);
                    std::cout << d << std::endl;
                    std::cout << g << std::endl
                              << h << std::endl;
                    break;
                }

                case 2: {
                    d_PP(v0, v3, d);
                    Eigen::Matrix<double, 6, 1> g;
                    Eigen::Matrix<double, 6, 6> h;
                    g_PP(v0, v3, g);
                    H_PP(h);
                    std::cout << d << std::endl;
                    std::cout << g << std::endl
                              << h << std::endl;
                    break;
                }

                case 3: {
                    d_PE(v0, v1, v2, d);
                    Eigen::Matrix<double, 9, 1> g;
                    Eigen::Matrix<double, 9, 9> h;
                    g_PE(v0, v1, v2, g);
                    H_PE(v0, v1, v2, h);
                    std::cout << d << std::endl;
                    std::cout << g << std::endl
                              << h << std::endl;
                    break;
                }

                case 4: {
                    d_PE(v0, v2, v3, d);
                    Eigen::Matrix<double, 9, 1> g;
                    Eigen::Matrix<double, 9, 9> h;
                    g_PE(v0, v2, v3, g);
                    H_PE(v0, v2, v3, h);
                    std::cout << d << std::endl;
                    std::cout << g << std::endl
                              << h << std::endl;
                    break;
                }

                case 5: {
                    d_PE(v0, v3, v1, d);
                    Eigen::Matrix<double, 9, 1> g;
                    Eigen::Matrix<double, 9, 9> h;
                    g_PE(v0, v3, v1, g);
                    H_PE(v0, v3, v1, h);
                    std::cout << d << std::endl;
                    std::cout << g << std::endl
                              << h << std::endl;
                    break;
                }

                case 6: {
                    d_PT(v0, v1, v2, v3, d);
                    Eigen::Matrix<double, 12, 1> g;
                    Eigen::Matrix<double, 12, 12> h;
                    g_PT(v0, v1, v2, v3, g);
                    H_PT(v0, v1, v2, v3, h);
                    std::cout << d << std::endl;
                    std::cout << g << std::endl
                              << h << std::endl;
                    break;
                }

                default:
                    d = -1.0;
                    break;
                }

                double b, g_b, h_b;
                compute_b(d, dHat * dHat, b);
                compute_g_b(d, dHat * dHat, g_b);
                compute_H_b(d, dHat * dHat, h_b);
                std::cout << "barrier function value and derivatives:" << std::endl;
                std::cout << b << std::endl
                          << g_b << std::endl
                          << h_b << std::endl;

                break;
            }

            case 101: { // ./src/Projects/Diagnostic/diagnostic 0 101 output/resultFolder 100
                if (argc < 5) {
                    std::cout << "need result folder path and number of time steps" << std::endl;
                    exit(0);
                }

                std::ofstream statsFile(std::string(argv[3]) + "/perTSStats.txt");
                if (!statsFile.is_open()) {
                    std::cout << "folder " << argv[3] << "does not exist" << std::endl;
                }

                int numOfTS = std::stoi(argv[4]);
                for (int tsI = 0; tsI < numOfTS; ++tsI) {
                    std::ifstream infoFileI(std::string(argv[3]) + "/info" + std::to_string(tsI) + ".txt");
                    if (infoFileI.is_open()) {
                        std::string line;
                        std::stringstream ss;
                        std::string bypass;

                        std::getline(infoFileI, line);
                        std::getline(infoFileI, line);
                        ss = std::stringstream(line);
                        int tsNum;
                        double iterAmt;
                        ss >> tsNum >> iterAmt;
                        assert(tsNum == tsI);
                        statsFile << tsNum << " " << iterAmt << " ";

                        std::getline(infoFileI, line);
                        ss = std::stringstream(line);
                        int curMem, peakMem;
                        ss >> bypass >> curMem >> bypass >> bypass >> peakMem;
                        std::getline(infoFileI, line);
                        std::getline(infoFileI, line);
                        std::getline(infoFileI, line);
                        ss = std::stringstream(line);
                        double curRunTime;
                        ss >> curRunTime;
                        statsFile << curRunTime << " ";
                        statsFile << curMem << " " << peakMem << " ";

                        statsFile << std::endl;
                        infoFileI.close();
                        std::cout << "finished reading info" << tsI << ".txt" << std::endl;
                    }
                    else {
                        std::cout << "info" << tsI << ".txt does not exist" << std::endl;
                    }
                }

                statsFile.close();
                break;
            }

            case 102: { // ./src/Projects/Diagnostic/diagnostic 0 102 input/ 0 200
                if (argc < 6) {
                    std::cout << "need input seg file folder path, start frame, and end frame" << std::endl;
                    exit(0);
                }

                int startFrame = std::stoi(argv[4]), endFrame = std::stoi(argv[5]);
                for (int i = startFrame; i <= endFrame; ++i) {
                    std::ifstream segFile((std::string(argv[3]) + "/" + std::to_string(i) + ".seg").c_str());
                    if (!segFile.is_open()) {
                        std::cout << argv[3] << " does not exist" << std::endl;
                        continue;
                    }

                    std::ofstream objFile((std::string(argv[3]) + "/seg" + std::to_string(i) + ".obj").c_str());
                    assert(objFile.is_open());

                    std::string line;
                    while (std::getline(segFile, line)) {
                        std::stringstream ss(line);
                        std::string key;
                        ss >> key;
                        if (key == "v") {
                            objFile << line << std::endl;
                        }
                        else if (key == "s") {
                            int v0I, v1I;
                            ss >> v0I >> v1I;
                            objFile << "f " << v0I << " " << v1I << " " << v1I << std::endl;
                        }
                    }

                    segFile.close();
                    objFile.close();
                }
                break;
            }

            // msh to vtk
            case 104: { // ./src/Projects/Diagnostic/diagnostic 0 104 input/tetMeshes/mat20x20.msh
                if (argc < 4) {
                    std::cout << "need input tetmesh file path" << std::endl;
                    exit(0);
                }

                Eigen::MatrixXd TV;
                Eigen::MatrixXi TT, F;
                IglUtils::readTetMesh(argv[3], TV, TT, F);
                std::string saveFilePath = std::string(argv[3]).substr(0, std::string(argv[3]).find_last_of('.')) + ".vtk";
                IglUtils::saveTetMesh_vtk(saveFilePath, TV, TT);

                break;
            }

            case 105: { // test linear solver with loaded matrix
                // ./src/Projects/Diagnostic/diagnostic 0 105 Mnan
                if (argc < 4) {
                    std::cout << "need linear system filepath" << std::endl;
                    exit(-1);
                }

                LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver
                    = LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>::create(IPC_DEFAULT_LINSYSSOLVER);
                Eigen::VectorXd rhs, x;

                // linSysSolver->load(argv[3], rhs);
                // IglUtils::writeSparseMatrixToFile("M_MATLAB", linSysSolver, true);

#ifdef IPC_WITH_AMGCL
                if (IPC_DEFAULT_LINSYSSOLVER == LinSysSolverType::AMGCL) {
                    dynamic_cast<AMGCLSolver<Eigen::VectorXi, Eigen::VectorXd>*>(linSysSolver)->load_AMGCL(argv[3], rhs);
                    dynamic_cast<AMGCLSolver<Eigen::VectorXi, Eigen::VectorXd>*>(linSysSolver)->write_AMGCL((std::string(argv[3]) + "check").c_str(), rhs);
                }
                else {
#endif
                    linSysSolver->load(argv[3], rhs);
                    IglUtils::writeSparseMatrixToFile(std::string(argv[3]) + "check", linSysSolver, true);
#ifdef IPC_WITH_AMGCL
                }
#endif

                linSysSolver->analyze_pattern();

                linSysSolver->factorize();

                linSysSolver->solve(rhs, x);

                std::cout << x.squaredNorm() << std::endl;

                break;
            }

            default:
                std::cout << "No diagMode " << diagMode << std::endl;
                break;
            }
        }
        else {
            std::cout << "Please enter diagMode!" << std::endl;
        }
    }

    static int poisson(
        int n,
        std::vector<int>& ptr,
        std::vector<int>& col,
        std::vector<double>& val,
        std::vector<double>& rhs)
    {
        int n2 = n * n; // Number of points in the grid.
        double h = 1.0 / (n - 1); // Grid spacing.

        ptr.clear();
        ptr.reserve(n2 + 1);
        ptr.push_back(0);
        col.clear();
        col.reserve(n2 * 5); // We use 5-point stencil, so the matrix
        val.clear();
        val.reserve(n2 * 5); // will have at most n2 * 5 nonzero elements.

        rhs.resize(n2);

        for (int j = 0, k = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i, ++k) {
                if (i == 0 || i == n - 1 || j == 0 || j == n - 1) {
                    // Boundary point. Use Dirichlet condition.
                    col.push_back(k);
                    val.push_back(1.0);

                    rhs[k] = 0.0;
                }
                else {
                    // Interior point. Use 5-point finite difference stencil.
                    col.push_back(k - n);
                    val.push_back(-1.0 / (h * h));

                    col.push_back(k - 1);
                    val.push_back(-1.0 / (h * h));

                    col.push_back(k);
                    val.push_back(4.0 / (h * h));

                    col.push_back(k + 1);
                    val.push_back(-1.0 / (h * h));

                    col.push_back(k + n);
                    val.push_back(-1.0 / (h * h));

                    rhs[k] = 1.0;
                }

                ptr.push_back(col.size());
            }
        }

        return n2;
    }
};
} // namespace IPC

int main(int argc, char* argv[])
{
    IPC::Diagnostic::run(argc, argv);
    return 0;
}

#endif /* Diagnostic_hpp */
