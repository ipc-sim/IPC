//
//  MeshProcessing.hpp
//  IPC
//
//  Created by Minchen Li on 1/31/18.
//

#include "IglUtils.hpp"

#include "Timer.hpp"

#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/euler_characteristic.h>
#include <igl/per_vertex_normals.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

#ifdef USE_TBB
#include <tbb/parallel_for.h>
#endif

#include <cstdio>

std::string outputFolderPath = "output/";

namespace IPC {
class MeshProcessing {
public:
    static void run(int argc, char* argv[])
    {
        if (argc > 2) {
            Eigen::MatrixXd V, UV, N, TV;
            Eigen::MatrixXi F, FUV, FN, TT;

            const std::string meshPath = std::string(argv[2]);
            const std::string meshFolderPath = meshPath.substr(0, meshPath.find_last_of('/') + 1);
            const std::string meshFileName = meshPath.substr(meshPath.find_last_of('/') + 1);
            const std::string meshName = meshFileName.substr(0, meshFileName.find_last_of('.'));
            const std::string suffix = meshFileName.substr(meshFileName.find_last_of('.'));

            if (suffix == ".off") {
                igl::readOFF(meshPath, V, F);
            }
            else if (suffix == ".obj") {
                igl::readOBJ(meshPath, V, UV, N, F, FUV, FN);
            }
            else if (suffix == ".msh") {
                int procMode = std::stoi(argv[3]);
                bool findSurface = true;
                if ((procMode == 3) || (procMode == 4)) {
                    findSurface = false;
                }
                IglUtils::readTetMesh(meshPath, TV, TT, F, findSurface);
            }
            else if (suffix == ".tet") {
                std::cout << "verschoor's scene format" << std::endl;
            }
            else {
                std::cout << "unkown mesh file format!" << std::endl;
                return;
            }

            if (argc > 3) {
                int procMode = 0;
                procMode = std::stoi(argv[3]);
                switch (procMode) {

                // ./src/Projects/MeshProcessing/meshprocessing 0 input/triMeshes/mat100x100.obj 3 1e10 0
                case 3: { // tetrahedralize the input surface mesh
                    // Tetrahedralized interior
                    Eigen::MatrixXd TV;
                    Eigen::MatrixXi TT;
                    Eigen::MatrixXi TF;

                    using namespace Eigen;
                    using namespace std;

                    // Tetrahedralize the interior
                    std::string flag("pq1.414");
                    if (argc > 4) {
                        double maxElemVol = std::stod(argv[4]);
                        if (maxElemVol > 0.0) {
                            flag += std::string("a") + std::string(argv[4]);
                        }
                    }
                    if (argc > 5) {
                        int addSteinerPoints = stod(argv[5]);
                        if (!addSteinerPoints) {
                            flag += 'Y';
                        }
                    }
                    else {
                        flag += 'Y';
                    }
                    igl::copyleft::tetgen::tetrahedralize(V, F, flag.c_str(), TV, TT, TF);

                    IglUtils::saveTetMesh((meshFolderPath + meshName + ".msh").c_str(),
                        TV, TT);

                    break;
                }

                case 19: { // ./src/Projects/MeshProcessing/meshprocessing 0 armadillo.tet 19
                    std::ifstream file(argv[2]);
                    if (!file.is_open()) {
                        std::cout << "file does not exist" << std::endl;
                        exit(-1);
                    }

                    int n_v, n_T, bypass;
                    file >> n_v >> n_T >> bypass >> bypass >> bypass;

                    Eigen::MatrixXd V(n_v, 3);
                    Eigen::MatrixXi F(n_T, 4);
                    for (int i = 0; i < n_v; ++i) {
                        file >> V(i, 0) >> V(i, 1) >> V(i, 2) >> bypass;
                    }
                    for (int i = 0; i < n_T; ++i) {
                        file >> F(i, 0) >> F(i, 1) >> F(i, 2) >> F(i, 3);
                    }

                    IglUtils::saveTetMesh("tetMesh.msh", V, F);

                    file.close();
                    break;
                }

                case 20: { // ./src/Projects/MeshProcessing/meshprocessing 0 armadillo.obj 20 armadillo.msh
                    Eigen::MatrixXd TV;
                    Eigen::MatrixXi TT, TF;
                    IglUtils::readTetMesh(argv[4], TV, TT, TF, true);

                    Eigen::Matrix3d rotMtr = (Eigen::AngleAxisd(0.0, Eigen::Vector3d::UnitX())
                        * Eigen::AngleAxisd(M_PI / 2.0, Eigen::Vector3d::UnitY())
                        * Eigen::AngleAxisd(0.0, Eigen::Vector3d::UnitZ()))
                                                 .toRotationMatrix();
                    for (int i = 0; i < TV.rows(); ++i) {
                        TV.row(i) = TV.row(i) * rotMtr.transpose();
                    }

                    Eigen::RowVector3d max_t = TV.colwise().maxCoeff();
                    Eigen::RowVector3d min_t = TV.colwise().minCoeff();
                    Eigen::RowVector3d center_t = (max_t + min_t) / 2.0;
                    Eigen::RowVector3d max_s = V.colwise().maxCoeff();
                    Eigen::RowVector3d min_s = V.colwise().minCoeff();
                    Eigen::RowVector3d center_s = (max_s + min_s) / 2.0;
                    double scale = (max_s[0] - min_s[0]) / (max_t[0] - min_t[0]);
                    for (int i = 0; i < TV.rows(); ++i) {
                        TV.row(i) = scale * (TV.row(i) - center_t) + center_s;
                    }

                    IglUtils::saveTetMesh("tetMesh.msh", TV, TT, TF);

                    break;
                }

                default:
                    std::cout << "No procMode " << procMode << std::endl;
                    break;
                }
            }
            else {
                std::cout << "Please enter procMode!" << std::endl;
            }
        }
        else {
            std::cout << "Please enter mesh file path!" << std::endl;
        }
    }
};

} // namespace IPC

int main(int argc, char* argv[])
{
    // mesh processing mode
    IPC::MeshProcessing::run(argc, argv);
    return 0;
}
