//
//  Mesh.cpp
//  IPC
//
//  Created by Minchen Li on 8/30/17.
//

#include "Mesh.hpp"
#include "IglUtils.hpp"
#include "Optimizer.hpp"
#include "Timer.hpp"
#include "CCDUtils.hpp"

#include <igl/triangle/triangulate.h>
#include <igl/cotmatrix.h>
#include <igl/avg_edge_length.h>
#include <igl/writeOBJ.h>
#include <igl/list_to_matrix.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>

#ifdef USE_TBB
#include <tbb/tbb.h>
#endif

#include <fstream>
#include <deque>
#include <unordered_map>

#include <spdlog/spdlog.h>

extern std::ofstream logFile;
extern Timer timer_step;

namespace IPC {

template <int dim>
Mesh<dim>::Mesh(void)
{
}

template <int dim>
Mesh<dim>::Mesh(const Eigen::MatrixXd& V_mesh,
    const Eigen::MatrixXi& F_mesh, const Eigen::MatrixXi& SF_mesh,
    const Eigen::MatrixXi& CE_mesh, const Eigen::MatrixXd& Vt_mesh,
    const std::vector<int>& p_componentNodeRange,
    const std::vector<int>& p_componentSFRange,
    const std::vector<int>& p_componentCERange,
    const std::vector<int>& p_componentCoDim,
    const std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>>& p_componentMaterial,
    const std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>>& p_componentLVels,
    const std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>>& p_componentAVels,
    const std::vector<std::pair<Eigen::Vector3i, std::array<Eigen::Vector3d, 2>>>& p_componentInitVels,
    const std::vector<std::pair<std::vector<int>, std::array<Eigen::Vector3d, 2>>>& p_DBCInfo,
    const std::map<int, Eigen::Matrix<double, 1, dim>>& p_NeumannBC,
    const std::vector<std::pair<int, std::string>>& p_meshSeqFolderPath,
    double YM, double PR, double rho)
{
    assert(V_mesh.rows() > 0);
    assert(F_mesh.rows() > 0);

    V_rest = V_mesh;
    F = F_mesh;
    SF = SF_mesh;
    CE = CE_mesh;
    componentNodeRange = p_componentNodeRange;
    componentSFRange = p_componentSFRange;
    componentCERange = p_componentCERange;
    componentCoDim = p_componentCoDim;
    componentMaterial = p_componentMaterial;
    componentLVels = p_componentLVels;
    componentAVels = p_componentAVels;
    componentInitVels = p_componentInitVels;
    DBCInfo = p_DBCInfo;
    NeumannBC = p_NeumannBC;
    meshSeqFolderPath = p_meshSeqFolderPath;
    if (Vt_mesh.rows() == V_mesh.rows()) {
        V = Vt_mesh;
    }
    else {
        assert(Vt_mesh.rows() == 0);
        V = Eigen::MatrixXd::Zero(V_rest.rows(), dim);
        spdlog::warn("No Vt provided, initialized to all 0");
    }

    density = rho;
    m_YM = YM;
    m_PR = PR;
    computeFeatures(false, true);

    // initialize default lame parameters per element
    setLameParam(YM, PR);
}

void initCylinder(double r1_x, double r1_y, double r2_x, double r2_y, double height, int circle_res, int height_resolution,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd* uv_coords_per_face = NULL,
    Eigen::MatrixXi* uv_coords_face_ids = NULL)
{
    assert(DIM == 2);
#if (DIM == 2)
    int nvertices = circle_res * (height_resolution + 1);
    int nfaces = 2 * circle_res * height_resolution;

    V.resize(nvertices, 3);
    if (uv_coords_per_face) {
        uv_coords_per_face->resize(nvertices, 2);
    }
    F.resize(nfaces, 3);
    for (int j = 0; j < height_resolution + 1; j++) {
        for (int i = 0; i < circle_res; i++) {
            double t = (double)j / (double)height_resolution;
            double h = height * t;
            double theta = i * 2 * M_PI / circle_res;
            double r_x = r1_x * t + r2_x * (1 - t);
            double r_y = r1_y * t + r2_y * (1 - t);
            V.row(j * circle_res + i) = Eigen::Vector3d(r_x * cos(theta), height - h, r_y * sin(theta));
            if (uv_coords_per_face) {
                uv_coords_per_face->row(j * circle_res + i) = Eigen::Vector2d(r_x * cos(theta), r_y * sin(theta));
            }

            if (j < height_resolution) {
                int vl0 = j * circle_res + i;
                int vl1 = j * circle_res + (i + 1) % circle_res;
                int vu0 = (j + 1) * circle_res + i;
                int vu1 = (j + 1) * circle_res + (i + 1) % circle_res;
                F.row(2 * (j * circle_res + i) + 0) = Eigen::Vector3i(vl0, vl1, vu1);
                F.row(2 * (j * circle_res + i) + 1) = Eigen::Vector3i(vu0, vl0, vu1);
            }
        }
    }
#endif
}

template <int dim>
void Mesh<dim>::computeMassMatrix(const igl::MassMatrixType type)
{
    const Eigen::MatrixXd& V = V_rest;
    Eigen::SparseMatrix<double>& M = massMatrix;

    typedef double Scalar;

    using namespace Eigen;
    using namespace std;
    using namespace igl;

    const int n = V.rows();
    const int m = F.rows();
    const int simplex_size = F.cols();

    MassMatrixType eff_type = type;
    // Use voronoi of for triangles by default, otherwise barycentric
    if (type == MASSMATRIX_TYPE_DEFAULT) {
        eff_type = (simplex_size == 3 ? MASSMATRIX_TYPE_VORONOI : MASSMATRIX_TYPE_BARYCENTRIC);
    }

    // Not yet supported
    assert(type != MASSMATRIX_TYPE_FULL);

    Matrix<int, Dynamic, 1> MI;
    Matrix<int, Dynamic, 1> MJ;
    Matrix<Scalar, Dynamic, 1> MV;
    if (simplex_size == 3) {
        // Triangles
        // edge lengths numbered same as opposite vertices
        Matrix<Scalar, Dynamic, 3> l(m, 3);
        // loop over faces
        for (int i = 0; i < m; i++) {
            l(i, 0) = (V.row(F(i, 1)) - V.row(F(i, 2))).norm();
            l(i, 1) = (V.row(F(i, 2)) - V.row(F(i, 0))).norm();
            l(i, 2) = (V.row(F(i, 0)) - V.row(F(i, 1))).norm();
        }
        Matrix<Scalar, Dynamic, 1> dblA;
        doublearea(l, 0., dblA);

        switch (eff_type) {
        case MASSMATRIX_TYPE_BARYCENTRIC:
            // diagonal entries for each face corner
            MI.resize(m * 3, 1);
            MJ.resize(m * 3, 1);
            MV.resize(m * 3, 1);
            MI.block(0 * m, 0, m, 1) = F.col(0);
            MI.block(1 * m, 0, m, 1) = F.col(1);
            MI.block(2 * m, 0, m, 1) = F.col(2);
            MJ = MI;
            repmat(dblA, 3, 1, MV);
            MV.array() /= 6.0;
            break;
        case MASSMATRIX_TYPE_VORONOI: {
            // diagonal entries for each face corner
            // http://www.alecjacobson.com/weblog/?p=874
            MI.resize(m * 3, 1);
            MJ.resize(m * 3, 1);
            MV.resize(m * 3, 1);
            MI.block(0 * m, 0, m, 1) = F.col(0);
            MI.block(1 * m, 0, m, 1) = F.col(1);
            MI.block(2 * m, 0, m, 1) = F.col(2);
            MJ = MI;

            // Holy shit this needs to be cleaned up and optimized
            Matrix<Scalar, Dynamic, 3> cosines(m, 3);
            cosines.col(0) = (l.col(2).array().pow(2) + l.col(1).array().pow(2) - l.col(0).array().pow(2)) / (l.col(1).array() * l.col(2).array() * 2.0);
            cosines.col(1) = (l.col(0).array().pow(2) + l.col(2).array().pow(2) - l.col(1).array().pow(2)) / (l.col(2).array() * l.col(0).array() * 2.0);
            cosines.col(2) = (l.col(1).array().pow(2) + l.col(0).array().pow(2) - l.col(2).array().pow(2)) / (l.col(0).array() * l.col(1).array() * 2.0);
            Matrix<Scalar, Dynamic, 3> barycentric = cosines.array() * l.array();
            normalize_row_sums(barycentric, barycentric);
            Matrix<Scalar, Dynamic, 3> partial = barycentric;
            partial.col(0).array() *= dblA.array() * 0.5;
            partial.col(1).array() *= dblA.array() * 0.5;
            partial.col(2).array() *= dblA.array() * 0.5;
            Matrix<Scalar, Dynamic, 3> quads(partial.rows(), partial.cols());
            quads.col(0) = (partial.col(1) + partial.col(2)) * 0.5;
            quads.col(1) = (partial.col(2) + partial.col(0)) * 0.5;
            quads.col(2) = (partial.col(0) + partial.col(1)) * 0.5;

            quads.col(0) = (cosines.col(0).array() < 0).select(0.25 * dblA, quads.col(0));
            quads.col(1) = (cosines.col(0).array() < 0).select(0.125 * dblA, quads.col(1));
            quads.col(2) = (cosines.col(0).array() < 0).select(0.125 * dblA, quads.col(2));

            quads.col(0) = (cosines.col(1).array() < 0).select(0.125 * dblA, quads.col(0));
            quads.col(1) = (cosines.col(1).array() < 0).select(0.25 * dblA, quads.col(1));
            quads.col(2) = (cosines.col(1).array() < 0).select(0.125 * dblA, quads.col(2));

            quads.col(0) = (cosines.col(2).array() < 0).select(0.125 * dblA, quads.col(0));
            quads.col(1) = (cosines.col(2).array() < 0).select(0.125 * dblA, quads.col(1));
            quads.col(2) = (cosines.col(2).array() < 0).select(0.25 * dblA, quads.col(2));

            MV.block(0 * m, 0, m, 1) = quads.col(0);
            MV.block(1 * m, 0, m, 1) = quads.col(1);
            MV.block(2 * m, 0, m, 1) = quads.col(2);

            break;
        }
        case MASSMATRIX_TYPE_FULL:
            assert(false && "Implementation incomplete");
            break;
        default:
            assert(false && "Unknown Mass matrix eff_type");
        }
    }
    else if (simplex_size == 4) {
        assert(V.cols() == 3);
        if (eff_type != MASSMATRIX_TYPE_BARYCENTRIC) {
            spdlog::info("switched to barycentric mass matrix for tetrahedra");
        }
        MI.resize(m * 4, 1);
        MJ.resize(m * 4, 1);
        MV.resize(m * 4, 1);
        MI.block(0 * m, 0, m, 1) = F.col(0);
        MI.block(1 * m, 0, m, 1) = F.col(1);
        MI.block(2 * m, 0, m, 1) = F.col(2);
        MI.block(3 * m, 0, m, 1) = F.col(3);
        MJ = MI;
        // loop over tets
        for (int i = 0; i < m; i++) {
            // http://en.wikipedia.org/wiki/Tetrahedron#Volume
            Matrix<Scalar, 3, 1> v0m3, v1m3, v2m3;
            v0m3.head(V.cols()) = V.row(F(i, 0)) - V.row(F(i, 3));
            v1m3.head(V.cols()) = V.row(F(i, 1)) - V.row(F(i, 3));
            v2m3.head(V.cols()) = V.row(F(i, 2)) - V.row(F(i, 3));
            Scalar v = fabs(v0m3.dot(v1m3.cross(v2m3))) / 6.0;
            MV(i + 0 * m) = v / 4.0;
            MV(i + 1 * m) = v / 4.0;
            MV(i + 2 * m) = v / 4.0;
            MV(i + 3 * m) = v / 4.0;
        }
    }
    else {
        // Unsupported simplex size
        assert(false && "Unsupported simplex size");
    }

    if constexpr (dim == 3) {
        // for surface mesh as kinematic object
        for (int compI = 0; compI < componentCoDim.size(); ++compI) {
            if (componentCoDim[compI] == 3) {
                continue;
            }
            else if (componentCoDim[compI] == 1) {
                int m = componentCERange[compI + 1] - componentCERange[compI];
                int oldRows = MI.rows();

                MI.conservativeResize(oldRows + m * 2, 1);
                MI.block(oldRows, 0, m, 1) = CE.col(0).segment(componentCERange[compI], m);
                MI.block(oldRows + m, 0, m, 1) = CE.col(1).segment(componentCERange[compI], m);

                MJ.conservativeResize(oldRows + m * 2, 1);
                MJ.bottomRows(m * 2) = MI.bottomRows(m * 2);

                MV.conservativeResize(oldRows + m * 2, 1);
                for (int i = 0; i < m; ++i) {
                    double l = (V.row(CE(componentCERange[compI] + i, 0)) - V.row(CE(componentCERange[compI] + i, 1))).norm();
                    MV[oldRows + i] = MV[oldRows + m + i] = l * l * l * M_PI / 12.0;
                }
            }
            else if (componentCoDim[compI] == 0) {
                int m = componentCERange[compI + 1] - componentCERange[compI];
                int oldRows = MI.rows();

                MI.conservativeResize(oldRows + m, 1);
                MI.bottomRows(m) = Eigen::VectorXi::LinSpaced(m, componentCERange[compI],
                    componentCERange[compI + 1] - 1);

                MJ.conservativeResize(oldRows + m, 1);
                MJ.bottomRows(m) = MI.bottomRows(m);

                MV.conservativeResize(oldRows + m, 1);
                MV.bottomRows(m).setZero(); // to set value later
            }
            else { // componentCoDim[compI] == 2

                int m = componentSFRange[compI + 1] - componentSFRange[compI];
                int oldRows = MI.rows();

                MI.conservativeResize(oldRows + m * 3, 1);
                MI.block(oldRows, 0, m, 1) = SF.col(0).segment(componentSFRange[compI], m);
                MI.block(oldRows + m, 0, m, 1) = SF.col(1).segment(componentSFRange[compI], m);
                MI.block(oldRows + m * 2, 0, m, 1) = SF.col(2).segment(componentSFRange[compI], m);

                MJ.conservativeResize(oldRows + m * 3, 1);
                MJ.bottomRows(m * 3) = MI.bottomRows(m * 3);

                MV.conservativeResize(oldRows + m * 3, 1);

                // edge lengths numbered same as opposite vertices
                Matrix<Scalar, Dynamic, 3> l(m, 3);
                // loop over faces
                for (int i = 0; i < m; i++) {
                    l(i, 0) = (V.row(SF(componentSFRange[compI] + i, 1)) - V.row(SF(componentSFRange[compI] + i, 2))).norm();
                    l(i, 1) = (V.row(SF(componentSFRange[compI] + i, 2)) - V.row(SF(componentSFRange[compI] + i, 0))).norm();
                    l(i, 2) = (V.row(SF(componentSFRange[compI] + i, 0)) - V.row(SF(componentSFRange[compI] + i, 1))).norm();
                }
                Matrix<Scalar, Dynamic, 1> dblA;
                doublearea(l, 0., dblA);

                switch (eff_type) {
                case MASSMATRIX_TYPE_BARYCENTRIC: {
                    // diagonal entries for each face corner
                    Eigen::MatrixXd curMV;
                    repmat(dblA, 3, 1, curMV);
                    curMV.array() /= 6.0;

                    MV.bottomRows(m * 3) = curMV;
                    break;
                }
                case MASSMATRIX_TYPE_VORONOI: {
                    // diagonal entries for each face corner
                    // http://www.alecjacobson.com/weblog/?p=874

                    // Holy shit this needs to be cleaned up and optimized
                    Matrix<Scalar, Dynamic, 3> cosines(m, 3);
                    cosines.col(0) = (l.col(2).array().pow(2) + l.col(1).array().pow(2) - l.col(0).array().pow(2)) / (l.col(1).array() * l.col(2).array() * 2.0);
                    cosines.col(1) = (l.col(0).array().pow(2) + l.col(2).array().pow(2) - l.col(1).array().pow(2)) / (l.col(2).array() * l.col(0).array() * 2.0);
                    cosines.col(2) = (l.col(1).array().pow(2) + l.col(0).array().pow(2) - l.col(2).array().pow(2)) / (l.col(0).array() * l.col(1).array() * 2.0);
                    Matrix<Scalar, Dynamic, 3> barycentric = cosines.array() * l.array();
                    normalize_row_sums(barycentric, barycentric);
                    Matrix<Scalar, Dynamic, 3> partial = barycentric;
                    partial.col(0).array() *= dblA.array() * 0.5;
                    partial.col(1).array() *= dblA.array() * 0.5;
                    partial.col(2).array() *= dblA.array() * 0.5;
                    Matrix<Scalar, Dynamic, 3> quads(partial.rows(), partial.cols());
                    quads.col(0) = (partial.col(1) + partial.col(2)) * 0.5;
                    quads.col(1) = (partial.col(2) + partial.col(0)) * 0.5;
                    quads.col(2) = (partial.col(0) + partial.col(1)) * 0.5;

                    quads.col(0) = (cosines.col(0).array() < 0).select(0.25 * dblA, quads.col(0));
                    quads.col(1) = (cosines.col(0).array() < 0).select(0.125 * dblA, quads.col(1));
                    quads.col(2) = (cosines.col(0).array() < 0).select(0.125 * dblA, quads.col(2));

                    quads.col(0) = (cosines.col(1).array() < 0).select(0.125 * dblA, quads.col(0));
                    quads.col(1) = (cosines.col(1).array() < 0).select(0.25 * dblA, quads.col(1));
                    quads.col(2) = (cosines.col(1).array() < 0).select(0.125 * dblA, quads.col(2));

                    quads.col(0) = (cosines.col(2).array() < 0).select(0.125 * dblA, quads.col(0));
                    quads.col(1) = (cosines.col(2).array() < 0).select(0.125 * dblA, quads.col(1));
                    quads.col(2) = (cosines.col(2).array() < 0).select(0.25 * dblA, quads.col(2));

                    MV.block(oldRows + 0 * m, 0, m, 1) = quads.col(0);
                    MV.block(oldRows + 1 * m, 0, m, 1) = quads.col(1);
                    MV.block(oldRows + 2 * m, 0, m, 1) = quads.col(2);

                    break;
                }
                case MASSMATRIX_TYPE_FULL:
                    assert(false && "Implementation incomplete");
                    break;
                default:
                    assert(false && "Unknown Mass matrix eff_type");
                }

                Eigen::VectorXd avgELen = l.rowwise().mean() / 3.0; // multiply it on pi*r^2 to approximate 2/3*pi*r^3
                MV.block(oldRows + 0 * m, 0, m, 1) = MV.block(oldRows + 0 * m, 0, m, 1).cwiseProduct(avgELen);
                MV.block(oldRows + 1 * m, 0, m, 1) = MV.block(oldRows + 1 * m, 0, m, 1).cwiseProduct(avgELen);
                MV.block(oldRows + 2 * m, 0, m, 1) = MV.block(oldRows + 2 * m, 0, m, 1).cwiseProduct(avgELen);
            }
        }
    }

    sparse(MI, MJ, MV, n, n, M);

    massMatrix *= density;
    // IglUtils::writeSparseMatrixToFile("mass", massMatrix, true);

    double avgMass = avgNodeMass(dim);
    for (int compI = 0; compI < componentCoDim.size(); ++compI) {
        if (componentCoDim[compI] == 0) {
            for (int vI = componentNodeRange[compI]; vI < componentNodeRange[compI + 1]; ++vI) {
                massMatrix.coeffRef(vI, vI) = avgMass;
            }
        }
    }
}

template <int dim>
void Mesh<dim>::computeFeatures(bool multiComp, bool resetFixedV)
{
    if (resetFixedV) {
        fixedVert.clear();
        fixedVert.insert(0);
        isFixedVert.resize(0);
        isFixedVert.resize(V.rows(), false);
        isFixedVert[0] = true;
    }

    restTriInv.resize(F.rows());
    triArea.resize(F.rows());
    vFLoc.resize(0);
    vFLoc.resize(V.rows());
    std::vector<Eigen::RowVector3d> vertNormals(V_rest.rows(), Eigen::RowVector3d::Zero());
    for (int triI = 0; triI < F.rows(); triI++) {
        const Eigen::Matrix<int, 1, dim + 1>& triVInd = F.row(triI);

        vFLoc[triVInd[0]].insert(std::pair<int, int>(triI, 0));
        vFLoc[triVInd[1]].insert(std::pair<int, int>(triI, 1));
        vFLoc[triVInd[2]].insert(std::pair<int, int>(triI, 2));
        if constexpr (dim == 3) {
            vFLoc[triVInd[3]].insert(std::pair<int, int>(triI, 3));
        }

        const Eigen::Vector3d& P1 = V_rest.row(triVInd[0]);
        const Eigen::Vector3d& P2 = V_rest.row(triVInd[1]);
        const Eigen::Vector3d& P3 = V_rest.row(triVInd[2]);

        Eigen::Matrix<double, dim, dim> X0;
        const Eigen::Vector3d& P4 = V_rest.row(triVInd[3]);
        X0.col(0) = P2 - P1;
        X0.col(1) = P3 - P1;
        X0.col(2) = P4 - P1;
        restTriInv[triI] = X0.inverse();

        const Eigen::Vector3d P2m1 = P2 - P1;
        const Eigen::Vector3d P3m1 = P3 - P1;
        const Eigen::RowVector3d normalVec = P2m1.cross(P3m1).transpose();

        triArea[triI] = X0.determinant() / dim / (dim - 1);
        vertNormals[triVInd[0]] += normalVec;
        vertNormals[triVInd[1]] += normalVec;
        vertNormals[triVInd[2]] += normalVec;
    }
    avgEdgeLen = igl::avg_edge_length(V_rest, F);
    for (auto& vNI : vertNormals) {
        vNI.normalize();
    }

    computeMassMatrix(igl::MASSMATRIX_TYPE_VORONOI);

    bbox.block(0, 0, 1, 3) = V_rest.row(0);
    bbox.block(1, 0, 1, 3) = V_rest.row(0);
    for (int vI = 1; vI < V_rest.rows(); vI++) {
        const Eigen::RowVector3d& v = V_rest.row(vI);
        for (int dimI = 0; dimI < 3; dimI++) {
            if (v[dimI] < bbox(0, dimI)) {
                bbox(0, dimI) = v[dimI];
            }
            if (v[dimI] > bbox(1, dimI)) {
                bbox(1, dimI) = v[dimI];
            }
        }
    }

    vNeighbor.resize(0);
    vNeighbor.resize(V_rest.rows());
    for (int elemI = 0; elemI < F.rows(); elemI++) {
        const Eigen::Matrix<int, 1, dim + 1>& elemVInd = F.row(elemI);
        for (int vI = 0; vI < dim + 1; vI++) {
            for (int vJ = vI + 1; vJ < dim + 1; vJ++) {
                vNeighbor[elemVInd[vI]].insert(elemVInd[vJ]);
                vNeighbor[elemVInd[vJ]].insert(elemVInd[vI]);
            }
        }
    }
    for (int sfI = 0; sfI < SF.rows(); ++sfI) {
        const Eigen::Matrix<int, 1, dim>& sfVInd = SF.row(sfI);
        for (int vI = 0; vI < dim; vI++) {
            for (int vJ = vI + 1; vJ < dim; vJ++) {
                vNeighbor[sfVInd[vI]].insert(sfVInd[vJ]);
                vNeighbor[sfVInd[vJ]].insert(sfVInd[vI]);
            }
        }
    }
    for (int ceI = 0; ceI < CE.rows(); ++ceI) {
        vNeighbor[CE(ceI, 0)].insert(CE(ceI, 1));
        vNeighbor[CE(ceI, 1)].insert(CE(ceI, 0));
    }

    std::set<std::pair<int, int>> SFEdges_set;
    for (int sfI = 0; sfI < SF.rows(); ++sfI) {
        auto finder = SFEdges_set.find(std::pair<int, int>(SF(sfI, 1), SF(sfI, 0)));
        if (finder == SFEdges_set.end()) {
            SFEdges_set.insert(std::pair<int, int>(SF(sfI, 0), SF(sfI, 1)));
        }

        finder = SFEdges_set.find(std::pair<int, int>(SF(sfI, 2), SF(sfI, 1)));
        if (finder == SFEdges_set.end()) {
            SFEdges_set.insert(std::pair<int, int>(SF(sfI, 1), SF(sfI, 2)));
        }

        finder = SFEdges_set.find(std::pair<int, int>(SF(sfI, 0), SF(sfI, 2)));
        if (finder == SFEdges_set.end()) {
            SFEdges_set.insert(std::pair<int, int>(SF(sfI, 2), SF(sfI, 0)));
        }
    }
    SFEdges = std::vector<std::pair<int, int>>(SFEdges_set.begin(), SFEdges_set.end());
    for (int ceI = 0; ceI < CE.rows(); ++ceI) {
        SFEdges.emplace_back(CE(ceI, 0), CE(ceI, 1));
    }

    // sanity check
    int isoNodeCount = 0;
    for (int vI = 0; vI < V_rest.rows(); ++vI) {
        if (vNeighbor[vI].empty()) {
            ++isoNodeCount;
        }
    }
    spdlog::info("{:d} isolated nodes detected", isoNodeCount);

    computeBoundaryVert(SF);
}

template <int dim>
void Mesh<dim>::resetFixedVert(const std::set<int>& p_fixedVert)
{
    isFixedVert.resize(0);
    isFixedVert.resize(V.rows(), false);
    for (const auto& vI : p_fixedVert) {
        assert(vI < V.rows());
        isFixedVert[vI] = true;
    }

    fixedVert = p_fixedVert;
}
template <int dim>
void Mesh<dim>::addFixedVert(int vI)
{
    assert(vI < V.rows());
    fixedVert.insert(vI);
    isFixedVert[vI] = true;
}
template <int dim>
void Mesh<dim>::addFixedVert(const std::vector<int>& p_fixedVert)
{
    for (const auto& vI : p_fixedVert) {
        assert(vI < V.rows());
        isFixedVert[vI] = true;
    }

    fixedVert.insert(p_fixedVert.begin(), p_fixedVert.end());
}

template <int dim>
void Mesh<dim>::removeFixedVert(int vI)
{
    assert(vI < V.rows());
    fixedVert.erase(vI);
    isFixedVert[vI] = false;
}

template <int dim>
double Mesh<dim>::avgNodeMass(int coDim) const
{
    assert(massMatrix.diagonal().size() == V.rows());

    if (coDim == -1) {
        // include all masses
        return Eigen::VectorXd(massMatrix.diagonal()).mean();
    }
    else {
        assert(coDim >= 0 && coDim <= 3);

        double mSum = 0.0;
        int totalVAmt = 0;
        for (int compI = 0; compI < componentCoDim.size(); ++compI) {
            if (componentCoDim[compI] == coDim) {
                mSum += massMatrix.diagonal().segment(componentNodeRange[compI], componentNodeRange[compI + 1] - componentNodeRange[compI]).sum();
                totalVAmt += componentNodeRange[compI + 1] - componentNodeRange[compI];
            }
        }

        if (totalVAmt) {
            return mSum / totalVAmt;
        }
        else {
            return 0.0;
        }
    }
}
template <int dim>
double Mesh<dim>::matSpaceBBoxSize2(int coDim) const
{
    if (coDim == -1) {
        // include all nodes
        return (V_rest.colwise().maxCoeff() - V_rest.colwise().minCoeff()).squaredNorm();
    }
    else {
        assert(coDim >= 0 && coDim <= 3);

        Eigen::Array<double, 1, dim> bottomLeft, topRight;
        bottomLeft = __DBL_MAX__;
        topRight = -__DBL_MAX__;
        for (int compI = 0; compI < componentCoDim.size(); ++compI) {
            if (componentCoDim[compI] == coDim) {
                bottomLeft = bottomLeft.min(V_rest.block(componentNodeRange[compI], 0, componentNodeRange[compI + 1] - componentNodeRange[compI], dim).colwise().minCoeff().array());
                topRight = topRight.max(V_rest.block(componentNodeRange[compI], 0, componentNodeRange[compI + 1] - componentNodeRange[compI], dim).colwise().maxCoeff().array());
            }
        }

        if ((bottomLeft - topRight > 0.0).any()) {
            return 0.0;
        }
        else {
            return (topRight - bottomLeft).matrix().squaredNorm();
        }
    }
}

template <int dim>
int Mesh<dim>::vICoDim(int vI) const
{
    int compI = 1;
    while (compI < componentNodeRange.size() && vI >= componentNodeRange[compI]) {
        ++compI;
    }
    assert(compI < componentNodeRange.size());
    return componentCoDim[compI - 1];
}
template <int dim>
int Mesh<dim>::sfICoDim(int sfI) const
{
    int compI = 1;
    while (compI < componentSFRange.size() && sfI >= componentSFRange[compI]) {
        ++compI;
    }
    assert(compI < componentSFRange.size());
    return componentCoDim[compI - 1];
}

template <int dim>
void Mesh<dim>::setLameParam(double YM, double PR)
{
    u = Eigen::VectorXd::Constant(F.rows(), YM / 2.0 / (1.0 + PR));
    lambda = Eigen::VectorXd::Constant(F.rows(), YM * PR / (1.0 + PR) / (1.0 - 2.0 * PR));
    for (const auto& matI : componentMaterial) {
        massMatrix.diagonal().segment(componentNodeRange[matI.first[0]],
            componentNodeRange[matI.first[0] + 1] - componentNodeRange[matI.first[0]])
            *= matI.second[0] / density;
        u.segment(matI.first[1], matI.first[2] - matI.first[1]).setConstant(matI.second[1] / 2.0 / (1.0 + matI.second[2]));
        lambda.segment(matI.first[1], matI.first[2] - matI.first[1]).setConstant(matI.second[1] * matI.second[2] / (1.0 + matI.second[2]) / (1.0 - 2.0 * matI.second[2]));
    }

    //        int nSeg = 8;
    //        double xMin = V_rest.col(0).minCoeff();
    //        double xSegLen = (V_rest.col(0).maxCoeff() - xMin) / nSeg;
    //        YM *= 1000, PR = 0.49;
    //        double u_stiffer = YM / 2.0 / (1.0 + PR);
    //        double lambda_stiffer = YM * PR / (1.0 + PR) / (1.0 - 2.0 * PR);
    //        for(int elemI = 0; elemI < F.rows(); ++elemI) {
    //            const Eigen::Matrix<int, 1, dim + 1>& elemVInd = F.row(elemI);
    //            Eigen::RowVector3d center;
    //            center.setZero();
    //            for(int i = 0; i < dim + 1; ++i) {
    //                center += V_rest.row(elemVInd[i]);
    //            }
    //            center /= dim + 1;
    //
    //            int segI = std::min(nSeg - 1, std::max(0, int((center[0] - xMin) / xSegLen)));
    //            if(segI % 2 == 0) {
    //                u[elemI] = u_stiffer;
    //                lambda[elemI] = lambda_stiffer;
    //            }
    //        }

    //        double yThres = (V_rest.col(1).maxCoeff() + V_rest.col(1).minCoeff()) / 2.0;
    //        YM *= 10;
    //        double u_stiffer = YM / 2.0 / (1.0 + PR);
    //        double lambda_stiffer = YM * PR / (1.0 + PR) / (1.0 - 2.0 * PR);
    //        for(int elemI = 0; elemI < F.rows(); ++elemI) {
    //            const Eigen::Matrix<int, 1, dim + 1>& elemVInd = F.row(elemI);
    //            Eigen::RowVector3d center;
    //            center.setZero();
    //            for(int i = 0; i < dim + 1; ++i) {
    //                center += V_rest.row(elemVInd[i]);
    //            }
    //            center /= dim + 1;
    //
    //            if(center[1] > yThres) {
    //                u[elemI] = u_stiffer;
    //                lambda[elemI] = lambda_stiffer;
    //            }
    //        }
}

template <int dim>
bool Mesh<dim>::checkInversion(int triI, bool mute) const
{
    assert(triI < F.rows());

    const double eps = 0.0; //1.0e-20 * avgEdgeLen * avgEdgeLen;

    const Eigen::Matrix<int, 1, dim + 1>& triVInd = F.row(triI);

    Eigen::Matrix<double, dim, dim> e_u;
    e_u.col(0) = (V.row(triVInd[1]) - V.row(triVInd[0])).transpose();
    e_u.col(1) = (V.row(triVInd[2]) - V.row(triVInd[0])).transpose();
    if constexpr (dim == 3) {
        e_u.col(2) = (V.row(triVInd[3]) - V.row(triVInd[0])).transpose();
    }

    const double dbArea = e_u.determinant();
    if (dbArea < eps) {
        if (!mute) {
            spdlog::error("Element inversion detected: {:g} < {:g}", dbArea, eps);
            spdlog::info("mesh triangle count: ", F.rows());
            logFile << "***Element inversion detected: " << dbArea << " < " << eps << std::endl;
        }
        return false;
    }
    else {
        return true;
    }
}
template <int dim>
bool Mesh<dim>::checkInversion(bool mute, const std::vector<int>& triangles) const
{
    if (triangles.empty()) {
        for (int triI = 0; triI < F.rows(); triI++) {
            if (u[triI] && lambda[triI] && !checkInversion(triI, mute)) {
                return false;
            }
        }
    }
    else {
        for (const auto& triI : triangles) {
            if (u[triI] && lambda[triI] && !checkInversion(triI, mute)) {
                return false;
            }
        }
    }

    return true;
}

template <int dim>
void Mesh<dim>::save(const std::string& filePath, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
    const Eigen::MatrixXd UV, const Eigen::MatrixXi& FUV) const
{
    assert(dim == 2);
#if (DIM == 2)
    std::ofstream out;
    out.open(filePath);
    assert(out.is_open());

    for (int vI = 0; vI < V.rows(); vI++) {
        const Eigen::RowVector3d& v = V.row(vI);
        out << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    for (int vI = 0; vI < UV.rows(); vI++) {
        const Eigen::RowVector2d& uv = UV.row(vI);
        out << "vt " << uv[0] << " " << uv[1] << std::endl;
    }

    if (FUV.rows() == F.rows()) {
        for (int triI = 0; triI < F.rows(); triI++) {
            const Eigen::RowVector3i& tri = F.row(triI);
            const Eigen::RowVector3i& tri_UV = FUV.row(triI);
            out << "f " << tri[0] + 1 << "/" << tri_UV[0] + 1 << " " << tri[1] + 1 << "/" << tri_UV[1] + 1 << " " << tri[2] + 1 << "/" << tri_UV[2] + 1 << std::endl;
        }
    }
    else {
        for (int triI = 0; triI < F.rows(); triI++) {
            const Eigen::RowVector3i& tri = F.row(triI);
            out << "f " << tri[0] + 1 << "/" << tri[0] + 1 << " " << tri[1] + 1 << "/" << tri[1] + 1 << " " << tri[2] + 1 << "/" << tri[2] + 1 << std::endl;
        }
    }

    out.close();
#endif
}

template <int dim>
void Mesh<dim>::save(const std::string& filePath) const
{
    save(filePath, V_rest, F, V);
}

template <int dim>
void Mesh<dim>::saveAsMesh(const std::string& filePath, bool scaleUV,
    const Eigen::MatrixXi& SF) const
{
    IglUtils::saveTetMesh(filePath, V, F, SF);
}

template <int dim>
void Mesh<dim>::saveSurfaceMesh(const std::string& filePath, bool useInvShift) const
{
    if constexpr (dim == 3) {
        Eigen::MatrixXd V_surf(SVI.size(), 3);
        std::unordered_map<int, int> vI2SVI;
        for (int svI = 0; svI < SVI.size(); ++svI) {
            vI2SVI[SVI[svI]] = svI;
            V_surf.row(svI) = V.row(SVI[svI]);
        }

        Eigen::MatrixXi F_surf(SF.rows(), 3);
        for (int sfI = 0; sfI < SF.rows(); ++sfI) {
            F_surf(sfI, 0) = vI2SVI[SF(sfI, 0)];
            F_surf(sfI, 1) = vI2SVI[SF(sfI, 1)];
            F_surf(sfI, 2) = vI2SVI[SF(sfI, 2)];
        }

        if (useInvShift) {
            V_surf.rowwise() += invShift.transpose();
        }
        igl::writeOBJ(filePath, V_surf, F_surf);
    }
    else {
        saveAsMesh(filePath);
    }
}

template <int dim>
void Mesh<dim>::saveBCNodes(const std::string& filePath) const
{
    FILE* out = fopen(filePath.c_str(), "w");
    assert(out);
    for (const auto& fixedVI : fixedVert) {
        fprintf(out, "%d\n", fixedVI);
    }
    for (const auto& NMI : NeumannBC) {
        fprintf(out, "%d %le %le", NMI.first, NMI.second[0], NMI.second[1]);
        if constexpr (dim == 2) {
            fprintf(out, "\n");
        }
        else {
            fprintf(out, " %le\n", NMI.second[2]);
        }
    }
    fclose(out);
}

template <int dim>
bool Mesh<dim>::isBoundaryVert(int vI) const
{
    assert(vI < m_isBoundaryVert.size());
    return m_isBoundaryVert[vI];
}

template <int dim>
void Mesh<dim>::computeBoundaryVert(const Eigen::MatrixXi& SF)
{
#if (DIM == 2)
    m_isBoundaryVert.resize(0);
    m_isBoundaryVert.resize(V.rows(), false);
    std::vector<std::vector<int>> Ls;
    igl::boundary_loop(F, Ls);
    for (auto& L : Ls)
        for (auto& vI : L)
            m_isBoundaryVert[vI] = true;
#else
    m_isBoundaryVert.resize(0);
    m_isBoundaryVert.resize(V.rows(), false);
    for (int sfI = 0; sfI < SF.rows(); ++sfI) {
        const Eigen::Matrix<int, 1, dim>& sfVInd = SF.row(sfI);
        m_isBoundaryVert[sfVInd[0]] = true;
        m_isBoundaryVert[sfVInd[1]] = true;
        if constexpr (dim == 3) {
            m_isBoundaryVert[sfVInd[2]] = true;
        }
    }
    for (int ceI = 0; ceI < CE.rows(); ++ceI) {
        m_isBoundaryVert[CE(ceI, 0)] = true;
        m_isBoundaryVert[CE(ceI, 1)] = true;
    }
    for (int vI = 0; vI < vNeighbor.size(); ++vI) {
        if (vNeighbor[vI].empty()) {
            m_isBoundaryVert[vI] = true;
        }
    }

    SVI.resize(0);
    for (int vI = 0; vI < V.rows(); ++vI) {
        if (m_isBoundaryVert[vI]) {
            SVI.conservativeResize(SVI.size() + 1);
            SVI[SVI.size() - 1] = vI;
        }
    }
#endif
}

template class Mesh<DIM>;

} // namespace IPC
