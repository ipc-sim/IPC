//
//  IglUtils.cpp
//  IPC
//
//  Created by Minchen Li on 8/30/17.
//

#include "IglUtils.hpp"
#include "Triplet.hpp"

#include "Timer.hpp"

#ifdef USE_TBB
#include <tbb/parallel_for.h>
#endif

#include <set>
#include <unordered_map>

#include <ghc/fs_std.hpp> // filesystem
#include <spdlog/spdlog.h>
#include <mshio/mshio.h>

extern Timer timer_temp, timer_temp2;

namespace IPC {

void IglUtils::writeSparseMatrixToFile(const std::string& filePath, const Eigen::SparseMatrix<double>& mtr, bool MATLAB)
{
    std::ofstream out;
    out.open(filePath);
    if (out.is_open()) {
        if (!MATLAB) {
            out << mtr.rows() << " " << mtr.cols() << " " << mtr.nonZeros() << std::endl;
        }
        for (int k = 0; k < mtr.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(mtr, k); it; ++it) {
                out << it.row() + MATLAB << " " << it.col() + MATLAB << " " << it.value() << std::endl;
            }
        }
        out.close();
    }
    else {
        spdlog::error("writeSparseMatrixToFile failed! file open error!");
    }
}

void IglUtils::writeSparseMatrixToFile(const std::string& filePath,
    const std::map<std::pair<int, int>, double>& mtr,
    bool MATLAB)
{
    std::ofstream out;
    out.open(filePath);
    if (out.is_open()) {
        if (!MATLAB) {
            int mtrSize = mtr.rbegin()->first.first + 1;
            out << mtrSize << " " << mtrSize << " " << mtr.size() << std::endl;
        }
        for (const auto& entryI : mtr) {
            out << entryI.first.first + MATLAB << " "
                << entryI.first.second + MATLAB << " "
                << entryI.second << std::endl;
        }
        out.close();
    }
    else {
        spdlog::error("writeSparseMatrixToFile failed! file open error!");
    }
}

void IglUtils::writeSparseMatrixToFile(const std::string& filePath,
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver,
    bool MATLAB)
{
    std::ofstream out;
    out.open(filePath);
    if (out.is_open()) {
        if (!MATLAB) {
            out << linSysSolver->getNumRows() << " " << linSysSolver->getNumRows() << " " << linSysSolver->getNumNonzeros() << std::endl;
        }
        for (int rowI = 0; rowI < linSysSolver->getNumRows(); rowI++) {
            for (const auto& colIter : linSysSolver->getIJ2aI()[rowI]) {
                out << rowI + MATLAB << " "
                    << colIter.first + MATLAB << " "
                    << linSysSolver->coeffMtr(rowI, colIter.first) << std::endl;
            }
        }
        out.close();
    }
    else {
        spdlog::error("writeSparseMatrixToFile failed! file open error!");
    }
}

void IglUtils::loadSparseMatrixFromFile(const std::string& filePath, Eigen::SparseMatrix<double>& mtr)
{
    std::ifstream in;
    in.open(filePath);
    if (in.is_open()) {
        int rows, cols, nonZeroAmt;
        in >> rows >> cols >> nonZeroAmt;
        mtr.resize(rows, cols);
        std::vector<Eigen::Triplet<double>> IJV;
        IJV.reserve(nonZeroAmt);
        int i, j;
        double v;
        for (int nzI = 0; nzI < nonZeroAmt; nzI++) {
            assert(!in.eof());
            in >> i >> j >> v;
            IJV.emplace_back(Eigen::Triplet<double>(i, j, v));
        }
        in.close();
        mtr.setFromTriplets(IJV.begin(), IJV.end());
    }
    else {
        spdlog::error("loadSparseMatrixToFile failed! file open error!");
    }
}

void IglUtils::writeVectorToFile(const std::string& filePath,
    const Eigen::VectorXd& vec)
{
    FILE* out = fopen(filePath.c_str(), "w");
    assert(out);
    for (int dI = 0; dI < vec.size(); ++dI) {
        fprintf(out, "%le\n", vec[dI]);
    }
    fclose(out);
}
void IglUtils::readVectorFromFile(const std::string& filePath,
    Eigen::VectorXd& vec)
{
    FILE* in = fopen(filePath.c_str(), "r");
    assert(in);

    vec.resize(0);
    char buf[BUFSIZ];
    while (!feof(in) && fgets(buf, BUFSIZ, in)) {
        vec.conservativeResize(vec.size() + 1);
        sscanf(buf, "%le", &vec[vec.size() - 1]);
    }

    fclose(in);
}

bool IglUtils::readSEG(const std::string& filePath, Eigen::MatrixXd& V, Eigen::MatrixXi& E)
{
    std::ifstream file(filePath);
    if (file.is_open()) {
        V.resize(0, 3);
        E.resize(0, 2);

        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string token;
            ss >> token;
            if (token == "v") {
                V.conservativeResize(V.rows() + 1, V.cols());
                ss >> V(V.rows() - 1, 0) >> V(V.rows() - 1, 1) >> V(V.rows() - 1, 2);
            }
            else if (token == "s") {
                E.conservativeResize(E.rows() + 1, E.cols());
                ss >> E(E.rows() - 1, 0) >> E(E.rows() - 1, 1);
            }
        }
        file.close();
        E.array() -= 1;
        return true;
    }
    else {
        spdlog::warn("SEG file does not exist");
        return false;
    }
}

void IglUtils::writeSEG(const std::string& filePath,
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& E)
{
    FILE* file = fopen(filePath.c_str(), "w");
    assert(file);

    for (int vI = 0; vI < V.rows(); ++vI) {
        fprintf(file, "v %.10le %.10le %.10le\n", V(vI, 0), V(vI, 1), V(vI, 2));
    }

    for (int ceI = 0; ceI < E.rows(); ++ceI) {
        fprintf(file, "s %d %d\n", E(ceI, 0) + 1, E(ceI, 1) + 1);
    }

    fclose(file);
}

/////////////////////////////////////////////////////////////////
// 2D line segments intersection checking code
// based on Real-Time Collision Detection by Christer Ericson
// (Morgan Kaufmaan Publishers, 2005 Elvesier Inc)
double Signed2DTriArea(const Eigen::RowVector2d& a, const Eigen::RowVector2d& b, const Eigen::RowVector2d& c)
{
    return (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0]);
}

void IglUtils::findSurfaceTris(const Eigen::MatrixXi& TT, Eigen::MatrixXi& F)
{
    // TODO: merge with below
    std::unordered_map<Triplet, int> tri2Tet(4 * TT.rows());
    for (int elemI = 0; elemI < TT.rows(); elemI++) {
        const Eigen::RowVector4i& elemVInd = TT.row(elemI);
        tri2Tet[Triplet(elemVInd[0], elemVInd[2], elemVInd[1])] = elemI;
        tri2Tet[Triplet(elemVInd[0], elemVInd[3], elemVInd[2])] = elemI;
        tri2Tet[Triplet(elemVInd[0], elemVInd[1], elemVInd[3])] = elemI;
        tri2Tet[Triplet(elemVInd[1], elemVInd[2], elemVInd[3])] = elemI;
    }

    // TODO: parallelize
    std::vector<Eigen::RowVector3i> tmpF;
    for (const auto& triI : tri2Tet) {
        const Triplet& triVInd = triI.first;
        // find dual triangle with reversed indices:
        bool isSurfaceTriangle = //
            tri2Tet.find(Triplet(triVInd[2], triVInd[1], triVInd[0])) == tri2Tet.end()
            && tri2Tet.find(Triplet(triVInd[1], triVInd[0], triVInd[2])) == tri2Tet.end()
            && tri2Tet.find(Triplet(triVInd[0], triVInd[2], triVInd[1])) == tri2Tet.end();
        if (isSurfaceTriangle) {
            tmpF.emplace_back(triVInd[0], triVInd[1], triVInd[2]);
        }
    }

    F.resize(tmpF.size(), 3);
    for (int i = 0; i < F.rows(); i++) {
        F.row(i) = tmpF[i];
    }
}
void IglUtils::buildSTri2Tet(const Eigen::MatrixXi& F, const Eigen::MatrixXi& SF,
    std::vector<int>& sTri2Tet)
{
    // TODO: merge with above
    std::map<Triplet, int> tri2Tet;
    for (int elemI = 0; elemI < F.rows(); elemI++) {
        const Eigen::RowVector4i& elemVInd = F.row(elemI);
        tri2Tet[Triplet(elemVInd[0], elemVInd[2], elemVInd[1])] = elemI;
        tri2Tet[Triplet(elemVInd[0], elemVInd[3], elemVInd[2])] = elemI;
        tri2Tet[Triplet(elemVInd[0], elemVInd[1], elemVInd[3])] = elemI;
        tri2Tet[Triplet(elemVInd[1], elemVInd[2], elemVInd[3])] = elemI;
    }

    sTri2Tet.resize(SF.rows());
#ifdef USE_TBB
    tbb::parallel_for(0, (int)SF.rows(), 1, [&](int triI)
#else
    for (int triI = 0; triI < SF.rows(); triI++)
#endif
        {
            const Eigen::RowVector3i& triVInd = SF.row(triI);
            auto finder = tri2Tet.find(Triplet(triVInd.data()));
            if (finder == tri2Tet.end()) {
                finder = tri2Tet.find(Triplet(triVInd[1], triVInd[2], triVInd[0]));
                if (finder == tri2Tet.end()) {
                    finder = tri2Tet.find(Triplet(triVInd[2], triVInd[0], triVInd[1]));
                }
            }

            if (finder == tri2Tet.end()) {
                sTri2Tet[triI] = -1; // surface mesh
            }
            else {
                sTri2Tet[triI] = finder->second;
            }
        }
#ifdef USE_TBB
    );
#endif
}
void IglUtils::saveTetMesh(const std::string& filePath,
    const Eigen::MatrixXd& TV, const Eigen::MatrixXi& TT,
    const Eigen::MatrixXi& p_F, bool findSurface)
{
    // saveTetMesh_msh4(filePath, TV, TT, p_F, findSurface);

    assert(TV.rows() > 0);
    assert(TV.cols() == 3);
    assert(TT.rows() > 0);
    assert(TT.cols() == 4);

    mshio::MshSpec spec;
    auto& format = spec.mesh_format;
    format.version = "4.1"; // Only version "2.2" and "4.1" are supported.
    format.file_type = 0; // 0: ASCII, 1: binary.
    format.data_size = sizeof(size_t); // Size of data, defined as sizeof(size_t) = 8.

    // Store the nodes
    auto& nodes = spec.nodes;
    nodes.num_entity_blocks = 1; // Number of node blocks.
    nodes.num_nodes = TV.rows(); // Total number of nodes.
    nodes.min_node_tag = 1;
    nodes.max_node_tag = TV.rows();
    nodes.entity_blocks.resize(1); // A std::vector of node blocks.

    auto& node_block = nodes.entity_blocks.front();
    node_block.entity_dim = 3; // The dimension of the entity.
    node_block.entity_tag = 0; // The entity these nodes belongs to.
    node_block.parametric = 0; // 0: non-parametric, 1: parametric.
    node_block.num_nodes_in_block = TV.rows(); // The number of nodes in block.
    node_block.tags.resize(TV.rows()); // A std::vector of unique, positive node tags.
    node_block.data.resize(TV.size()); // A std::vector of coordinates (x,y,z,<u>,<v>,<w>,...)
    for (size_t i = 0; i < TV.rows(); i++) {
        node_block.tags[i] = i + 1;
        for (size_t j = 0; j < TV.cols(); j++) {
            node_block.data[TV.cols() * i + j] = TV(i, j);
        }
    }

    // Store the elements
    auto& elements = spec.elements;
    elements.num_entity_blocks = 1; // Number of element blocks.
    elements.num_elements = TT.rows(); // Total number of elements.
    elements.min_element_tag = 1;
    elements.max_element_tag = TT.rows();
    elements.entity_blocks.resize(1); // A std::vector of element blocks.

    auto& elm_block = elements.entity_blocks.front();
    elm_block.entity_dim = 3; // The dimension of the elements.
    elm_block.entity_tag = 0; // The entity these elements belongs to.
    elm_block.element_type = 4; // Linear tet.
    elm_block.num_elements_in_block = TT.rows(); // The number of elements in this block.
    elm_block.data.resize(TT.size() + TT.rows(), 0);
    size_t data_index = 0;
    for (size_t i = 0; i < TT.rows(); i++) {
        elm_block.data[data_index++] = i + 1;
        for (size_t j = 0; j < TT.cols(); j++) {
            elm_block.data[data_index++] = TT(i, j) + 1;
        }
    }

    Eigen::MatrixXi F_found;
    if (p_F.rows() > 0) {
        assert(p_F.cols() == 3);
    }
    else if (findSurface) {
        findSurfaceTris(TT, F_found);
    }
    const Eigen::MatrixXi& F = ((p_F.rows() > 0) ? p_F : F_found);

    std::ofstream out(filePath, std::ios::out);
    if (out.is_open()) {
        // Save the msh file
        mshio::save_msh(out, spec);
        // Append the file with the surface
        out << "$Surface\n";
        out << F.rows() << "\n";
        for (int triI = 0; triI < F.rows(); triI++) {
            const Eigen::RowVector3i& triVInd = F.row(triI);
            out << (triVInd[0] + 1) << " " << (triVInd[1] + 1) << " " << (triVInd[2] + 1) << "\n";
        }
        out << "$EndSurface\n";
        out.close();
    }
    else {
        spdlog::error("Unable to save mesh to {}", filePath);
    }
}
void IglUtils::saveTetMesh_msh4(const std::string& filePath,
    const Eigen::MatrixXd& TV, const Eigen::MatrixXi& TT,
    const Eigen::MatrixXi& p_F, bool findSurface)
{
    assert(TV.rows() > 0);
    assert(TV.cols() == 3);
    assert(TT.rows() > 0);
    assert(TT.cols() == 4);

    Eigen::MatrixXi F_found;
    if (p_F.rows() > 0) {
        assert(p_F.cols() == 3);
    }
    else if (findSurface) {
        findSurfaceTris(TT, F_found);
    }
    const Eigen::MatrixXi& F = ((p_F.rows() > 0) ? p_F : F_found);

    FILE* out = fopen(filePath.c_str(), "w");
    assert(out);

    fprintf(out, "$MeshFormat\n4 0 8\n$EndMeshFormat\n");

    fprintf(out, "$Entities\n0 0 0 1\n");
    fprintf(out, "0 %le %le %le %le %le %le 0 0\n$EndEntities\n",
        TV.col(0).minCoeff(), TV.col(1).minCoeff(), TV.col(2).minCoeff(),
        TV.col(0).maxCoeff(), TV.col(1).maxCoeff(), TV.col(2).maxCoeff());

    fprintf(out, "$Nodes\n1 %lu\n0 3 0 %lu\n", TV.rows(), TV.rows());
    for (int vI = 0; vI < TV.rows(); vI++) {
        const Eigen::RowVector3d& v = TV.row(vI);
        fprintf(out, "%d %le %le %le\n", vI + 1, v[0], v[1], v[2]);
    }
    fprintf(out, "$EndNodes\n");

    fprintf(out, "$Elements\n1 %lu\n0 3 4 %lu\n", TT.rows(), TT.rows());
    for (int elemI = 0; elemI < TT.rows(); elemI++) {
        const Eigen::RowVector4i& tetVInd = TT.row(elemI);
        fprintf(out, "%d %d %d %d %d\n", elemI + 1,
            tetVInd[0] + 1, tetVInd[1] + 1, tetVInd[2] + 1, tetVInd[3] + 1);
    }
    fprintf(out, "$EndElements\n");

    fprintf(out, "$Surface\n");
    fprintf(out, "%lu\n", F.rows());
    for (int triI = 0; triI < F.rows(); triI++) {
        const Eigen::RowVector3i& triVInd = F.row(triI);
        fprintf(out, "%d %d %d\n", triVInd[0] + 1, triVInd[1] + 1, triVInd[2] + 1);
    }
    fprintf(out, "$EndSurface\n");

    fclose(out);
}
void IglUtils::saveTetMesh_vtk(const std::string& filePath,
    const Eigen::MatrixXd& X, const Eigen::MatrixXi& indices)
{
    std::ofstream os(filePath);

    assert(X.rows() != 0 && "The X array for writing tetmesh vtk is empty.");
    assert(indices.rows() != (size_t)0 && "The tet mesh data structure for writing tetmesh vtk is empty.");

    os << "# vtk DataFile Version 2.0\n";
    os << "Unstructured Grid\n";
    os << "ASCII\n";
    os << "DATASET UNSTRUCTURED_GRID\n";

    os << "POINTS " << X.rows() << " ";
    os << "double\n";

    for (size_t i = 0; i < X.rows(); i++) {
        os << X(i, 0) << " " << X(i, 1) << " " << X(i, 2) << "\n";
    }

    os << std::endl;

    os << "CELLS " << indices.rows() << " " << 5 * indices.rows() << "\n";
    for (int i = 0; i < indices.rows(); ++i) {
        const Eigen::RowVector4i& m = indices.row(i);
        os << 4 << " " << m(0) << " " << m(1) << " " << m(2) << " " << m(3) << "\n";
    }
    os << std::endl;

    os << "CELL_TYPES " << indices.rows() << "\n";
    for (size_t i = 0; i < indices.rows(); i++) {
        os << 10 << std::endl;
    }

    os.close();
}
bool IglUtils::readTetMesh(const std::string& filePath,
    Eigen::MatrixXd& TV, Eigen::MatrixXi& TT,
    Eigen::MatrixXi& F, bool findSurface)
{
    if (!fs::exists(filePath)) {
        return false;
    }

    mshio::MshSpec spec;
    try {
        spec = mshio::load_msh(filePath);
    }
    catch (...) {
        // MshIO only supports MSH 2.2 and 4.1 not 4.0
        return readTetMesh_msh4(filePath, TV, TT, F, findSurface);
    }

    const auto& nodes = spec.nodes;
    const auto& els = spec.elements;
    const int vAmt = nodes.num_nodes;
    int elemAmt = 0;
    for (const auto& e : els.entity_blocks) {
        assert(e.entity_dim == 3);
        assert(e.element_type == 4); // linear tet
        elemAmt += e.num_elements_in_block;
    }

    TV.resize(vAmt, 3);
    int index = 0;
    for (const auto& n : nodes.entity_blocks) {
        for (int i = 0; i < n.num_nodes_in_block * 3; i += 3) {
            TV.row(index) << n.data[i], n.data[i + 1], n.data[i + 2];
            ++index;
        }
    }

    TT.resize(elemAmt, 4);
    int elm_index = 0;
    for (const auto& e : els.entity_blocks) {
        for (int i = 0; i < e.data.size(); i += 5) {
            index = 0;
            for (int j = i + 1; j <= i + 4; ++j) {
                TT(elm_index, index++) = e.data[j] - 1;
            }
            ++elm_index;
        }
    }

    if (!findSurface) {
        spdlog::warn("readTetMesh is finding the surface because $Surface is not supported by MshIO");
    }
    spdlog::info("Finding the surface triangle mesh for {:s}", filePath);
    findSurfaceTris(TT, F);

    spdlog::info(
        "tet mesh loaded with {:d} nodes, {:d} tets, and {:d} surface triangles.",
        TV.rows(), TT.rows(), F.rows());

    return true;
}
bool IglUtils::readTetMesh_msh4(const std::string& filePath,
    Eigen::MatrixXd& TV, Eigen::MatrixXi& TT,
    Eigen::MatrixXi& F, bool findSurface)
{
    FILE* in = fopen(filePath.c_str(), "r");
    if (!in) {
        return false;
    }

    TV.resize(0, 3);
    TT.resize(0, 4);
    F.resize(0, 3);

    char buf[BUFSIZ];
    while ((!feof(in)) && fgets(buf, BUFSIZ, in)) {
        if (strncmp("$Nodes", buf, 6) == 0) {
            fgets(buf, BUFSIZ, in);
            int vAmt;
            sscanf(buf, "1 %d", &vAmt);
            TV.resize(vAmt, 3);
            fgets(buf, BUFSIZ, in);
            break;
        }
    }
    assert(TV.rows() > 0);
    int bypass;
    for (int vI = 0; vI < TV.rows(); vI++) {
        fscanf(in, "%d %le %le %le\n", &bypass, &TV(vI, 0), &TV(vI, 1), &TV(vI, 2));
    }

    while ((!feof(in)) && fgets(buf, BUFSIZ, in)) {
        if (strncmp("$Elements", buf, 9) == 0) {
            fgets(buf, BUFSIZ, in);
            int elemAmt;
            sscanf(buf, "1 %d", &elemAmt);
            TT.resize(elemAmt, 4);
            fgets(buf, BUFSIZ, in);
            break;
        }
    }
    assert(TT.rows() > 0);
    for (int elemI = 0; elemI < TT.rows(); elemI++) {
        fscanf(in, "%d %d %d %d %d\n", &bypass,
            &TT(elemI, 0), &TT(elemI, 1), &TT(elemI, 2), &TT(elemI, 3));
    }
    TT.array() -= 1;

    while ((!feof(in)) && fgets(buf, BUFSIZ, in)) {
        if (strncmp("$Surface", buf, 7) == 0) {
            fgets(buf, BUFSIZ, in);
            int elemAmt;
            sscanf(buf, "%d", &elemAmt);
            F.resize(elemAmt, 3);
            break;
        }
    }
    for (int triI = 0; triI < F.rows(); triI++) {
        fscanf(in, "%d %d %d\n", &F(triI, 0), &F(triI, 1), &F(triI, 2));
    }
    if (F.rows() > 0) {
        F.array() -= 1;
    }
    else if (findSurface) {
        // if no surface triangles information provided, then find
        spdlog::info("Finding the surface triangle mesh for {:s}", filePath);
        findSurfaceTris(TT, F);
    }

    spdlog::info(
        "tet mesh loaded with {:d} nodes, {:d} tets, and {:d} surface triangles.",
        TV.rows(), TT.rows(), F.rows());

    fclose(in);

    return true;
}
void IglUtils::readNodeEle(const std::string& filePath,
    Eigen::MatrixXd& TV, Eigen::MatrixXi& TT,
    Eigen::MatrixXi& F)
{
    FILE* in = fopen((filePath + ".node").c_str(), "r");
    if (!in) {
        spdlog::error("Unable to open node file: {:s}", filePath + ".node");
    }
    assert(in);

    int nN, nDim;
    fscanf(in, "%d %d 0 0", &nN, &nDim);
    spdlog::info("{:d} {:d}", nN, nDim);
    assert(nN >= 4);
    assert(nDim == 3);

    int bypass;
    TV.conservativeResize(nN, nDim);
    for (int vI = 0; vI < nN; ++vI) {
        fscanf(in, "%d %le %le %le", &bypass,
            &TV(vI, 0), &TV(vI, 1), &TV(vI, 2));
    }

    fclose(in);

    in = fopen((filePath + ".ele").c_str(), "r");
    assert(in);

    int nE, nDimp1;
    fscanf(in, "%d %d 0", &nE, &nDimp1);
    spdlog::info("{:d} {:d}", nE, nDimp1);
    assert(nE >= 0);
    assert(nDimp1 == 4);

    TT.conservativeResize(nE, nDimp1);
    for (int tI = 0; tI < nE; ++tI) {
        fscanf(in, "%d %d %d %d %d", &bypass,
            &TT(tI, 0), &TT(tI, 1), &TT(tI, 2), &TT(tI, 3));
    }

    fclose(in);

    findSurfaceTris(TT, F);

    std::cout << "tet mesh loaded with " << TV.rows() << " nodes, "
              << TT.rows() << " tets, and " << F.rows() << " surface tris." << std::endl;
}

void IglUtils::dF_div_dx_mult(const Eigen::Matrix<double, DIM, DIM>& right,
    const Eigen::Matrix<double, DIM, DIM>& A,
    Eigen::Matrix<double, DIM*(DIM + 1), 1>& result)
{
#if (DIM == 2)
    const double _0000 = right(0, 0) * A(0, 0);
    const double _0010 = right(0, 0) * A(1, 0);
    const double _0101 = right(0, 1) * A(0, 1);
    const double _0111 = right(0, 1) * A(1, 1);
    const double _1000 = right(1, 0) * A(0, 0);
    const double _1010 = right(1, 0) * A(1, 0);
    const double _1101 = right(1, 1) * A(0, 1);
    const double _1111 = right(1, 1) * A(1, 1);

    result[2] = _0000 + _0101;
    result[3] = _1000 + _1101;
    result[4] = _0010 + _0111;
    result[5] = _1010 + _1111;
    result[0] = -result[2] - result[4];
    result[1] = -result[3] - result[5];
#else
    result[3] = A.row(0).dot(right.row(0));
    result[4] = A.row(0).dot(right.row(1));
    result[5] = A.row(0).dot(right.row(2));
    result[6] = A.row(1).dot(right.row(0));
    result[7] = A.row(1).dot(right.row(1));
    result[8] = A.row(1).dot(right.row(2));
    result[9] = A.row(2).dot(right.row(0));
    result[10] = A.row(2).dot(right.row(1));
    result[11] = A.row(2).dot(right.row(2));
    result[0] = -result[3] - result[6] - result[9];
    result[1] = -result[4] - result[7] - result[10];
    result[2] = -result[5] - result[8] - result[11];
#endif
}

void IglUtils::findBorderVerts(const Eigen::MatrixXd& V,
    std::vector<std::vector<int>>& borderVerts,
    double ratio)
{
    // resize to match size
    Eigen::RowVectorXd bottomLeft = V.colwise().minCoeff();
    Eigen::RowVectorXd topRight = V.colwise().maxCoeff();
    Eigen::RowVectorXd range = topRight - bottomLeft;

    borderVerts.resize(2);
    for (int vI = 0; vI < V.rows(); vI++) {
        if (V(vI, 0) < bottomLeft[0] + range[0] * ratio) {
            borderVerts[0].emplace_back(vI);
        }
        else if (V(vI, 0) > topRight[0] - range[0] * ratio) {
            borderVerts[1].emplace_back(vI);
        }
    }
}

void IglUtils::Init_Dirichlet(
    const Eigen::MatrixXd& X,
    const Eigen::Vector3d& relBoxMin,
    const Eigen::Vector3d& relBoxMax,
    const std::vector<bool>& isNodeOnBoundary,
    std::vector<int>& selectedVerts)
{
    if (!X.rows()) {
        return;
    }

    Eigen::Array3d bboxMin = X.colwise().minCoeff().array();
    Eigen::Array3d bboxMax = X.colwise().maxCoeff().array();

    Eigen::Array3d rangeMin = (bboxMax - bboxMin) * relBoxMin.array() + bboxMin;
    Eigen::Array3d rangeMax = (bboxMax - bboxMin) * relBoxMax.array() + bboxMin;

    for (int id = 0; id < X.rows(); ++id) {
        if (!isNodeOnBoundary[id]) {
            continue;
        }

        const Eigen::Array3d x = X.row(id).array();
        if ((x >= rangeMin).all() && (x <= rangeMax).all()) {
            selectedVerts.emplace_back(id);
        }
    }
}

} // namespace IPC
