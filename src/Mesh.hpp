//
//  Mesh.hpp
//  IPC
//
//  Created by Minchen Li on 8/30/17.
//

#ifndef Mesh_hpp
#define Mesh_hpp

#include <Eigen/Eigen>

#include <igl/massmatrix.h>

#include <set>
#include <array>
#include <unordered_map>

#include "Types.hpp"

namespace IPC {

struct DirichletBC {
    DirichletBC(const std::vector<int>& vertIds,
        const Eigen::Vector3d& linearVelocity,
        const Eigen::Vector3d& angularVelocity,
        const std::array<double, 2>& timeRange = { 0.0, std::numeric_limits<double>::infinity() })
        : vertIds(vertIds), linearVelocity(linearVelocity), angularVelocity(angularVelocity), timeRange(timeRange) {}

    bool isZero() const
    {
        return linearVelocity.isZero() && angularVelocity.isZero();
    }

    std::vector<int> vertIds;
    Eigen::Vector3d linearVelocity;
    Eigen::Vector3d angularVelocity;
    std::array<double, 2> timeRange = { 0.0, std::numeric_limits<double>::infinity() };
};

enum class DirichletBCType {
    NOT_DBC = 0,
    ZERO = 1,
    NONZERO = 2
};

struct NeumannBC {
    NeumannBC(const std::vector<int>& vertIds,
        const Eigen::Vector3d& force,
        const std::array<double, 2>& timeRange = { 0.0, std::numeric_limits<double>::infinity() })
        : vertIds(vertIds), force(force), timeRange(timeRange) {}

    std::vector<int> vertIds;
    Eigen::Vector3d force;
    std::array<double, 2> timeRange = { 0.0, std::numeric_limits<double>::infinity() };
};

template <int dim>
class Mesh {
public: // owned data
    Eigen::MatrixXd V_rest; ///< @brief node positions in material space (rest shape)
    Eigen::MatrixXd V; ///< @brief node positions in current world space (current shape)
    Eigen::MatrixXi F; ///< @brief node indices of tetrahedra/triangle element
    Eigen::MatrixXi SF; ///< @brief node indices of surface triangles
    Eigen::MatrixXi CE; ///< @brief node indices of codimensional segments

    Eigen::MatrixXd V_prev; ///< @brief Vertices of the mesh at the previous time-step.
    Eigen::VectorXi SVI; ///< @brief node indices of surface nodes

    std::vector<int> componentNodeRange; ///< @brief V range of each loaded component [0(start0), end0(start1), end1(start2), ...]
    std::vector<int> componentSFRange; ///< @brief SF range of each loaded component [0(start0), end0(start1), end1(start2), ...]
    std::vector<int> componentCERange; ///< @brief CE range of each loaded component [0(start0), end0(start1), end1(start2), ...]
    std::vector<int> componentCoDim; ///< @brief codimension of each loaded component [coDim0, coDim1, ...]
    std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>> componentMaterial; ///< @brief E and nu of each loaded component [compI, startElemI, endElemI, density, E, nu]
    std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>> componentLVels; ///< @brief scripted linear velocity if any of each loaded component [compI, startElemI, endElemI, density, E, nu]
    std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>> componentAVels; ///< @brief scripted angular velocity if any of each loaded component [compI, startElemI, endElemI, density, E, nu]
    std::vector<std::pair<Eigen::Vector3i, std::array<Eigen::Vector3d, 2>>> componentInitVels; ///< @brief initial linear and angular velocity for any of each loaded component [compI, startElemI, endElemI, density, E, nu]
    std::vector<DirichletBC> DirichletBCs;
    std::vector<NeumannBC> NeumannBCs;
    std::vector<std::pair<int, std::string>> meshSeqFolderPath;

public: // owned features
    Eigen::SparseMatrix<double> massMatrix; // V.rows() wide
    double density, m_YM, m_PR;
    Eigen::VectorXd u, lambda;
    Eigen::VectorXd triArea; // triangle rest area
    double avgEdgeLen;
    std::set<int> DBCVertexIds; // for linear solve
    std::vector<DirichletBCType> vertexDBCType;
    Eigen::Matrix<double, 2, 3> bbox;
    std::vector<std::vector<int>> borderVerts_primitive;
    std::vector<Eigen::Matrix<double, dim, dim>> restTriInv;
    std::vector<bool> m_isBoundaryVert;

    // indices for fast access
    std::vector<std::set<int>> vNeighbor;
    std::vector<std::set<std::pair<int, int>>> vFLoc;
    std::vector<std::pair<int, int>> SFEdges;

public: // constructor
    // default constructor that doesn't do anything
    Mesh(void);

    // initialize from a triangle mesh, V will be constructed from UV_mesh in 2D,
    // V_mesh will be used to initialize restShape
    Mesh(const Eigen::MatrixXd& V_mesh,
        const Eigen::MatrixXi& F_mesh,
        const Eigen::MatrixXi& SF_mesh,
        const Eigen::MatrixXi& CE_mesh,
        const Eigen::MatrixXd& Vt_mesh,
        const std::vector<int>& p_componentNodeRange,
        const std::vector<int>& p_componentSFRange,
        const std::vector<int>& p_componentCERange,
        const std::vector<int>& p_componentCoDim,
        const std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>>& componentMaterial,
        const std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>>& componentLVels,
        const std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>>& componentAVels,
        const std::vector<std::pair<Eigen::Vector3i, std::array<Eigen::Vector3d, 2>>>& componentInitVels,
        const std::vector<DirichletBC>& DirichletBCs,
        const std::vector<NeumannBC>& NeumannBCs,
        const std::vector<std::pair<int, std::string>>& meshSeqFolderPath,
        double YM, double PR, double rho);

public: // API
    void computeMassMatrix(const igl::MassMatrixType type = igl::MASSMATRIX_TYPE_VORONOI);
    void computeFeatures(bool multiComp = false, bool resetDBCV = false);

    void resetDBCVertices(const std::map<int, DirichletBCType>& p_DBCVert = std::map<int, DirichletBCType>());
    void addDBCVertex(int vI, DirichletBCType type);
    void addDBCVertices(const std::vector<int>& p_DBCVert, DirichletBCType type);
    void addDBCVertices(const std::vector<std::pair<int, DirichletBCType>>& p_DBCVert);
    void removeDBCVertex(int vI);

    bool isDBCVertex(int vI) const { return vertexDBCType[vI] != DirichletBCType::NOT_DBC; }
    bool isProjectDBCVertex(int vI, bool projectDBC) const
    {
        assert(vI < vertexDBCType.size());
#ifdef EXPORT_FRICTION_DATA
        return vertexDBCType[vI] != DirichletBCType::NOT_DBC && projectDBC;
#else
        return vertexDBCType[vI] == DirichletBCType::ZERO
            || (vertexDBCType[vI] == DirichletBCType::NONZERO && projectDBC);
#endif
    }

    double avgNodeMass(int coDim = -1) const;
    double matSpaceBBoxSize2(int coDim = -1) const;

    int vICoDim(int vI) const;
    int sfICoDim(int sfI) const;

    void setLameParam(double YM, double PR);

    bool checkInversion(int triI, bool mute) const;
    bool checkInversion(bool mute = false, const std::vector<int>& triangles = std::vector<int>()) const;

    void save(const std::string& filePath, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        const Eigen::MatrixXd UV, const Eigen::MatrixXi& FUV = Eigen::MatrixXi()) const;
    void save(const std::string& filePath) const;

    void saveAsMesh(const std::string& filePath, bool scaleUV = false,
        const Eigen::MatrixXi& SF = Eigen::MatrixXi()) const;

    void saveSurfaceMesh(const std::string& filePath, bool use_V_prev = false, bool useInvShift = true) const;

    void saveBCNodes(const std::string& filePath) const;

public: // helper function
    bool isBoundaryVert(int vI) const;
    void computeBoundaryVert(const Eigen::MatrixXi& SF);
};

} // namespace IPC

#endif /* Mesh_hpp */
