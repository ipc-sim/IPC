//
//  Config.hpp
//  IPC
//
//  Created by Minchen Li on 7/12/18.
//

#ifndef Config_hpp
#define Config_hpp

#include "AnimScripter.hpp"
#include "CollisionObject.h"
#include "CollisionConstraints.hpp"
#include <ccd.hpp>

#include <iostream>
#include <map>
#include <memory>

namespace IPC {

enum EnergyType {
    ET_NH,
    ET_FCR
};

enum TimeIntegrationType {
    TIT_BE,
    TIT_NM
};

enum TimeStepperType {
    TST_NEWTON
};

enum ConstraintSolverType {
    CST_QP,
    CST_SQP,
    CST_IP
};

enum QPSolverType {
    QP_OSQP,
    QP_GUROBI
};

struct InputDBC {
    InputDBC(const Eigen::Vector3d& minBBox, const Eigen::Vector3d& maxBBox,
        const Eigen::Vector3d& linearVelocity, const Eigen::Vector3d& angularVelocity,
        const std::array<double, 2>& timeRange)
        : minBBox(minBBox), maxBBox(maxBBox), linearVelocity(linearVelocity), angularVelocity(angularVelocity), timeRange(timeRange)
    {
    }

    Eigen::Vector3d minBBox;
    Eigen::Vector3d maxBBox;
    Eigen::Vector3d linearVelocity;
    Eigen::Vector3d angularVelocity;
    std::array<double, 2> timeRange;
};

struct InputNBC {
    InputNBC(const Eigen::Vector3d& minBBox, const Eigen::Vector3d& maxBBox,
        const Eigen::Vector3d& force, const std::array<double, 2>& timeRange)
        : minBBox(minBBox), maxBBox(maxBBox), force(force), timeRange(timeRange)
    {
    }

    Eigen::Vector3d minBBox;
    Eigen::Vector3d maxBBox;
    Eigen::Vector3d force;
    std::array<double, 2> timeRange;
};

class Config {
public:
    std::string filePath;

    double size = -1;
    Eigen::Vector3d rotAxis;
    double rotDeg = 0.0;

    AnimScriptType animScriptType = AST_NULL;
    double handleRatio = 0.01;
    std::string meshSeqFolderPath;
    std::array<double, 2> DBCTimeRange = { 0.0, std::numeric_limits<double>::infinity() };
    std::array<double, 2> NBCTimeRange = { 0.0, std::numeric_limits<double>::infinity() };

    double rho = 1000.0;
    bool withGravity = true;
    EnergyType energyType = ET_NH;
    double YM = 1e5, PR = 0.4;

    TimeIntegrationType timeIntegrationType = TIT_BE;
    double beta = 0.25, gamma = 0.5;
    double duration = 5.0, dt = 0.025;
    double dampingRatio = 0.0, dampingStiff = 0.0;
    bool isConstrained = true, isSelfCollision = true;
    double selfFric = 0.0;

    std::vector<double> tol;
    int warmStart = 0;
    ConstraintSolverType constraintSolverType = CST_IP;
    double constraintOffset = 0.0;
    bool useActiveSetConvergence = false;

    std::vector<std::string> inputShapePaths;
    std::vector<Eigen::Vector3d> inputShapeTranslates;
    std::vector<Eigen::Matrix3d> inputShapeRotates;
    std::vector<Eigen::Vector3d> inputShapeScales;
    std::vector<Eigen::Vector3d> inputShapeMaterials;
    std::vector<Eigen::Vector3d> inputShapeLVels;
    std::vector<Eigen::Vector3d> inputShapeAVels;
    std::vector<std::array<Eigen::Vector3d, 2>> inputShapeInitVels;
    std::vector<std::pair<int, std::string>> inputShapeMeshSeqFolderPath;
    std::vector<std::pair<int, InputDBC>> inputShapeDBC;
    std::vector<std::pair<int, InputNBC>> inputShapeNBC;

    bool orthographic = false;
    double zoom = 1.0;
    bool cameraTracking = false;
    double playBackSpeed = 1.0;

    bool restart = false;
    std::string statusPath;
    bool disableCout = false;

    // collision objects
    std::vector<std::shared_ptr<CollisionObject<DIM>>> collisionObjects;
    std::vector<std::shared_ptr<CollisionObject<DIM>>> meshCollisionObjects;

    std::string appendStr;

    std::vector<double> tuning; // the parameter that is currently tuning
    bool useAbsParameters = false;
    double kappaMinMultiplier = 1e11;
    int fricIterAmt = 1;
    std::vector<double> scriptParams; // the parameters for setting script if any,
    // like initial velocity, position, etc

    ccd::CCDMethod ccdMethod = ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER;
    double ccdTolerance = 1e-6;

    /// @brief Constraint type for SQP method of handling collisions.
    CollisionConstraintType constraintType = CollisionConstraintType::VOLUME;

    /// @brief QP solver for use with SQP collision handling.
    QPSolverType qpSolverType = QPSolverType::QP_OSQP;

public:
    static const std::vector<std::string> energyTypeStrs;
    static const std::vector<std::string> timeIntegrationTypeStrs;
    static const std::vector<std::string> constraintSolverTypeStrs;
    static const std::vector<std::string> constraintTypeStrs;
    static const std::vector<std::string> QPSolverTypeStrs;

public:
    Config(void);
    ~Config(void);
    int loadFromFile(const std::string& filePath);
    void backUpConfig(const std::string& filePath);

public:
    void appendInfoStr(std::string& inputStr) const;

public:
    static EnergyType getEnergyTypeByStr(const std::string& str);
    static std::string getStrByEnergyType(EnergyType energyType);
    static TimeIntegrationType getTimeIntegrationTypeByStr(const std::string& str);
    static std::string getStrByTimeIntegrationType(TimeIntegrationType timeIntegrationType);
    static ConstraintSolverType getConstraintSolverTypeByStr(const std::string& str);
    static std::string getStrByConstraintSolverType(ConstraintSolverType constraintSolverType);
    static ccd::CCDMethod getCCDMethodTypeByStr(const std::string& str);
    static CollisionConstraintType getConstraintTypeByStr(const std::string& str);
    static QPSolverType getQPSolverTypeByStr(const std::string& str);
};

} // namespace IPC

#endif /* Config_hpp */
