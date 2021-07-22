//
//  AnimScripter.hpp
//  IPC
//
//  Created by Minchen Li on 6/20/18.
//

#ifndef AnimScripter_hpp
#define AnimScripter_hpp

#include "Energy.hpp"
#include "Mesh.hpp"
#include "LinSysSolver.hpp"
#include "CollisionObject.h"

#include <cstdio>
#include <map>
#include <unordered_map>

namespace IPC {

enum AnimScriptType {
    AST_NULL,
    AST_SCALEF,
    AST_HANG,
    AST_HANG2,
    AST_HANGTOPLEFT,
    AST_HANGLEFT,
    AST_SWING,
    AST_STAMP,
    AST_STAMPTOPLEFT,
    AST_STAMPBOTH,
    AST_UNDSTAMP,
    AST_STAMPINV,
    AST_STAND,
    AST_STANDINV,
    AST_TOPBOTTOMFIX,
    AST_FIXLOWERHALF,
    AST_CORNER,
    AST_PUSH,
    AST_TEAR,
    AST_UPNDOWN,
    AST_STRETCH,
    AST_STRETCHNPAUSE,
    AST_SQUASH,
    AST_STRETCHNSQUASH,
    AST_BEND,
    AST_TWIST,
    AST_TWISTNSTRETCH,
    AST_TWISTNSNS,
    AST_TWISTNSNS_OLD,
    AST_RUBBERBANDPULL,
    AST_FOURLEGPULL,
    AST_HEADTAILPULL,
    AST_ONEPOINT,
    AST_RANDOM,
    AST_FALL,
    AST_FALL_NOSHIFT,
    AST_DRAGDOWN,
    AST_DRAGRIGHT,
    AST_TOGGLETOP,
    AST_LEFTHITRIGHT,
    AST_DROP,
    AST_XYROTATE,
    AST_CURTAIN,
    AST_FIXRIGHTMOST1,
    AST_PUSHRIGHTMOST1,

    AST_NMFIXBOTTOMDRAGLEFT,
    AST_NMFIXBOTTOMDRAGFORWARD,

    AST_ACOSQUASH,
    AST_ACOSQUASH6,
    AST_ACOSQUASHSHEAR,
    AST_MCOSQUASH,
    AST_MCOROTSQUASH,
    AST_MCOROTCYLINDERS,

    AST_DCOROTCYLINDERS,
    AST_DCOSQUASH,
    AST_DCOSQUASH6,
    AST_DCOFIX,
    AST_DCOSEGBEDSQUASH,
    AST_DCOBALLHITWALL,
    AST_DCOSQUEEZEOUT,
    AST_DCOVERSCHOORROLLER,
    AST_DCOHAMMERWALNUT,
    AST_DCOCUT,

    AST_UTOPIA_COMPARISON,

    AST_MESHSEQ_FROMFILE
};

template <int dim>
class AnimScripter {
protected:
    AnimScriptType animScriptType;

    std::vector<std::vector<int>> handleVerts;

    std::map<int, Eigen::Matrix<double, dim, 1>> velocity_handleVerts;
    std::pair<int, Eigen::Matrix<double, dim, 2>> velocityTurningPoints;

    std::map<int, double> angVel_handleVerts;
    std::map<int, Eigen::Matrix<double, dim, 1>> rotCenter_handleVerts;

    Eigen::VectorXd searchDir;
    double stepSize;
    std::unordered_map<int, std::pair<Eigen::Matrix<double, 1, dim>, Eigen::Matrix<double, 1, dim>>> targetPos;
    double dist2Tol;

    std::vector<Eigen::Matrix<double, dim, 1>> ACOVelocity;
    std::vector<Eigen::Matrix<double, dim, 1>> MCOVelocity;
    std::vector<double> ACOStepSizeLeft;
    std::vector<double> MCOStepSizeLeft;
    std::vector<Eigen::Matrix<double, dim, 1>> MCORotCenter;
    std::vector<Eigen::Matrix<double, dim, 1>> MCOAngVel;

    std::string meshSeqFolderPath;
    int meshI;

    double stepStartTime = 0, stepEndTime = 0;
    std::array<double, 2> DBCTimeRange = { 0.0, std::numeric_limits<double>::infinity() };
    std::array<double, 2> NBCTimeRange = { 0.0, std::numeric_limits<double>::infinity() };
    std::set<int> DBCVertexIdsBK;

protected:
    static const std::vector<std::string> animScriptTypeStrs;

public:
    AnimScripter(AnimScriptType p_animScriptType = AST_NULL);

public:
    void initAnimScript(Mesh<dim>& mesh,
        const std::vector<std::shared_ptr<CollisionObject<dim>>>& ACO,
        const std::vector<std::shared_ptr<CollisionObject<dim>>>& MCO,
        const std::array<double, 2>& DBCTimeRange,
        const std::array<double, 2>& NBCTimeRange);
    void initVelocity(const Mesh<dim>& mesh, const std::vector<double>& params,
        Eigen::VectorXd& velocity) const;
    int stepAnimScript(Mesh<dim>& mesh,
        const SpatialHash<dim>& sh,
        std::vector<std::shared_ptr<CollisionObject<dim>>>& ACO,
        std::vector<std::shared_ptr<CollisionObject<dim>>>& MCO,
        const ccd::CCDMethod method,
        double dt, double dHat, const std::vector<Energy<dim>*>& energyTerms,
        bool isSelfCollision, bool forceIntersectionLineSearch);
    double getCompletedStepSize(void) const;
    double computeCompletedStepSize(const Mesh<dim>& mesh);
    void augmentMDBCEnergy(const Mesh<dim>& data, double& E, double rho) const;
    void augmentMDBCGradient(const Mesh<dim>& data, Eigen::VectorXd& g, double rho) const;
    void augmentMDBCHessian(const Mesh<dim>& data, LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver, double rho) const;
    void updateLambda(const Mesh<dim>& data, double rho);
    void resetLambda(void);

    double getCOCompletedStepSize(void) const;

    bool isNBCActive(void) const;
    bool isNBCActive(const Mesh<dim>& mesh, int NBCi) const;

public:
    void setAnimScriptType(AnimScriptType p_animScriptType);
    void setMeshSeqFolderPath(const std::string& p_meshSeqFolderPath);

public:
    static AnimScriptType getAnimScriptTypeByStr(const std::string& str);
    static std::string getStrByAnimScriptType(AnimScriptType animScriptType);

protected:
    void setDBCVertices(Mesh<dim>& mesh) const;
};

} // namespace IPC

#endif /* AnimScripter_hpp */
