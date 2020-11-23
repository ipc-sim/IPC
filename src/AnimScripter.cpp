//
//  AnimScripter.cpp
//  IPC
//
//  Created by Minchen Li on 6/20/18.
//

#include "AnimScripter.hpp"

#include "MeshCO.hpp"
#include "SelfCollisionHandler.hpp"

#include "BarrierFunctions.hpp"
#include "get_feasible_steps.hpp"
#include "IglUtils.hpp"

#ifdef USE_TBB
#include <tbb/tbb.h>
#endif

#include <spdlog/spdlog.h>

#include <igl/readOBJ.h>

extern Timer timer_temp3;
extern std::string outputFolderPath;
extern int numOfCCDFail;
extern std::ofstream logFile;

namespace IPC {

template <int dim>
const std::vector<std::string> AnimScripter<dim>::animScriptTypeStrs = {
    "null", "scaleF", "hang", "hang2", "hangTopLeft", "hangLeft", "swing",
    "stamp", "stampTopLeft", "stampBoth", "undstamp", "stampInv",
    "stand", "standInv", "topbottomfix", "fixLowerHalf",
    "corner", "push", "tear", "upndown", "stretch", "stretchAndPause", "squash", "stretchnsquash",
    "bend", "twist", "twistnstretch", "twistnsns", "twistnsns_old",
    "rubberBandPull", "fourLegPull", "headTailPull",
    "onepoint", "random", "fall", "fallNoShift", "dragdown", "dragright", "toggleTop",
    "leftHitRight", "drop", "XYRotate", "curtain",
    "fixRightMost1", "pushRightMost1",
    "NMFixBottomDragLeft", "NMFixBottomDragForward",
    "ACOSquash", "ACOSquash6", "ACOSquashShear",
    "MCOSquash", "MCORotSquash", "MCORotCylinders",
    "DCORotCylinders", "DCOSquash", "DCOSquash6", "DCOFix", "DCOSegBedSquash", "DCOBallHitWall",
    "DCOSqueezeOut", "DCOVerschoorRoller", "DCOHammerWalnut", "DCOCut",
    "utopiaComparison",
    "meshSeqFromFile"
};

template <int dim>
AnimScripter<dim>::AnimScripter(AnimScriptType p_animScriptType)
    : animScriptType(p_animScriptType), meshI(0)
{
}

template <int dim>
void AnimScripter<dim>::initAnimScript(Mesh<dim>& mesh,
    const std::vector<CollisionObject<dim>*>& ACO,
    const std::vector<CollisionObject<dim>*>& MCO,
    double DBCTimeRange[2],
    double NBCTimeRange[2])
{
    this->DBCTimeRange[0] = DBCTimeRange[0];
    this->DBCTimeRange[1] = DBCTimeRange[1];
    this->NBCTimeRange[0] = NBCTimeRange[0];
    this->NBCTimeRange[1] = NBCTimeRange[1];

    switch (animScriptType) {
    case AST_NULL: {
        mesh.resetFixedVert();
        std::vector<int> fixedVerts;

        // fix codim objects
        for (int compI = 0; compI < mesh.componentCoDim.size(); ++compI) {
            if (mesh.componentCoDim[compI] < 3) {
                for (int i = mesh.componentNodeRange[compI]; i < mesh.componentNodeRange[compI + 1]; ++i) {
                    fixedVerts.emplace_back(i);
                }
            }
        }

        // fix objects with scripted linear and angular velocity
        for (const auto& scriptedLVelI : mesh.componentLVels) {
            for (int i = scriptedLVelI.first[1]; i < scriptedLVelI.first[2]; ++i) {
                fixedVerts.emplace_back(i);
            }
        }
        for (const auto& scriptedAVelI : mesh.componentAVels) {
            for (int i = scriptedAVelI.first[1]; i < scriptedAVelI.first[2]; ++i) {
                fixedVerts.emplace_back(i);
            }
        }

        // fix Dirchlet nodes
        for (const auto& DBCInfoI : mesh.DBCInfo) {
            for (const auto& vI : DBCInfoI.first) {
                fixedVerts.emplace_back(vI);
            }
        }

        mesh.addFixedVert(fixedVerts);
        break;
    }

    case AST_SCALEF: {
        mesh.resetFixedVert();
        //                handleVerts.resize(0);
        //                int bI = 0;
        //                for(const auto borderI : mesh.borderVerts_primitive) {
        //                    mesh.addFixedVert(borderI);
        //                    handleVerts.emplace_back(borderI);
        //                    bI++;
        //                }

        Eigen::Matrix3d M;
        //                M << 1.5, 0.5, 0.0,
        //                0.0, 0.5, -0.5,
        //                0.0, 0.0, 1.0;
        M << 1.5, 0.0, 0.0,
            0.0, 1.5, 0.0,
            0.0, 0.0, 1.5;
        for (int i = 0; i < mesh.V.rows(); ++i) {
            mesh.V.row(i) = (M * mesh.V.row(i).transpose()).transpose();
        }

        break;
    }

    case AST_HANG:
        mesh.resetFixedVert();
        for (const auto borderI : mesh.borderVerts_primitive) {
            mesh.addFixedVert(borderI.back());
        }
        break;

    case AST_HANG2: {
        mesh.resetFixedVert();

        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;

        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.01) {
                mesh.addFixedVert(vI);
                velocity_handleVerts[vI].setZero();
            }
        }
        break;
    }

    case AST_HANGTOPLEFT: {
        mesh.resetFixedVert();

        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;

        for (const auto& vI : mesh.borderVerts_primitive[0]) {
            if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.01 && (mesh.V(vI, 2) > topRight[2] - range[2] * 0.01 || mesh.V(vI, 2) < bottomLeft[2] + range[2] * 0.01)) {
                mesh.addFixedVert(vI);
            }
        }
        break;
    }

    case AST_HANGLEFT: {
        mesh.resetFixedVert();

        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;

        for (const auto& vI : mesh.borderVerts_primitive[0]) {
            if ((mesh.V(vI, 2) <= topRight[2]) && (mesh.V(vI, 2) >= bottomLeft[2])) {
                mesh.addFixedVert(vI);
            }
        }
        break;
    }

    case AST_SWING: {
        mesh.V.col(1).array() += 1.3 * (mesh.V.col(1).maxCoeff() - mesh.V.col(1).minCoeff());

        mesh.resetFixedVert();

        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;

        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 0) < bottomLeft[0] + range[0] * 0.05) {
                mesh.addFixedVert(vI);
                velocity_handleVerts[vI].setZero();
            }
        }
        break;
    }

    case AST_STAMP:
        mesh.resetFixedVert();
        mesh.addFixedVert(mesh.borderVerts_primitive[0]);
        break;

    case AST_STAMPTOPLEFT: {
        mesh.resetFixedVert();

        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;

        for (const auto& vI : mesh.borderVerts_primitive[0]) {
            if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.01) {
                mesh.addFixedVert(vI);
            }
        }
        break;
    }

    case AST_STAMPBOTH:
        mesh.resetFixedVert();
        mesh.addFixedVert(mesh.borderVerts_primitive[0]);
        mesh.addFixedVert(mesh.borderVerts_primitive[1]);
        break;

    case AST_UNDSTAMP:
        mesh.resetFixedVert();
        mesh.addFixedVert(mesh.borderVerts_primitive[0]);
        handleVerts.emplace_back(mesh.borderVerts_primitive[0]);

        for (const auto bVI : mesh.borderVerts_primitive[0]) {
            velocity_handleVerts[bVI].setZero();
            velocity_handleVerts[bVI][1] = 1.8;
        }

        velocityTurningPoints.first = mesh.borderVerts_primitive[0].front();
        velocityTurningPoints.second(1, 0) = mesh.V(velocityTurningPoints.first, 1) - 0.6;
        velocityTurningPoints.second(1, 1) = mesh.V(velocityTurningPoints.first, 1) + 0.6;
        break;

    case AST_STAMPINV: {
        mesh.resetFixedVert();

        // fix bottom
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 0) < bottomLeft[0] + range[0] * 0.01) {
                mesh.addFixedVert(vI);
            }
        }

        double xOffset = 1.1 * mesh.V(*mesh.fixedVert.begin(), 0);
        mesh.V.col(0) *= -0.1;
        mesh.V.col(0).array() += xOffset;
        break;
    }

    case AST_STAND: {
        mesh.resetFixedVert();

        // fix bottom
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.01) {
                mesh.addFixedVert(vI);
            }
        }
        break;
    }

    case AST_STANDINV: {
        mesh.resetFixedVert();

        // fix bottom
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.01) {
                mesh.addFixedVert(vI);
            }
        }

        double yOffset = 1.1 * mesh.V(*mesh.fixedVert.begin(), 1);
        mesh.V.col(1) *= -0.1;
        mesh.V.col(1).array() += yOffset;
        break;
    }

    case AST_TOPBOTTOMFIX: {
        mesh.resetFixedVert();

        // grab top, waist, and bottom:
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.02) {
                mesh.addFixedVert(vI);
            }
            else if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.02) {
                mesh.addFixedVert(vI);
            }
        }

        break;
    }

    case AST_FIXLOWERHALF: {
        mesh.resetFixedVert();

        // grab lower half:
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.5) {
                mesh.addFixedVert(vI);
            }
        }

        break;
    }

    case AST_CORNER: {
        mesh.resetFixedVert();

        // grab top, waist, and bottom:
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 0) < bottomLeft[0] + range[0] * 0.01) {
                mesh.addFixedVert(vI);
            }
            else if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.01) {
                mesh.addFixedVert(vI);
            }
            else if (mesh.V(vI, 2) < bottomLeft[2] + range[2] * 0.01) {
                mesh.addFixedVert(vI);
            }
        }

        break;
    }

    case AST_PUSH: {
        mesh.resetFixedVert();
        handleVerts.resize(0);
        handleVerts.resize(1);

        // fix bottom
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        bool turningPointSet = false;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.01) {
                mesh.addFixedVert(vI);
                velocity_handleVerts[vI].setZero();
            }
            else if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.01) {
                mesh.addFixedVert(vI);
                handleVerts[0].emplace_back(vI);

                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][1] = -1.0;
                if (!turningPointSet) {
                    turningPointSet = true;
                    velocityTurningPoints.first = vI;
                    velocityTurningPoints.second(0, 1) = mesh.V(vI, 1) - 0.5;
                }
            }
        }
        break;
    }

    case AST_TEAR: {
        mesh.resetFixedVert();
        handleVerts.resize(0);
        handleVerts.resize(1);

        // fix bottom
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        bool turningPointSet = false;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.01) {
                mesh.addFixedVert(vI);
                velocity_handleVerts[vI].setZero();
            }
            else if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.01) {
                mesh.addFixedVert(vI);
                handleVerts[0].emplace_back(vI);

                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = -5;
                if (!turningPointSet) {
                    turningPointSet = true;
                    velocityTurningPoints.first = vI;
                    velocityTurningPoints.second(0, 0) = mesh.V(vI, 0) - 4;
                }
            }
        }
        break;
    }

    case AST_UPNDOWN: {
        mesh.resetFixedVert();

        int bI = 0;
        for (const auto borderI : mesh.borderVerts_primitive) {
            mesh.addFixedVert(borderI);
            handleVerts.emplace_back(borderI);

            for (const auto bVI : borderI) {
                velocity_handleVerts[bVI].setZero();
                velocity_handleVerts[bVI][1] = std::pow(-1.0, bI) * 1.8;
            }
            ++bI;
        }

        velocityTurningPoints.first = mesh.borderVerts_primitive[0].front();
        velocityTurningPoints.second(1, 0) = mesh.V(velocityTurningPoints.first, 1) - 0.6;
        velocityTurningPoints.second(1, 1) = mesh.V(velocityTurningPoints.first, 1) + 0.6;

        break;
    }

    case AST_STRETCH: {
        mesh.resetFixedVert();
        handleVerts.resize(0);
        int bI = 0;
        for (const auto borderI : mesh.borderVerts_primitive) {
            mesh.addFixedVert(borderI);
            handleVerts.emplace_back(borderI);
            for (const auto bVI : borderI) {
                velocity_handleVerts[bVI].setZero();
                velocity_handleVerts[bVI][0] = std::pow(-1.0, bI) * -0.1;
            }
            bI++;
        }
        break;
    }

    case AST_STRETCHNPAUSE: {
        mesh.resetFixedVert();
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        bool turningPointAdded = false;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 0) < bottomLeft[0] + range[0] * 0.01) {
                mesh.addFixedVert(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = -1.0;
                if (!turningPointAdded) {
                    velocityTurningPoints.first = vI;
                    velocityTurningPoints.second.setZero();
                    velocityTurningPoints.second(0, 0) = -0.28;
                }
            }
            else if (mesh.V(vI, 0) > topRight[0] - range[0] * 0.01) {
                mesh.addFixedVert(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = 1.0;
            }
        }

        break;
    }

    case AST_SQUASH: {
        mesh.resetFixedVert();
        handleVerts.resize(0);
        int bI = 0;
        for (const auto borderI : mesh.borderVerts_primitive) {
            mesh.addFixedVert(borderI);
            handleVerts.emplace_back(borderI);
            for (const auto bVI : borderI) {
                velocity_handleVerts[bVI].setZero();
                velocity_handleVerts[bVI][0] = std::pow(-1.0, bI) * 0.03;
            }
            bI++;
        }
        break;
    }

    case AST_STRETCHNSQUASH: {
        mesh.resetFixedVert();
        handleVerts.resize(0);
        int bI = 0;
        for (const auto borderI : mesh.borderVerts_primitive) {
            mesh.addFixedVert(borderI);
            handleVerts.emplace_back(borderI);
            for (const auto bVI : borderI) {
                velocity_handleVerts[bVI].setZero();
                velocity_handleVerts[bVI][0] = std::pow(-1.0, bI) * -0.9;
            }
            bI++;
        }

        velocityTurningPoints.first = mesh.borderVerts_primitive[0].front();
        velocityTurningPoints.second(0, 0) = mesh.V(velocityTurningPoints.first, 0) - 0.8;
        velocityTurningPoints.second(0, 1) = mesh.V(velocityTurningPoints.first, 0) + 0.4;

        break;
    }

    case AST_BEND: {
        mesh.resetFixedVert();
        handleVerts.resize(0);
        int bI = 0;
        for (const auto borderI : mesh.borderVerts_primitive) {
            mesh.addFixedVert(borderI);
            handleVerts.emplace_back(borderI);
            for (int bVI = 0; bVI + 1 < borderI.size(); bVI++) {
                angVel_handleVerts[borderI[bVI]] = std::pow(-1.0, bI) * -0.05 * M_PI;
                rotCenter_handleVerts[borderI[bVI]] = mesh.V.row(borderI.back()).transpose();
            }
            bI++;
        }
        break;
    }

    case AST_TWIST: {
        mesh.resetFixedVert();

        const Eigen::RowVector3d rotCenter = mesh.bbox.colwise().mean();

        handleVerts.resize(0);
        int bI = 0;
        for (const auto borderI : mesh.borderVerts_primitive) {
            mesh.addFixedVert(borderI);
            handleVerts.emplace_back(borderI);
            for (int bVI = 0; bVI < borderI.size(); bVI++) {
                angVel_handleVerts[borderI[bVI]] = std::pow(-1.0, bI) * -0.4 * M_PI;
                rotCenter_handleVerts[borderI[bVI]] = rotCenter.transpose().topRows(dim);
            }
            bI++;
        }
        break;
    }

    case AST_TWISTNSTRETCH: {
        mesh.resetFixedVert();

        const Eigen::RowVector3d rotCenter = mesh.bbox.colwise().mean();

        handleVerts.resize(0);
        int bI = 0;
        for (const auto borderI : mesh.borderVerts_primitive) {
            mesh.addFixedVert(borderI);
            handleVerts.emplace_back(borderI);
            for (int bVI = 0; bVI < borderI.size(); bVI++) {
                angVel_handleVerts[borderI[bVI]] = std::pow(-1.0, bI) * -0.1 * M_PI;
                rotCenter_handleVerts[borderI[bVI]] = rotCenter.transpose().topRows(dim);

                velocity_handleVerts[borderI[bVI]].setZero();
                velocity_handleVerts[borderI[bVI]][0] = std::pow(-1.0, bI) * -0.1;
            }
            bI++;
        }
        break;
    }

    case AST_TWISTNSNS_OLD:
    case AST_TWISTNSNS: {
        mesh.resetFixedVert();

        const Eigen::RowVector3d rotCenter = mesh.bbox.colwise().mean();

        handleVerts.resize(0);
        int bI = 0;
        for (const auto borderI : mesh.borderVerts_primitive) {
            mesh.addFixedVert(borderI);
            handleVerts.emplace_back(borderI);
            for (int bVI = 0; bVI < borderI.size(); bVI++) {
                angVel_handleVerts[borderI[bVI]] = std::pow(-1.0, bI) * -0.4 * M_PI;
                rotCenter_handleVerts[borderI[bVI]] = rotCenter.transpose().topRows(dim);

                velocity_handleVerts[borderI[bVI]].setZero();
                if (animScriptType == AST_TWISTNSNS) {
                    velocity_handleVerts[borderI[bVI]][0] = std::pow(-1.0, bI) * -1.2;
                }
                else if (animScriptType == AST_TWISTNSNS_OLD) {
                    velocity_handleVerts[borderI[bVI]][0] = std::pow(-1.0, bI) * -0.9;
                }
            }
            bI++;
        }

        velocityTurningPoints.first = mesh.borderVerts_primitive[0].front();
        if (animScriptType == AST_TWISTNSNS) {
            velocityTurningPoints.second(0, 0) = mesh.V(velocityTurningPoints.first, 0) - 1.2;
        }
        else {
            velocityTurningPoints.second(0, 0) = mesh.V(velocityTurningPoints.first, 0) - 0.8;
        }
        velocityTurningPoints.second(0, 1) = mesh.V(velocityTurningPoints.first, 0) + 0.4;
        break;
    }

    case AST_RUBBERBANDPULL: {
        mesh.resetFixedVert();
        handleVerts.resize(0);
        handleVerts.resize(2);

        // grab top, waist, and bottom:
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        bool turningPointSet = false;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.02) {
                mesh.addFixedVert(vI);
                handleVerts[1].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][1] = -0.2;
            }
            else if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.02) {
                mesh.addFixedVert(vI);
                handleVerts[1].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][1] = 0.2;
            }
            else if ((mesh.V(vI, 1) < topRight[1] - range[1] * 0.48) && (mesh.V(vI, 1) > bottomLeft[1] + range[1] * 0.48)) {
                mesh.addFixedVert(vI);
                handleVerts[0].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = -2.5; // previously -2.0
                if (!turningPointSet) {
                    turningPointSet = true;
                    velocityTurningPoints.first = vI;
                    velocityTurningPoints.second(0, 0) = mesh.V(vI, 0) - 5.0;
                }
            }
        }

        break;
    }

    case AST_FOURLEGPULL: {
        mesh.resetFixedVert();
        handleVerts.resize(0);
        handleVerts.resize(2);

        // grab top, waist, and bottom:
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        bool turningPointSet = false;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if ((mesh.V(vI, 1) > topRight[1] - range[1] * 0.129) && ((mesh.V(vI, 0) < bottomLeft[0] + range[0] * 0.16))) {
                // left top
                mesh.addFixedVert(vI);
                handleVerts[1].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
            }
            else if ((mesh.V(vI, 1) > topRight[1] - range[1] * 0.16) && (mesh.V(vI, 0) > topRight[0] - range[0] * 0.16)) {
                // right top
                mesh.addFixedVert(vI);
                handleVerts[0].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = 2.5;
            }
            else if ((mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.02) && (mesh.V(vI, 0) > topRight[0] - range[0] * 0.25)) {
                // bottom right
                mesh.addFixedVert(vI);
                handleVerts[0].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = 2.5;
                velocity_handleVerts[vI][1] = -3.5;
            }
            else if ((mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.02) && (mesh.V(vI, 0) < bottomLeft[0] + range[0] * 0.25)) {
                // bottom left
                mesh.addFixedVert(vI);
                handleVerts[0].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][1] = -3.5;

                if (!turningPointSet) {
                    turningPointSet = true;
                    velocityTurningPoints.first = vI;
                    velocityTurningPoints.second(0, 1) = mesh.V(vI, 1) - 5.0;
                }
            }
        }
        break;
    }

    case AST_HEADTAILPULL: {
        mesh.resetFixedVert();
        handleVerts.resize(0);
        handleVerts.resize(2);

        // grab top, waist, and bottom:
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        bool turningPointSet = false;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 2) < bottomLeft[2] + range[2] * 0.02) {
                // left
                mesh.addFixedVert(vI);
                handleVerts[0].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = 3.5;

                if (!turningPointSet) {
                    turningPointSet = true;
                    velocityTurningPoints.first = vI;
                    velocityTurningPoints.second(0, 0) = mesh.V(vI, 0) + 4.5;
                }
            }
            else if (mesh.V(vI, 2) > topRight[2] - range[2] * 0.02) {
                // right
                mesh.addFixedVert(vI);
                handleVerts[0].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = 3.5;
            }
            else if ((mesh.V(vI, 2) > bottomLeft[2] + range[2] * 0.46) && (mesh.V(vI, 2) < bottomLeft[2] + range[2] * 0.54)) {
                // bottom right
                mesh.addFixedVert(vI);
                handleVerts[1].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
            }
        }
        break;
    }

    case AST_ONEPOINT: {
        const Eigen::RowVector3d center = mesh.bbox.colwise().mean();
        mesh.V.rowwise() = center.leftCols(dim);
        mesh.V.col(1).array() += (mesh.bbox(1, 1) - mesh.bbox(0, 1)) / 2.0;
        break;
    }

    case AST_RANDOM: {
        mesh.V.setRandom();
        mesh.V /= 2.0;
        Eigen::RowVector3d offset = mesh.bbox.colwise().mean();
        offset[1] += (mesh.bbox(1, 1) - mesh.bbox(0, 1)) / 2.0;
        offset.leftCols(dim) -= mesh.V.row(0);
        mesh.V.rowwise() += offset.leftCols(dim);
        break;
    }

    case AST_FALL: {
        mesh.V.col(1).array() += 0.5 * (mesh.V.colwise().maxCoeff() - mesh.V.colwise().minCoeff()).norm();
        mesh.resetFixedVert();
        break;
    }

    case AST_FALL_NOSHIFT: {
        mesh.resetFixedVert();
        break;
    }

    case AST_DRAGDOWN: {
        mesh.V.col(1).array() += 0.5 * (mesh.V.colwise().maxCoeff() - mesh.V.colwise().minCoeff()).norm();
        mesh.resetFixedVert();

        // grab bottom:
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if ((mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.1) && (mesh.V(vI, 0) < bottomLeft[0] + range[0] * 0.52) && (mesh.V(vI, 0) > bottomLeft[0] + range[0] * 0.42)) {
                // bottom middle
                mesh.addFixedVert(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][1] = -1.5;
            }
        }
        break;
    }

    case AST_DRAGRIGHT: {
        mesh.V.col(1).array() += 0.5 * (mesh.V.colwise().maxCoeff() - mesh.V.colwise().minCoeff()).norm();
        mesh.resetFixedVert();

        // grab right:
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 0) > topRight[0] - range[0] * 0.04) { // 0.001 for dolphin rightmost 2 nodes
                mesh.addFixedVert(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = 0.5;
            }
        }
        break;
    }

    case AST_TOGGLETOP: {
        mesh.resetFixedVert();
        handleVerts.resize(0);
        handleVerts.resize(1);

        // grab top
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        bool turningPointSet = false;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.02) {
                mesh.addFixedVert(vI);
                handleVerts[0].emplace_back(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = -0.5;
                if (!turningPointSet) {
                    turningPointSet = true;
                    velocityTurningPoints.first = vI;
                    velocityTurningPoints.second(0, 0) = mesh.V(vI, 0) - 0.1;
                }
            }
        }
        break;
    }

    case AST_LEFTHITRIGHT:
    case AST_DROP:
    case AST_XYROTATE: {
        mesh.resetFixedVert();
        break;
    }

    case AST_CURTAIN: {
        mesh.resetFixedVert();

        // grab top nodes
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            for (int pinI = 0; pinI < 8; ++pinI) {
                if (mesh.V(vI, 0) > bottomLeft[0] + range[0] / 7.0 * pinI - range[0] * 0.0025 && mesh.V(vI, 0) < bottomLeft[0] + range[0] / 7.0 * pinI + range[0] * 0.0025 && mesh.V(vI, 1) > topRight[1] - range[1] * 0.005) {
                    mesh.addFixedVert(vI);
                    velocity_handleVerts[vI].setZero();
                    velocity_handleVerts[vI][0] = 0.04 * (7.0 - pinI) / 7.0;
                    break;
                }
            }
        }
        break;
    }

    case AST_FIXRIGHTMOST1: {
        mesh.resetFixedVert();

        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 0) > topRight[0] - 1.0e-3 * range[0]) {
                mesh.addFixedVert(vI);
                break;
            }
        }
        break;
    }

    case AST_PUSHRIGHTMOST1: {
        mesh.resetFixedVert();

        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 0) > topRight[0] - 1.0e-3 * range[0]) {
                mesh.addFixedVert(vI);
                velocity_handleVerts[vI].setZero();
                velocity_handleVerts[vI][0] = -0.15;
                break;
            }
        }
        break;
    }

    case AST_NMFIXBOTTOMDRAGLEFT: {
        mesh.resetFixedVert();

        // grab bottom and top nodes
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.05) {
                mesh.addFixedVert(vI);
            }
            else if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.05) {
                mesh.NeumannBC[vI].setZero();
                mesh.NeumannBC[vI][0] = -600.0;
            }
        }
        break;
    }

    case AST_NMFIXBOTTOMDRAGFORWARD: {
        mesh.resetFixedVert();

        // grab bottom and top nodes
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.05) {
                mesh.addFixedVert(vI);
            }
            else if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.05) {
                mesh.NeumannBC[vI].setZero();
                mesh.NeumannBC[vI][2] = 600.0;
            }
        }
        break;
    }

    case AST_ACOSQUASH: {
        mesh.resetFixedVert();
        ACOVelocity.resize(2);
        ACOVelocity[0].setZero();
        ACOVelocity[0][0] = 1.0;
        ACOVelocity[1].setZero();
        ACOVelocity[1][0] = -1.0;
        break;
    }

    case AST_ACOSQUASH6: {
        mesh.resetFixedVert();
        ACOVelocity.resize(6);

        ACOVelocity[0].setZero();
        ACOVelocity[0][0] = 1.0;
        ACOVelocity[1].setZero();
        ACOVelocity[1][0] = -1.0;

        ACOVelocity[2].setZero();
        ACOVelocity[2][1] = 1.0;
        ACOVelocity[3].setZero();
        ACOVelocity[3][1] = -1.0;

        ACOVelocity[4].setZero();
        ACOVelocity[4][2] = 1.0;
        ACOVelocity[5].setZero();
        ACOVelocity[5][2] = -1.0;
        break;
    }

    case AST_ACOSQUASHSHEAR: {
        mesh.resetFixedVert();
        ACOVelocity.resize(2);
        ACOVelocity[0].setZero();
        ACOVelocity[0][0] = 1.0;
        ACOVelocity[1].setZero();
        ACOVelocity[1][0] = -1.0;
        break;
    }

    case AST_MCOSQUASH: {
        mesh.resetFixedVert();
        MCOVelocity.resize(2);
        MCOVelocity[0].setZero();
        MCOVelocity[0][0] = 1.0;
        MCOVelocity[1].setZero();
        MCOVelocity[1][0] = -1.0;
        break;
    }

    case AST_MCOROTSQUASH: {
        mesh.resetFixedVert();

        MCOVelocity.resize(2);
        MCOVelocity[0].setZero();
        MCOVelocity[1].setZero();

        assert(MCO.size() >= 2);

        MCORotCenter.resize(2);
        MCORotCenter[0] = MCO[0]->origin;
        MCORotCenter[0][1] += dynamic_cast<MeshCO<dim>*>(MCO[0])->scale / 2.0;
        MCORotCenter[1] = MCO[1]->origin;
        MCORotCenter[1][1] -= dynamic_cast<MeshCO<dim>*>(MCO[1])->scale / 2.0;

        MCOAngVel.resize(2);
        MCOAngVel[0].setZero();
        MCOAngVel[0][2] = M_PI;
        MCOAngVel[1].setZero();
        MCOAngVel[1][2] = M_PI;
        break;
    }

    case AST_MCOROTCYLINDERS: {
        mesh.resetFixedVert();

        MCOVelocity.resize(4);
        MCOVelocity[0].setZero();
        MCOVelocity[1].setZero();
        MCOVelocity[2].setZero();
        MCOVelocity[3].setZero();

        assert(MCO.size() >= 4);

        MCORotCenter.resize(4);
        MCORotCenter[0] = MCO[0]->origin;
        MCORotCenter[1] = MCO[1]->origin;
        MCORotCenter[2] = MCO[2]->origin;
        MCORotCenter[3] = MCO[3]->origin;

        MCOAngVel.resize(4);
        MCOAngVel[0].setZero();
        MCOAngVel[0][0] = M_PI / 2.0;
        MCOAngVel[1].setZero();
        MCOAngVel[1][0] = -M_PI / 2.0;
        MCOAngVel[2].setZero();
        MCOAngVel[2][2] = -M_PI / 2.0;
        MCOAngVel[3].setZero();
        MCOAngVel[3][2] = M_PI / 2.0;

        break;
    }

    case AST_DCOROTCYLINDERS: {
        assert(mesh.componentNodeRange.size() >= 6);

        mesh.resetFixedVert();
        std::vector<int> fixedVerts;
        MCORotCenter.resize(4);
        for (int compI = 0; compI < 4; ++compI) {
            for (int i = mesh.componentNodeRange[compI]; i < mesh.componentNodeRange[compI + 1]; ++i) {
                fixedVerts.emplace_back(i);
            }

            MCORotCenter[compI] = 0.5 * (mesh.V.block(mesh.componentNodeRange[compI], 0, mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI], 3).colwise().maxCoeff() + mesh.V.block(mesh.componentNodeRange[compI], 0, mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI], 3).colwise().minCoeff());
        }
        mesh.addFixedVert(fixedVerts);

        MCOAngVel.resize(4);
        MCOAngVel[0].setZero();
        MCOAngVel[0][0] = M_PI / 2.0;
        MCOAngVel[1].setZero();
        MCOAngVel[1][0] = -M_PI / 2.0;
        MCOAngVel[2].setZero();
        MCOAngVel[2][2] = -M_PI / 2.0;
        MCOAngVel[3].setZero();
        MCOAngVel[3][2] = M_PI / 2.0;

        break;
    }

    case AST_DCOVERSCHOORROLLER: {
        assert(mesh.componentNodeRange.size() >= 8);

        mesh.resetFixedVert();
        std::vector<int> fixedVerts;
        MCORotCenter.resize(6);
        for (int compI = 0; compI < 6; ++compI) {
            for (int i = mesh.componentNodeRange[compI]; i < mesh.componentNodeRange[compI + 1]; ++i) {
                fixedVerts.emplace_back(i);
            }

            MCORotCenter[compI] = 0.5 * (mesh.V.block(mesh.componentNodeRange[compI], 0, mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI], 3).colwise().maxCoeff() + mesh.V.block(mesh.componentNodeRange[compI], 0, mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI], 3).colwise().minCoeff());
        }
        mesh.addFixedVert(fixedVerts);

        MCOAngVel.resize(6);
        MCOAngVel[0].setZero();
        MCOAngVel[0][2] = -4.0;
        MCOAngVel[1].setZero();
        MCOAngVel[1][2] = -2.0;
        MCOAngVel[2].setZero();
        MCOAngVel[2][2] = 2.0;
        MCOAngVel[3].setZero();
        MCOAngVel[3][2] = 4.0;
        MCOAngVel[4].setZero();
        MCOAngVel[4][0] = 2.0;
        MCOAngVel[5].setZero();
        MCOAngVel[5][0] = -2.0;

        break;
    }

    case AST_DCOHAMMERWALNUT: {
        assert(mesh.componentCoDim.size() >= 2);

        mesh.resetFixedVert();
        std::vector<int> fixedVerts;
        Eigen::Matrix<double, 1, dim> leftBottom, rightTop;
        leftBottom.setConstant(__DBL_MAX__);
        rightTop.setConstant(-__DBL_MAX__);
        for (int i = mesh.componentNodeRange[1]; i < mesh.componentNodeRange[2]; ++i) {
            fixedVerts.emplace_back(i);

            leftBottom = leftBottom.array().min(mesh.V.row(i).array());
            rightTop = rightTop.array().max(mesh.V.row(i).array());
        }
        mesh.addFixedVert(fixedVerts);

        MCORotCenter.resize(1);
        MCORotCenter[0][0] = rightTop[0];
        MCORotCenter[0][1] = leftBottom[1];
        MCORotCenter[0][2] = (rightTop[2] + leftBottom[2]) / 2.0;

        MCOAngVel.resize(1);
        MCOAngVel[0].setZero();
        MCOAngVel[0][2] = M_PI / 6;

        break;
    }

    case AST_DCOCUT: {
        assert(mesh.componentCoDim.size() >= 2);

        mesh.resetFixedVert();
        std::vector<int> fixedVerts;
        Eigen::Matrix<double, 1, dim> leftBottom, rightTop;
        leftBottom.setConstant(__DBL_MAX__);
        rightTop.setConstant(-__DBL_MAX__);
        for (int i = mesh.componentNodeRange[1]; i < mesh.componentNodeRange[2]; ++i) {
            fixedVerts.emplace_back(i);

            leftBottom = leftBottom.array().min(mesh.V.row(i).array());
            rightTop = rightTop.array().max(mesh.V.row(i).array());
        }
        mesh.addFixedVert(fixedVerts);

        MCOVelocity.resize(1);
        MCOVelocity.back()[0] = 0;
        MCOVelocity.back()[1] = -1;
        MCOVelocity.back()[2] = -1;

        break;
    }

    case AST_DCOSQUASH: {
        assert(mesh.componentNodeRange.size() >= 4);

        mesh.resetFixedVert();
        std::vector<int> fixedVerts;
        for (int compI = 0; compI < 2; ++compI) {
            for (int i = mesh.componentNodeRange[compI]; i < mesh.componentNodeRange[compI + 1]; ++i) {
                fixedVerts.emplace_back(i);
            }
        }
        mesh.addFixedVert(fixedVerts);

        MCOVelocity.resize(2);
        MCOVelocity[0].setZero();
        MCOVelocity[0][0] = 1.0;
        MCOVelocity[1].setZero();
        MCOVelocity[1][0] = -1.0;

        break;
    }

    case AST_DCOSQUASH6: {
        assert(mesh.componentNodeRange.size() >= 8);

        mesh.resetFixedVert();
        std::vector<int> fixedVerts;
        for (int compI = 0; compI < 6; ++compI) {
            for (int i = mesh.componentNodeRange[compI]; i < mesh.componentNodeRange[compI + 1]; ++i) {
                fixedVerts.emplace_back(i);
            }
        }
        mesh.addFixedVert(fixedVerts);

        MCOVelocity.resize(6);
        MCOVelocity[0].setZero();
        MCOVelocity[0][0] = 1.0;
        MCOVelocity[1].setZero();
        MCOVelocity[1][0] = -1.0;
        MCOVelocity[2].setZero();
        MCOVelocity[2][1] = 1.0;
        MCOVelocity[3].setZero();
        MCOVelocity[3][1] = -1.0;
        MCOVelocity[4].setZero();
        MCOVelocity[4][2] = 1.0;
        MCOVelocity[5].setZero();
        MCOVelocity[5][2] = -1.0;

        break;
    }

    case AST_MESHSEQ_FROMFILE:
    case AST_DCOBALLHITWALL:
    case AST_DCOFIX: {
        mesh.resetFixedVert();
        std::vector<int> fixedVerts;
        for (int compI = 0; compI < mesh.componentCoDim.size(); ++compI) {
            if (mesh.componentCoDim[compI] < 3) {
                for (int i = mesh.componentNodeRange[compI]; i < mesh.componentNodeRange[compI + 1]; ++i) {
                    fixedVerts.emplace_back(i);
                }
            }
        }
        mesh.addFixedVert(fixedVerts);

        break;
    }

    case AST_DCOSEGBEDSQUASH: {
        mesh.resetFixedVert();
        std::vector<int> fixedVerts;
        for (int compI = 0; compI < mesh.componentCoDim.size(); ++compI) {
            if (mesh.componentCoDim[compI] < 3) {
                for (int i = mesh.componentNodeRange[compI]; i < mesh.componentNodeRange[compI + 1]; ++i) {
                    fixedVerts.emplace_back(i);
                }

                if (compI >= (mesh.componentCoDim.size() + 1) / 2) {
                    for (int i = mesh.componentNodeRange[compI]; i < mesh.componentNodeRange[compI + 1]; ++i) {
                        velocity_handleVerts[i].setZero();
                        velocity_handleVerts[i][1] = -1.0;
                    }
                }
            }
        }
        mesh.addFixedVert(fixedVerts);

        break;
    }

    case AST_DCOSQUEEZEOUT: {
        mesh.resetFixedVert();
        std::vector<int> fixedVerts;
        for (int compI = 0; compI < mesh.componentCoDim.size(); ++compI) {
            if (mesh.componentCoDim[compI] < 3) {
                for (int i = mesh.componentNodeRange[compI]; i < mesh.componentNodeRange[compI + 1]; ++i) {
                    fixedVerts.emplace_back(i);
                }

                if (compI == 0) {
                    for (int i = mesh.componentNodeRange[compI]; i < mesh.componentNodeRange[compI + 1]; ++i) {
                        velocity_handleVerts[i].setZero();
                        velocity_handleVerts[i][1] = -0.3;
                    }
                }
            }
        }
        mesh.addFixedVert(fixedVerts);

        break;
    }

    case AST_UTOPIA_COMPARISON: {
        mesh.resetFixedVert();
        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        // fix bottom and push top down
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            bool isOnBottom = mesh.V(vI, 1) < bottomLeft[1] + range[0] * 1e-4;
            bool isOnTop = mesh.V(vI, 1) > topRight[1] - range[0] * 1e-4;
            if (isOnTop) {
                mesh.NeumannBC[vI].setZero();
                mesh.NeumannBC[vI][1] = -1.5;
            }
            else if (isOnBottom) {
                mesh.addFixedVert(vI);
            }
        }
        break;
    }

    default:
        spdlog::error("Invalid animScriptType: {:s}",
            getStrByAnimScriptType(animScriptType));
        assert(0 && "invalid animScriptType");
        break;
    }

    fixedVertBK = mesh.fixedVert;
}

template <int dim>
void AnimScripter<dim>::initVelocity(const Mesh<dim>& mesh,
    const std::vector<double>& params, Eigen::VectorXd& velocity) const
{
    switch (animScriptType) {
    case AST_NULL: {
        velocity.setZero(mesh.V.rows() * dim);
        for (const auto& initVelInfo : mesh.componentInitVels) {
            int compI = initVelInfo.first[0];
            if (mesh.componentCoDim[compI] == dim) {
                Eigen::Vector3d rotCenter = 0.5 * (mesh.V.block(mesh.componentNodeRange[compI], 0, mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI], 3).colwise().maxCoeff() + mesh.V.block(mesh.componentNodeRange[compI], 0, mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI], 3).colwise().minCoeff());

                for (int vI = initVelInfo.first[1]; vI < initVelInfo.first[2]; ++vI) {
                    if (!mesh.isFixedVert[vI]) {
                        velocity.template segment<dim>(vI * dim) = initVelInfo.second[0] + initVelInfo.second[1].cross(Eigen::Vector3d(mesh.V.row(vI).transpose()) - rotCenter);
                    }
                }
            }
        }
        break;
    }

    case AST_LEFTHITRIGHT: {
        velocity.setZero(mesh.V.rows() * dim);

        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 0) < bottomLeft[0] + range[0] / 2.0) {
                velocity[vI * dim] = 1.0;
            }
        }

        break;
    }

    case AST_DROP: {
        velocity.setZero(mesh.V.rows() * dim);
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            velocity[vI * dim + 1] = -1.0;
        }
        break;
    }

    case AST_XYROTATE: {
        velocity.setZero(mesh.V.rows() * dim);

        Eigen::RowVectorXd bottomLeft = mesh.V.colwise().minCoeff();
        Eigen::RowVectorXd topRight = mesh.V.colwise().maxCoeff();
        Eigen::RowVectorXd range = topRight - bottomLeft;
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.V(vI, 1) < bottomLeft[1] + range[1] * 0.01) {
                velocity[vI * dim] = 1.0;
            }
            else if (mesh.V(vI, 1) > topRight[1] - range[1] * 0.01) {
                velocity[vI * dim] = -1.0;
            }
        }
        break;
    }

    case AST_DCOBALLHITWALL: {
        double v_x = 1000.0;
        if (params.size() > 0 && params[0] == params[0]) {
            v_x = params[0];
        }

        velocity.setZero(mesh.V.rows() * dim);
        for (int vI = 0; vI < mesh.V.rows(); ++vI) {
            if (mesh.vICoDim(vI) == 3) {
                velocity[vI * dim] = v_x;
            }
        }
        break;
    }

    default:
        velocity.setZero(mesh.V.rows() * dim);
        break;
    }
}

template <int dim>
int AnimScripter<dim>::stepAnimScript(Mesh<dim>& mesh,
    const SpatialHash<dim>& sh,
    std::vector<CollisionObject<dim>*>& ACO,
    std::vector<CollisionObject<dim>*>& MCO,
    const ExactCCD::Method method,
    double dt, double dHat, const std::vector<Energy<dim>*>& energyTerms,
    bool isSelfCollision, bool forceIntersectionLineSearch)
{
    curTime += dt;

    searchDir.setZero(mesh.V.rows() * dim);
    int returnFlag = 0;
    switch (animScriptType) {
    case AST_NULL: {
        // move objects with scripted angular velocity
        for (const auto& scriptedAVelI : mesh.componentAVels) {
            Eigen::Matrix3d rotMtr = (Eigen::AngleAxisd(scriptedAVelI.second[0] * dt, Eigen::Vector3d::UnitX())
                * Eigen::AngleAxisd(scriptedAVelI.second[1] * dt, Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(scriptedAVelI.second[2] * dt, Eigen::Vector3d::UnitZ()))
                                         .toRotationMatrix();

            int compI = scriptedAVelI.first[0];
            Eigen::Vector3d MCORotCenter = 0.5 * (mesh.V.block(mesh.componentNodeRange[compI], 0, mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI], 3).colwise().maxCoeff() + mesh.V.block(mesh.componentNodeRange[compI], 0, mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI], 3).colwise().minCoeff());

            for (int vI = scriptedAVelI.first[1]; vI < scriptedAVelI.first[2]; ++vI) {
                searchDir.template segment<dim>(vI * dim)
                    = rotMtr.template block<dim, dim>(0, 0) * (mesh.V.row(vI).transpose() - MCORotCenter) + MCORotCenter - mesh.V.row(vI).transpose();
            }
        }

        // move objects with scripted linear velocity
        for (const auto& scriptedLVelI : mesh.componentLVels) {
            for (int i = scriptedLVelI.first[1]; i < scriptedLVelI.first[2]; ++i) {
                searchDir.template segment<dim>(i * dim) += scriptedLVelI.second * dt;
            }
        }

        // move Dirichlet nodes
        if (curTime > DBCTimeRange[0] && curTime <= DBCTimeRange[1]) {
            if (mesh.fixedVert.empty()) {
                mesh.resetFixedVert(fixedVertBK);
            }
            for (const auto& DBCInfoI : mesh.DBCInfo) {
                Eigen::Matrix3d rotMtr = (Eigen::AngleAxisd(DBCInfoI.second[1][0] * dt, Eigen::Vector3d::UnitX())
                    * Eigen::AngleAxisd(DBCInfoI.second[1][1] * dt, Eigen::Vector3d::UnitY())
                    * Eigen::AngleAxisd(DBCInfoI.second[1][2] * dt, Eigen::Vector3d::UnitZ()))
                                             .toRotationMatrix();

                Eigen::Vector3d min = mesh.V.row(DBCInfoI.first[0]), max = mesh.V.row(DBCInfoI.first[0]);
                for (const auto& vI : DBCInfoI.first) {
                    for (int d = 0; d < 3; ++d) {
                        if (mesh.V(vI, d) < min[d]) {
                            min[d] = mesh.V(vI, d);
                        }
                        if (mesh.V(vI, d) > max[d]) {
                            max[d] = mesh.V(vI, d);
                        }
                    }
                }
                Eigen::Vector3d rotCenter = (min + max) / 2;

                for (const auto& vI : DBCInfoI.first) {
                    searchDir.template segment<dim>(vI * dim)
                        = rotMtr.template block<dim, dim>(0, 0) * (mesh.V.row(vI).transpose() - rotCenter) + rotCenter + DBCInfoI.second[0] * dt - mesh.V.row(vI).transpose();
                }
            }
        }
        else {
            mesh.resetFixedVert();
        }

        // mesh sequence DBC
        for (const auto& meshSeqI : mesh.meshSeqFolderPath) {
            int compI = meshSeqI.first;
            Eigen::MatrixXd Vms;
            switch (mesh.componentCoDim[compI]) {
            case 3: {
                // msh
                Eigen::MatrixXi Fms, SFms;
                if (!IglUtils::readTetMesh(meshSeqI.second + "/" + std::to_string(meshI) + ".msh", Vms, Fms, SFms)) {
                    spdlog::error("Unable to read input file: {:s}", meshSeqI.second + "/" + std::to_string(meshI) + ".msh");
                    exit(-1);
                }
                break;
            }

            case 2: {
                // obj
                Eigen::MatrixXi Fms;
                if (!igl::readOBJ(meshSeqI.second + "/" + std::to_string(meshI) + ".obj", Vms, Fms)) {
                    spdlog::error("Unable to read input file: {:s}", meshSeqI.second + "/" + std::to_string(meshI) + ".obj");
                    exit(-1);
                }
                break;
            }

            case 1: {
                // seg
                Eigen::MatrixXi Ems;
                if (!IglUtils::readSEG(meshSeqI.second + "/" + std::to_string(meshI) + ".seg", Vms, Ems)) {
                    // if not found, from obj
                    Eigen::MatrixXi Fms;
                    if (!igl::readOBJ(meshSeqI.second + "/" + std::to_string(meshI) + ".obj", Vms, Fms)) {
                        spdlog::error("Unable to read input file: {:s}", meshSeqI.second + "/" + std::to_string(meshI) + ".seg or .obj");
                        exit(-1);
                    }
                }
                break;
            }

            case 0: {
                // pt
                Eigen::MatrixXi Fms;
                if (!igl::readOBJ(meshSeqI.second + "/" + std::to_string(meshI) + ".pt", Vms, Fms)) {
                    // if not found, from obj
                    if (!igl::readOBJ(meshSeqI.second + "/" + std::to_string(meshI) + ".obj", Vms, Fms)) {
                        spdlog::error("Unable to read input file: {:s}", meshSeqI.second + "/" + std::to_string(meshI) + ".pt or .obj");
                        exit(-1);
                    }
                }
                break;
            }

            default: {
                break;
            }
            }

            if (Vms.rows() != mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI]) {
                spdlog::warn("meshSeq #V does not match with mesh!");
                exit(-1);
            }
            for (int vI = mesh.componentNodeRange[compI]; vI < mesh.componentNodeRange[compI + 1]; ++vI) {
                searchDir.template segment<dim>(vI * dim) = (Vms.row(vI - mesh.componentNodeRange[compI]) - mesh.V.row(vI)).transpose();
            }
        }
        if (mesh.meshSeqFolderPath.size()) {
            ++meshI;
        }
        break;
    }

    case AST_HANG:
    case AST_HANG2:
    case AST_HANGTOPLEFT:
    case AST_HANGLEFT:
        break;

    case AST_STAMP:
    case AST_STAMPTOPLEFT:
    case AST_STAMPBOTH:
    case AST_STAMPINV:
        break;

    case AST_STAND:
    case AST_STANDINV:
    case AST_TOPBOTTOMFIX:
    case AST_FIXLOWERHALF:
        break;

    case AST_PUSH: {
        if (mesh.V(velocityTurningPoints.first, 1) <= velocityTurningPoints.second(0, 1)) {
            velocityTurningPoints.second(0, 1) = -__DBL_MAX__;
            for (const auto& vI : handleVerts[0]) {
                velocity_handleVerts[vI].setZero();
            }
        }

        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_TEAR: {
        if (mesh.V(velocityTurningPoints.first, 0) <= velocityTurningPoints.second(0, 0)) {
            velocityTurningPoints.second(0, 1) = -__DBL_MAX__;
            for (const auto& vI : handleVerts[0]) {
                velocity_handleVerts[vI] = -velocity_handleVerts[vI];
            }
        }

        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_UNDSTAMP:
    case AST_UPNDOWN: {
        bool flip = false;
        if ((mesh.V(velocityTurningPoints.first, 1) <= velocityTurningPoints.second(1, 0)) || (mesh.V(velocityTurningPoints.first, 1) >= velocityTurningPoints.second(1, 1))) {
            flip = true;
        }
        for (auto& movingVerts : velocity_handleVerts) {
            if (flip) {
                movingVerts.second[1] *= -1.0;
            }
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_STRETCH:
    case AST_SQUASH:
    case AST_DRAGDOWN:
        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;

    case AST_STRETCHNPAUSE: {
        if (mesh.V(velocityTurningPoints.first, 0) >= velocityTurningPoints.second(0, 0)) {
            for (const auto& movingVerts : velocity_handleVerts) {
                searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
            }
        }
        break;
    }

    case AST_DRAGRIGHT: {
        double MCORightMost = -__DBL_MAX__;
        for (const auto& mcoI : MCO) {
            double maxI = mcoI->V.col(0).maxCoeff();
            if (MCORightMost < maxI) {
                MCORightMost = maxI;
            }
        }
        if (mesh.V.col(0).minCoeff() > MCORightMost && mesh.fixedVert.size()) {
            mesh.resetFixedVert();
            for (auto& movingVerts : velocity_handleVerts) {
                movingVerts.second.setZero();
            }
            returnFlag = 1;
        }

        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_UTOPIA_COMPARISON: {
        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_STRETCHNSQUASH: {
        bool flip = false;
        if ((mesh.V(velocityTurningPoints.first, 0) <= velocityTurningPoints.second(0, 0)) || (mesh.V(velocityTurningPoints.first, 0) >= velocityTurningPoints.second(0, 1))) {
            flip = true;
        }
        for (auto& movingVerts : velocity_handleVerts) {
            if (flip) {
                movingVerts.second[0] *= -1.0;
            }
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_BEND:
        for (const auto& movingVerts : angVel_handleVerts) {
            const Eigen::Matrix3d rotMtr = Eigen::AngleAxis<double>(movingVerts.second * dt,
                Eigen::Vector3d::UnitZ())
                                               .toRotationMatrix();
            const auto rotCenter = rotCenter_handleVerts.find(movingVerts.first);
            assert(rotCenter != rotCenter_handleVerts.end());

            searchDir.segment<dim>(movingVerts.first * dim) = (rotMtr.block<dim, dim>(0, 0) * (mesh.V.row(movingVerts.first).transpose() - rotCenter->second) + rotCenter->second) - mesh.V.row(movingVerts.first).transpose();
        }
        break;

    case AST_TWIST:
        for (const auto& movingVerts : angVel_handleVerts) {
            const Eigen::Matrix3d rotMtr = Eigen::AngleAxis<double>(movingVerts.second * dt,
                Eigen::Vector3d::UnitX())
                                               .toRotationMatrix();
            const auto rotCenter = rotCenter_handleVerts.find(movingVerts.first);
            assert(rotCenter != rotCenter_handleVerts.end());

            searchDir.segment<dim>(movingVerts.first * dim) = (rotMtr.block<dim, dim>(0, 0) * (mesh.V.row(movingVerts.first).transpose() - rotCenter->second) + rotCenter->second) - mesh.V.row(movingVerts.first).transpose();
        }
        break;

    case AST_TWISTNSTRETCH: {
        for (const auto& movingVerts : angVel_handleVerts) {
            const Eigen::Matrix3d rotMtr = Eigen::AngleAxis<double>(movingVerts.second * dt,
                Eigen::Vector3d::UnitX())
                                               .toRotationMatrix();
            const auto rotCenter = rotCenter_handleVerts.find(movingVerts.first);
            assert(rotCenter != rotCenter_handleVerts.end());

            searchDir.segment<dim>(movingVerts.first * dim) = (rotMtr.block<dim, dim>(0, 0) * (mesh.V.row(movingVerts.first).transpose() - rotCenter->second) + rotCenter->second) - mesh.V.row(movingVerts.first).transpose();
        }
        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) += movingVerts.second * dt;
        }
        break;
    }

    case AST_TWISTNSNS_OLD:
    case AST_TWISTNSNS: {
        bool flip = false;
        if ((mesh.V(velocityTurningPoints.first, 0) <= velocityTurningPoints.second(0, 0)) || (mesh.V(velocityTurningPoints.first, 0) >= velocityTurningPoints.second(0, 1))) {
            flip = true;
        }

        for (auto& movingVerts : angVel_handleVerts) {
            //                    if(flip) {
            //                        movingVerts.second *= -1.0;
            //                    }

            const Eigen::Matrix3d rotMtr = Eigen::AngleAxis<double>(movingVerts.second * dt,
                Eigen::Vector3d::UnitX())
                                               .toRotationMatrix();
            const auto rotCenter = rotCenter_handleVerts.find(movingVerts.first);
            assert(rotCenter != rotCenter_handleVerts.end());

            searchDir.segment<dim>(movingVerts.first * dim) = (rotMtr.block<dim, dim>(0, 0) * (mesh.V.row(movingVerts.first).transpose() - rotCenter->second) + rotCenter->second) - mesh.V.row(movingVerts.first).transpose();
        }
        for (auto& movingVerts : velocity_handleVerts) {
            if (flip) {
                movingVerts.second[0] *= -1.0;
            }
            searchDir.segment<dim>(movingVerts.first * dim) += movingVerts.second * dt;
        }
        break;
    }

    case AST_RUBBERBANDPULL: {
        if (mesh.V(velocityTurningPoints.first, 0) <= velocityTurningPoints.second(0, 0)) {
            velocityTurningPoints.second(0, 0) = -__DBL_MAX__;
            for (const auto& vI : handleVerts[0]) {
                mesh.removeFixedVert(vI);
                velocity_handleVerts[vI].setZero();
            }
            for (const auto& vI : handleVerts[1]) {
                velocity_handleVerts[vI].setZero();
            }
            returnFlag = 1;
        }
        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_FOURLEGPULL: {
        if (mesh.V(velocityTurningPoints.first, 1) <= velocityTurningPoints.second(0, 1)) {
            velocityTurningPoints.second(0, 1) = -__DBL_MAX__;
            for (const auto& vI : handleVerts[0]) {
                mesh.removeFixedVert(vI);
                velocity_handleVerts[vI].setZero();
            }
            for (const auto& vI : handleVerts[1]) {
                velocity_handleVerts[vI].setZero();
            }
            returnFlag = 1;
        }

        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_HEADTAILPULL: {
        if (mesh.V(velocityTurningPoints.first, 0) >= velocityTurningPoints.second(0, 0)) {
            velocityTurningPoints.second(0, 0) = __DBL_MAX__;
            for (const auto& vI : handleVerts[0]) {
                mesh.removeFixedVert(vI);
                velocity_handleVerts[vI].setZero();
            }
            for (const auto& vI : handleVerts[1]) {
                velocity_handleVerts[vI].setZero();
            }
            returnFlag = 1;
        }

        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_ONEPOINT:
        break;

    case AST_RANDOM:
        break;

    case AST_FALL:
        break;

    case AST_FALL_NOSHIFT:
        break;

    case AST_TOGGLETOP: {
        if (mesh.V(velocityTurningPoints.first, 0) <= velocityTurningPoints.second(0, 0)) {
            velocityTurningPoints.second(0, 0) = -__DBL_MAX__;
            for (const auto& vI : handleVerts[0]) {
                mesh.removeFixedVert(vI);
                velocity_handleVerts[vI].setZero();
            }
            returnFlag = 1;
        }
        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_LEFTHITRIGHT:
    case AST_DROP:
    case AST_XYROTATE:
    case AST_FIXRIGHTMOST1:
        break;

    case AST_PUSHRIGHTMOST1:
    case AST_CURTAIN: {
        for (const auto& movingVerts : velocity_handleVerts) {
            searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
        }
        break;
    }

    case AST_NMFIXBOTTOMDRAGLEFT:
    case AST_NMFIXBOTTOMDRAGFORWARD:
        break;

    case AST_ACOSQUASH: {
        if (ACO.size() < 2) {
            throw "not enough ACO for ACOSquash script";
        }

        if (ACO[1]->origin[0] - ACO[0]->origin[0] < 0.1) {
            ACOVelocity[0][0] *= -1.0;
            ACOVelocity[1][0] *= -1.0;
        }

        ACOStepSizeLeft.resize(2);
        ACO[0]->move(ACOVelocity[0] * dt, mesh, sh, 0.5, ACOStepSizeLeft[0]);
        ACO[1]->move(ACOVelocity[1] * dt, mesh, sh, 0.5, ACOStepSizeLeft[1]);
        break;
    }

    case AST_ACOSQUASH6: {
        if (ACO.size() < 6) {
            throw "not enough ACO for ACOSquash6 script";
        }

        if (ACO[1]->origin[0] - ACO[0]->origin[0] < 0.2) {
            ACOVelocity[0][0] *= -1.0;
            ACOVelocity[1][0] *= -1.0;
        }
        if (ACO[3]->origin[1] - ACO[2]->origin[1] < 0.2) {
            ACOVelocity[2][1] *= -1.0;
            ACOVelocity[3][1] *= -1.0;
        }
        if (ACO[5]->origin[2] - ACO[4]->origin[2] < 0.2) {
            ACOVelocity[4][2] *= -1.0;
            ACOVelocity[5][2] *= -1.0;
        }

        ACOStepSizeLeft.resize(6);
        for (int acoI = 0; acoI < ACOStepSizeLeft.size(); ++acoI) {
            ACO[acoI]->move(ACOVelocity[acoI] * dt, mesh, sh, 0.5, ACOStepSizeLeft[acoI]);
        }
        break;
    }

    case AST_ACOSQUASHSHEAR: {
        if (ACO.size() < 2) {
            throw "not enough ACO for ACOSquashShear script";
        }

        if (ACO[1]->origin[0] - ACO[0]->origin[0] < 0.8) {
            ACOVelocity[0][0] = 0.0;
            ACOVelocity[1][0] = 0.0;
            ACOVelocity[1][1] = 1.0;
        }

        ACOStepSizeLeft.resize(2);
        ACO[0]->velocitydt = ACOVelocity[0] * dt;
        ACO[0]->move(ACO[0]->velocitydt, mesh, sh, 0.5, ACOStepSizeLeft[0]);
        ACO[1]->velocitydt = ACOVelocity[1] * dt;
        ACO[1]->move(ACO[1]->velocitydt, mesh, sh, 0.5, ACOStepSizeLeft[1]);
        break;
    }

    case AST_MCOSQUASH: {
        if (MCO.size() < 2) {
            throw "not enough MCO for MCOSquash script";
        }

        if (MCO[1]->origin[0] - MCO[0]->origin[0] < 1.5) {
            MCOVelocity[0][0] *= -1.0;
            MCOVelocity[1][0] *= -1.0;
        }

        MCO[0]->V_target = MCO[0]->V;
        MCO[0]->V_target.rowwise() += MCOVelocity[0].transpose() * dt;
        MCO[1]->V_target = MCO[1]->V;
        MCO[1]->V_target.rowwise() += MCOVelocity[1].transpose() * dt;

        MCOStepSizeLeft.resize(2);
        MCO[0]->move(MCOVelocity[0] * dt, mesh, sh, 0.5, MCOStepSizeLeft[0]);
        MCO[1]->move(MCOVelocity[1] * dt, mesh, sh, 0.5, MCOStepSizeLeft[1]);
        break;
    }

    case AST_MCOROTSQUASH: {
        if (MCO.size() < 2) {
            throw "not enough MCO for MCORotSquash script";
        }

        for (int coI = 0; coI < MCO.size(); ++coI) {
            MCO[coI]->V_target.resize(MCO[coI]->V.rows(), MCO[coI]->V.cols());
            Eigen::Matrix3d rotMtr = (Eigen::AngleAxisd(MCOAngVel[coI][0] * dt, Eigen::Vector3d::UnitX())
                * Eigen::AngleAxisd(MCOAngVel[coI][1] * dt, Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(MCOAngVel[coI][2] * dt, Eigen::Vector3d::UnitZ()))
                                         .toRotationMatrix();
            for (int vI = 0; vI < MCO[coI]->V.rows(); ++vI) {
                MCO[coI]->V_target.row(vI) = (rotMtr.template block<dim, dim>(0, 0) * (MCO[coI]->V.row(vI).transpose() - MCORotCenter[coI]) + MCORotCenter[coI]).transpose();
            }
            MCOVelocity[coI] = (rotMtr.template block<dim, dim>(0, 0) * (MCO[coI]->origin - MCORotCenter[coI]) + MCORotCenter[coI] - MCO[coI]->origin) / dt;
        }

        MCOStepSizeLeft.resize(2);
        MCO[0]->move(MCOVelocity[0] * dt, mesh, sh, 0.5, MCOStepSizeLeft[0]);
        MCO[1]->move(MCOVelocity[1] * dt, mesh, sh, 0.5, MCOStepSizeLeft[1]);
        break;
    }

    case AST_MCOROTCYLINDERS: {
        if (MCO.size() < 4) {
            throw "not enough MCO for MCOROTCYLINDERS script";
        }

        for (int coI = 0; coI < MCO.size(); ++coI) {
            MCO[coI]->V_target.resize(MCO[coI]->V.rows(), MCO[coI]->V.cols());
            Eigen::Matrix3d rotMtr = (Eigen::AngleAxisd(MCOAngVel[coI][0] * dt, Eigen::Vector3d::UnitX())
                * Eigen::AngleAxisd(MCOAngVel[coI][1] * dt, Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(MCOAngVel[coI][2] * dt, Eigen::Vector3d::UnitZ()))
                                         .toRotationMatrix();
            for (int vI = 0; vI < MCO[coI]->V.rows(); ++vI) {
                MCO[coI]->V_target.row(vI) = (rotMtr.template block<dim, dim>(0, 0) * (MCO[coI]->V.row(vI).transpose() - MCORotCenter[coI]) + MCORotCenter[coI]).transpose();
            }
            MCOVelocity[coI] = (rotMtr.template block<dim, dim>(0, 0) * (MCO[coI]->origin - MCORotCenter[coI]) + MCORotCenter[coI] - MCO[coI]->origin) / dt;
        }

        MCOStepSizeLeft.resize(4);
        MCO[0]->move(MCOVelocity[0] * dt, mesh, sh, 0.5, MCOStepSizeLeft[0]);
        MCO[1]->move(MCOVelocity[1] * dt, mesh, sh, 0.5, MCOStepSizeLeft[1]);
        MCO[2]->move(MCOVelocity[2] * dt, mesh, sh, 0.5, MCOStepSizeLeft[2]);
        MCO[3]->move(MCOVelocity[3] * dt, mesh, sh, 0.5, MCOStepSizeLeft[3]);
        break;
    }

    case AST_DCOROTCYLINDERS: {
        assert(mesh.componentNodeRange.size() >= 6);

        for (int compI = 0; compI < 4; ++compI) {
            Eigen::Matrix3d rotMtr = (Eigen::AngleAxisd(MCOAngVel[compI][0] * dt, Eigen::Vector3d::UnitX())
                * Eigen::AngleAxisd(MCOAngVel[compI][1] * dt, Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(MCOAngVel[compI][2] * dt, Eigen::Vector3d::UnitZ()))
                                         .toRotationMatrix();
            for (int vI = mesh.componentNodeRange[compI]; vI < mesh.componentNodeRange[compI + 1]; ++vI) {
                searchDir.template segment<dim>(vI * dim)
                    = rotMtr.template block<dim, dim>(0, 0) * (mesh.V.row(vI).transpose() - MCORotCenter[compI]) + MCORotCenter[compI] - mesh.V.row(vI).transpose();
            }
        }
        break;
    }

    case AST_DCOVERSCHOORROLLER: {
        assert(mesh.componentNodeRange.size() >= 8);

        for (int compI = 0; compI < 6; ++compI) {
            Eigen::Matrix3d rotMtr = (Eigen::AngleAxisd(MCOAngVel[compI][0] * dt, Eigen::Vector3d::UnitX())
                * Eigen::AngleAxisd(MCOAngVel[compI][1] * dt, Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(MCOAngVel[compI][2] * dt, Eigen::Vector3d::UnitZ()))
                                         .toRotationMatrix();
            for (int vI = mesh.componentNodeRange[compI]; vI < mesh.componentNodeRange[compI + 1]; ++vI) {
                searchDir.template segment<dim>(vI * dim)
                    = rotMtr.template block<dim, dim>(0, 0) * (mesh.V.row(vI).transpose() - MCORotCenter[compI]) + MCORotCenter[compI] - mesh.V.row(vI).transpose();
            }
        }
        break;
    }

    case AST_DCOHAMMERWALNUT: {
        assert(mesh.componentCoDim.size() >= 2);

        double yMin = __DBL_MAX__;
        for (int vI = mesh.componentNodeRange[1]; vI < mesh.componentNodeRange[2]; ++vI) {
            if (yMin > mesh.V(vI, 1)) {
                yMin = mesh.V(vI, 1);
            }
        }
        if (yMin > 0.05) {
            Eigen::Matrix3d rotMtr = (Eigen::AngleAxisd(MCOAngVel[0][0] * dt, Eigen::Vector3d::UnitX())
                * Eigen::AngleAxisd(MCOAngVel[0][1] * dt, Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(MCOAngVel[0][2] * dt, Eigen::Vector3d::UnitZ()))
                                         .toRotationMatrix();
            for (int vI = mesh.componentNodeRange[1]; vI < mesh.componentNodeRange[2]; ++vI) {
                searchDir.template segment<dim>(vI * dim)
                    = rotMtr.template block<dim, dim>(0, 0) * (mesh.V.row(vI).transpose() - MCORotCenter[0]) + MCORotCenter[0] - mesh.V.row(vI).transpose();
            }
        }

        break;
    }

    case AST_DCOCUT: {
        assert(mesh.componentCoDim.size() >= 2);

        double yMin = __DBL_MAX__;
        for (int vI = mesh.componentNodeRange[1]; vI < mesh.componentNodeRange[2]; ++vI) {
            if (yMin > mesh.V(vI, 1)) {
                yMin = mesh.V(vI, 1);
            }
        }
        if (yMin > 0.001) {
            for (int vI = mesh.componentNodeRange[1]; vI < mesh.componentNodeRange[2]; ++vI) {
                searchDir.template segment<dim>(vI * dim) = MCOVelocity[0] * dt;
            }
        }

        break;
    }

    case AST_DCOSQUASH: {
        assert(mesh.componentNodeRange.size() >= 4);

        double co0right = mesh.V.col(0).segment(mesh.componentNodeRange[0], mesh.componentNodeRange[1] - mesh.componentNodeRange[0]).maxCoeff();
        double co1left = mesh.V.col(0).segment(mesh.componentNodeRange[1], mesh.componentNodeRange[2] - mesh.componentNodeRange[1]).minCoeff();
        if (co1left - co0right < 0.1) {
            MCOVelocity[0][0] *= -1.0;
            MCOVelocity[1][0] *= -1.0;
        }

        for (int compI = 0; compI < 2; ++compI) {
            for (int vI = mesh.componentNodeRange[compI]; vI < mesh.componentNodeRange[compI + 1]; ++vI) {
                searchDir.template segment<dim>(vI * dim) = MCOVelocity[compI] * dt;
            }
        }

        break;
    }

    case AST_DCOSQUASH6: {
        assert(mesh.componentNodeRange.size() >= 8);

        double co0right = mesh.V.col(0).segment(mesh.componentNodeRange[0], mesh.componentNodeRange[1] - mesh.componentNodeRange[0]).maxCoeff();
        double co1left = mesh.V.col(0).segment(mesh.componentNodeRange[1], mesh.componentNodeRange[2] - mesh.componentNodeRange[1]).minCoeff();
        if (co1left - co0right < 0.1) {
            MCOVelocity[0][0] *= -1.0;
            MCOVelocity[1][0] *= -1.0;
            MCOVelocity[2][1] *= -1.0;
            MCOVelocity[3][1] *= -1.0;
            MCOVelocity[4][2] *= -1.0;
            MCOVelocity[5][2] *= -1.0;
        }

        for (int compI = 0; compI < 6; ++compI) {
            for (int vI = mesh.componentNodeRange[compI]; vI < mesh.componentNodeRange[compI + 1]; ++vI) {
                searchDir.template segment<dim>(vI * dim) = MCOVelocity[compI] * dt;
            }
        }

        break;
    }

    case AST_DCOBALLHITWALL:
    case AST_DCOFIX:
        break;

    case AST_DCOSEGBEDSQUASH: {
        double topMin = __DBL_MAX__, bottomMax = -__DBL_MAX__;
        for (int compI = 0; compI < mesh.componentCoDim.size(); ++compI) {
            if (mesh.componentCoDim[compI] < 3) {
                if (compI >= (mesh.componentCoDim.size() + 1) / 2) {
                    topMin = std::min(topMin, mesh.V.col(1).segment(mesh.componentNodeRange[compI], mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI]).minCoeff());
                }
                else {
                    bottomMax = std::max(bottomMax, mesh.V.col(1).segment(mesh.componentNodeRange[compI], mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI]).maxCoeff());
                }
            }
        }

        if (topMin - bottomMax > 0.1) {
            for (const auto& movingVerts : velocity_handleVerts) {
                searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
            }
        }

        break;
    }

    case AST_DCOSQUEEZEOUT: {
        double topMax = -__DBL_MAX__, bottomMax = -__DBL_MAX__, bottomMin = __DBL_MAX__;
        for (int compI = 0; compI < mesh.componentCoDim.size(); ++compI) {
            if (mesh.componentCoDim[compI] < 3) {
                if (compI == 0) {
                    topMax = std::max(topMax, mesh.V.col(1).segment(mesh.componentNodeRange[compI], mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI]).maxCoeff());
                }
                else if (compI == 1) {
                    bottomMin = std::min(bottomMin, mesh.V.col(1).segment(mesh.componentNodeRange[compI], mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI]).minCoeff());
                    bottomMax = std::max(bottomMax, mesh.V.col(1).segment(mesh.componentNodeRange[compI], mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI]).maxCoeff());
                }
            }
        }

        if (topMax > bottomMin + (bottomMax - bottomMin) / 3.8 * 0.9) {
            for (const auto& movingVerts : velocity_handleVerts) {
                searchDir.segment<dim>(movingVerts.first * dim) = movingVerts.second * dt;
            }
        }

        break;
    }

    case AST_MESHSEQ_FROMFILE: {
        ++meshI;
        Eigen::MatrixXd Vms;
        Eigen::MatrixXi Fms;
        igl::readOBJ(meshSeqFolderPath + "/" + std::to_string(meshI) + ".obj", Vms, Fms);
        for (int compI = 0; compI < mesh.componentCoDim.size(); ++compI) {
            if (mesh.componentCoDim[compI] < 3) {
                if (Vms.rows() != mesh.componentNodeRange[compI + 1] - mesh.componentNodeRange[compI]) {
                    spdlog::warn("meshSeq #V does not match with DCO!");
                    exit(-1);
                }
                for (int vI = mesh.componentNodeRange[compI]; vI < mesh.componentNodeRange[compI + 1]; ++vI) {
                    searchDir.segment<dim>(vI * dim) = (Vms.row(vI - mesh.componentNodeRange[compI]) - mesh.V.row(vI)).transpose();
                }
                break; // only for the first DCO
            }
        }
        break;
    }

    default:
        assert(0 && "invalid animScriptType");
        break;
    }

    targetPos.clear();
    for (int vI = 0; vI < mesh.V.rows(); ++vI) {
        if (mesh.isFixedVert[vI] || !searchDir.template segment<dim>(vI * dim).isZero()) {
            targetPos[vI].first = mesh.V.row(vI) + searchDir.template segment<dim>(vI * dim).transpose();
            targetPos[vI].second.setZero();
        }
    }
    dist2Tol = searchDir.squaredNorm() * 1.0e-6;

    stepSize = 1.0;
    for (const auto& energyTermI : energyTerms) {
        energyTermI->filterStepSize(mesh, searchDir, stepSize);
    }

    SpatialHash<dim> sh_CCD;
    if (isSelfCollision && !searchDir.isZero()) { // assume MBC does not collide with any other CO
#ifdef USE_SH_LFSS
        timer_temp3.start(11);
        sh_CCD.build(mesh, searchDir, stepSize, mesh.avgEdgeLen / 3.0);
        timer_temp3.stop();
#endif
        if (method == ExactCCD::Method::NONE) {
            std::vector<std::pair<int, int>> newCandidates;
            SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD(mesh, sh_CCD, searchDir, 0.5, newCandidates, stepSize);
        }
        else {
            SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD_exact(mesh, sh_CCD, searchDir, method, stepSize);
        }
    }

#ifdef OUTPUT_CCD_FAIL
    bool output = false;
    if (isSelfCollision && !searchDir.isZero() && method != ExactCCD::Method::NONE) {
        mesh.saveSurfaceMesh(outputFolderPath + "before" + std::to_string(numOfCCDFail) + ".obj");
        output = true;
    }
#endif

    // Check if the mesh is self intersecting before moving it
    bool isInitSelfIntersecting = false;
    if (!forceIntersectionLineSearch) {
        sh_CCD.build(mesh, mesh.avgEdgeLen / 3.0);
        isInitSelfIntersecting = !SelfCollisionHandler<dim>::checkEdgeTriIntersectionIfAny(mesh, sh_CCD);
    }

    // Move the vertices
    Eigen::MatrixXd V0 = mesh.V;
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.V.rows(), 1, [&](int vI)
#else
    for (int vI = 0; vI < mesh.V.rows(); ++vI)
#endif
        {
            mesh.V.row(vI) = V0.row(vI) + stepSize * searchDir.segment<dim>(vI * dim).transpose();
        }
#ifdef USE_TBB
    );
#endif

    // Line-search to prevent elment inversions
    if (energyTerms[0]->getNeedElemInvSafeGuard() && !searchDir.isZero()) {
        while (!mesh.checkInversion(true)) {
            logFile << "element inversion detected during step anim, backtrack!" << std::endl;
            stepSize /= 2.0;
#ifdef USE_TBB
            tbb::parallel_for(0, (int)mesh.V.rows(), 1, [&](int vI)
#else
            for (int vI = 0; vI < mesh.V.rows(); ++vI)
#endif
                {
                    mesh.V.row(vI) = V0.row(vI) + stepSize * searchDir.segment<dim>(vI * dim).transpose();
                }
#ifdef USE_TBB
            );
#endif
        }
    }

#ifndef NO_CCD_FAILSAFE
    // Line-search to prevent intersections (skip if there are already intersections)
    if ((forceIntersectionLineSearch || !isInitSelfIntersecting) && isSelfCollision && !searchDir.isZero()) {
        sh_CCD.build(mesh, mesh.avgEdgeLen / 3.0);
        while (!SelfCollisionHandler<dim>::checkEdgeTriIntersectionIfAny(mesh, sh_CCD)) {
#ifdef OUTPUT_CCD_FAIL
            if (output) {
                mesh.saveSurfaceMesh(outputFolderPath + "after" + std::to_string(numOfCCDFail) + ".obj");
                ++numOfCCDFail;
                output = false;
            }
#endif
            logFile << "intersection detected during step anim, backtrack!" << std::endl;
            stepSize /= 2.0;
#ifdef USE_TBB
            tbb::parallel_for(0, (int)mesh.V.rows(), 1, [&](int vI)
#else
            for (int vI = 0; vI < mesh.V.rows(); ++vI)
#endif
                {
                    mesh.V.row(vI) = V0.row(vI) + stepSize * searchDir.segment<dim>(vI * dim).transpose();
                }
#ifdef USE_TBB
            );
#endif
            sh_CCD.build(mesh, mesh.avgEdgeLen / 3.0);
        }
    }
#endif // NO_CCD_FAILSAFE

#ifdef CHECK_RATIONAL_CCD_GLOBAL
    Eigen::MatrixXd V_bk = mesh.V;

    mesh.V = V0;
    sh_CCD.build(mesh, searchDir, stepSize, mesh.avgEdgeLen / 3.0);
    double stepSize_rational = stepSize;
    SelfCollisionHandler<dim>::largestFeasibleStepSize_CCD_exact(mesh, sh_CCD, searchDir, ExactCCD::Method::RATIONAL_ROOT_PARITY, stepSize_rational);
    if (stepSize_rational != stepSize) {
        std::cout << "rational CCD detects interpenetration but inexact didn't" << std::endl;
    }

    mesh.V = V_bk;
#endif

    spdlog::debug("moved Dirichlet BC by {:g}", stepSize);
    spdlog::debug("moved CO by {:g}", getCOCompletedStepSize());

    return returnFlag;
}

template <int dim>
double AnimScripter<dim>::getCompletedStepSize(void) const
{
    return stepSize;
}

template <int dim>
double AnimScripter<dim>::computeCompletedStepSize(const Mesh<dim>& mesh)
{
    if (dist2Tol == 0.0) {
        stepSize = 1.0;
        return 1.0;
    }

    double sqNorm = 0.0;
    for (const auto& tpI : targetPos) {
        int vI = tpI.first;
        sqNorm += (mesh.V.row(vI) - tpI.second.first).squaredNorm();
    }
    stepSize = 1.0 - std::sqrt(sqNorm / (dist2Tol * 1.0e6));
    return stepSize;
}

template <int dim>
void AnimScripter<dim>::augmentMDBCEnergy(const Mesh<dim>& data, double& E, double rho) const
{
    for (const auto& tpI : targetPos) {
        int vI = tpI.first;
        E -= std::sqrt(data.massMatrix.coeff(vI, vI)) * tpI.second.second.cwiseProduct(data.V.row(vI) - tpI.second.first).sum();
        E += rho / 2.0 * data.massMatrix.coeff(vI, vI) * (data.V.row(vI) - tpI.second.first).squaredNorm();
    }
}

template <int dim>
void AnimScripter<dim>::augmentMDBCGradient(const Mesh<dim>& data, Eigen::VectorXd& g, double rho) const
{
    for (const auto& tpI : targetPos) {
        int vI = tpI.first;
        g.template segment<dim>(vI * dim) -= std::sqrt(data.massMatrix.coeff(vI, vI)) * tpI.second.second.transpose();
        g.template segment<dim>(vI * dim) += rho * data.massMatrix.coeff(vI, vI) * (data.V.row(vI) - tpI.second.first).transpose();
    }
}

template <int dim>
void AnimScripter<dim>::augmentMDBCHessian(const Mesh<dim>& data, LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSysSolver, double rho) const
{
    for (const auto& tpI : targetPos) {
        int vI = tpI.first;
        int startInd = vI * dim;
        double val = rho * data.massMatrix.coeff(vI, vI);
        linSysSolver->addCoeff(startInd, startInd, val);
        linSysSolver->addCoeff(startInd + 1, startInd + 1, val);
        if constexpr (dim == 3) {
            linSysSolver->addCoeff(startInd + 2, startInd + 2, val);
        }
    }
}

template <int dim>
void AnimScripter<dim>::updateLambda(const Mesh<dim>& data, double rho)
{
    for (auto& tpI : targetPos) {
        int vI = tpI.first;
        tpI.second.second -= rho * std::sqrt(data.massMatrix.coeff(vI, vI)) * (data.V.row(vI) - tpI.second.first);
    }
}

template <int dim>
void AnimScripter<dim>::resetLambda(void)
{
    for (auto& tpI : targetPos) {
        tpI.second.second.setZero();
    }
}

template <int dim>
double AnimScripter<dim>::getCOCompletedStepSize(void) const
{
    double maxStepSizeLeft = 0.0;
    for (const auto& i : ACOStepSizeLeft) {
        if (maxStepSizeLeft < i) {
            maxStepSizeLeft = i;
        }
    }
    for (const auto& i : MCOStepSizeLeft) {
        if (maxStepSizeLeft < i) {
            maxStepSizeLeft = i;
        }
    }
    return 1.0 - maxStepSizeLeft;
}

template <int dim>
bool AnimScripter<dim>::isNBCActive(void) const
{
    return curTime > NBCTimeRange[0] && curTime <= NBCTimeRange[1];
}

template <int dim>
void AnimScripter<dim>::setAnimScriptType(AnimScriptType p_animScriptType)
{
    animScriptType = p_animScriptType;
}

template <int dim>
void AnimScripter<dim>::setMeshSeqFolderPath(const std::string& p_meshSeqFolderPath)
{
    meshSeqFolderPath = p_meshSeqFolderPath;
}

template <int dim>
AnimScriptType AnimScripter<dim>::getAnimScriptTypeByStr(const std::string& str)
{
    for (int i = 0; i < animScriptTypeStrs.size(); i++) {
        if (str == animScriptTypeStrs[i]) {
            return AnimScriptType(i);
        }
    }
    spdlog::error("Invalid animScriptType: {:s}", str);
    return AST_NULL;
}
template <int dim>
std::string AnimScripter<dim>::getStrByAnimScriptType(AnimScriptType animScriptType)
{
    assert(animScriptType < animScriptTypeStrs.size());
    return animScriptTypeStrs[animScriptType];
}

template class AnimScripter<DIM>;

} // namespace IPC
