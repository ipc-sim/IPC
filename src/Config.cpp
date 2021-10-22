//
//  Config.cpp
//  IPC
//
//  Created by Minchen Li on 7/12/18.
//

#include "Config.hpp"
#include "IglUtils.hpp"

#include "HalfSpace.hpp"
#include "MeshCO.hpp"

#include <fstream>
#include <sstream>
#include <ctime>

#include <ghc/fs_std.hpp> // filesystem
#include <spdlog/spdlog.h>

namespace IPC {

const std::vector<std::string> Config::energyTypeStrs = {
    "NH", "FCR"
};
const std::vector<std::string> Config::timeIntegrationTypeStrs = {
    "BE", "NM"
};
const std::vector<std::string> Config::constraintSolverTypeStrs = {
    "QP", "SQP", "interiorPoint"
};
const std::vector<std::string> Config::constraintTypeStrs = {
    "volume", "graphics", "nonsmoothNewmark", "gapFunction", "CMR", "Verschoor", "STIV"
};
const std::vector<std::string> Config::QPSolverTypeStrs = {
    "OSQP", "Gurobi"
};

Config::Config(void)
{
    tuning.resize(6);
    tuning[0] = 0; // automatic kappa
    tuning[1] = tuning[2] = 1e-3; // dHat
    tuning[3] = 1e-9; // dHatEps
    tuning[4] = tuning[5] = 1e-3; // eps_v
}

Config::~Config(void)
{
}

std::string getRootDirectory()
{
    return fs::path(__FILE__).parent_path().parent_path().string();
}

std::string resolvePath(
    const std::string& path,
    const std::string& input_file_path)
{
    if (path.size() == 0) {
        throw fmt::format("empty path detected ({})", path);
    }

    if (path.at(0) == '/') {
        return path;
    }

    // Make the path relative to the input
    fs::path resolved_path = fs::path(input_file_path).parent_path() / path;
    if (fs::exists(resolved_path)) {
        return resolved_path.string();
    }

    // Return a path relative to the root directory of IPC
    return (fs::path(getRootDirectory()) / path).string();
}

template <typename T>
bool getOptional(std::stringstream& sstream, T& val, const T& defaultVal)
{
    if (sstream.eof()) {
        val = defaultVal;
        return false;
    }
    long pos = sstream.tellg();
    sstream >> val;
    if (sstream.fail()) {
        val = defaultVal;
        sstream.clear();
        sstream.seekg(pos);
        return false;
    }
    return true;
}

int Config::loadFromFile(const std::string& p_filePath)
{
    filePath = p_filePath;
    std::ifstream file(filePath);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string token;
            ss >> token;
            if (token == "energy") {
                std::string type;
                ss >> type;
                energyType = getEnergyTypeByStr(type);
            }
            else if (token == "timeIntegration") {
                std::string type;
                ss >> type;
                if (type == "NM") {
                    ss >> beta >> gamma;
                }
                timeIntegrationType = getTimeIntegrationTypeByStr(type);
            }
            else if (token == "linearSolver" || token == "linSysSolver") {
                std::string type;
                ss >> type;
                linSysSolverType = getLinSysSolverTypeByStr(type);
            }
            else if (token == "size") {
                ss >> size;
            }
            else if (token == "time") {
                ss >> duration >> dt;
            }
            else if (token == "density") {
                ss >> rho;
                assert(rho > 0.0);
            }
            else if (token == "stiffness") {
                ss >> YM >> PR;
            }
            else if (token == "turnOffGravity") {
                withGravity = false;
            }
            else if (token == "dampingStiff") {
                ss >> dampingStiff;
                if (dampingStiff < 0.0) {
                    dampingStiff = 0.0;
                }
                spdlog::info("dampingStiff = {}", dampingStiff);
            }
            else if (token == "dampingRatio") {
                ss >> dampingRatio;
                if (dampingRatio < 0.0) {
                    dampingRatio = 0.0;
                }
                else if (dampingRatio > 1.0) {
                    dampingRatio = 1.0;
                }
            }
            else if (token == "script") {
                std::string type;
                ss >> type;
                animScriptType = AnimScripter<DIM>::getAnimScriptTypeByStr(type);
                if (animScriptType == AST_MESHSEQ_FROMFILE) {
                    ss >> meshSeqFolderPath;
                    meshSeqFolderPath = resolvePath(meshSeqFolderPath, p_filePath);
                }

                int params = 0;
                ss >> params;
                if (params > 0) {
                    scriptParams.resize(params);
                    for (int i = 0; i < params; ++i) {
                        ss >> scriptParams[i];
                    }
                }
            }
            else if (token == "DBCTimeRange") {
                ss >> DBCTimeRange[0] >> DBCTimeRange[1];
            }
            else if (token == "NBCTimeRange") {
                ss >> NBCTimeRange[0] >> NBCTimeRange[1];
            }
            else if (token == "shape") {
                std::string type;
                ss >> type;
                std::string path;
                ss >> path;
                path = resolvePath(path, p_filePath);
                inputShapePaths.push_back(path);
                inputShapeTranslates.push_back(Eigen::Vector3d::Zero());
                inputShapeRotates.push_back(Eigen::Matrix3d::Identity());
                inputShapeScales.push_back(Eigen::Vector3d::Ones());
                inputShapeMaterials.push_back(Eigen::Vector3d::Constant(NAN));
                inputShapeLVels.emplace_back(NAN, NAN, NAN);
                inputShapeAVels.emplace_back(NAN, NAN, NAN);
                inputShapeInitVels.push_back(
                    { Eigen::Vector3d(NAN, NAN, NAN),
                        Eigen::Vector3d(NAN, NAN, NAN) });
            }
            else if (token == "shapes") {
                std::string type;
                ss >> type;
                int shapeNum;
                ss >> shapeNum;
                for (int shapeI = 0; shapeI < shapeNum; ++shapeI) {
                    std::string line_shapes;
                    std::getline(file, line_shapes);
                    std::stringstream ss_shapes(line_shapes);

                    std::string path;
                    ss_shapes >> path;
                    if (path.empty() || path[0] == '#') {
                        shapeI--; // This decrease will be canceled out
                        continue;
                    }
                    path = resolvePath(path, p_filePath);
                    inputShapePaths.push_back(path);
                    double x, y, z;
                    ss_shapes >> x >> y >> z;
                    inputShapeTranslates.push_back(Eigen::Vector3d(x, y, z));
                    ss_shapes >> x >> y >> z;
                    Eigen::Matrix3d rotMat = (Eigen::AngleAxisd(x / 180.0 * M_PI, Eigen::Vector3d::UnitX())
                        * Eigen::AngleAxisd(y / 180.0 * M_PI, Eigen::Vector3d::UnitY())
                        * Eigen::AngleAxisd(z / 180.0 * M_PI, Eigen::Vector3d::UnitZ()))
                                                 .toRotationMatrix();
                    inputShapeRotates.push_back(rotMat);
                    ss_shapes >> x >> y >> z;
                    inputShapeScales.push_back(Eigen::Vector3d(x, y, z));

                    double density = NAN, E = NAN, nu = NAN;
                    double vx = NAN, vy = NAN, vz = NAN, avx = NAN, avy = NAN, avz = NAN;
                    double initLV[3] = { NAN, NAN, NAN }, initAV[3] = { NAN, NAN, NAN };
                    std::string extra;
                    while (!ss_shapes.eof()) {
                        ss_shapes >> extra;
                        if (extra == "material") {
                            ss_shapes >> density >> E >> nu;
                        }
                        else if (extra == "linearVelocity") {
                            ss_shapes >> vx >> vy >> vz;
                        }
                        else if (extra == "angularVelocity") {
                            ss_shapes >> avx >> avy >> avz;
                            avx *= M_PI / 180;
                            avy *= M_PI / 180;
                            avz *= M_PI / 180;
                        }
                        else if (extra == "DBC") {
                            double DBC[12];
                            for (int i = 0; i < 12; ++i) {
                                ss_shapes >> DBC[i];
                            }
                            // Optional start and stop time of DBC
                            std::array<double, 2> timeRange;
                            getOptional(ss_shapes, timeRange[0], 0.0);
                            getOptional(ss_shapes, timeRange[1], std::numeric_limits<double>::infinity());

                            inputShapeDBC.emplace_back(shapeI,
                                InputDBC(
                                    Eigen::Vector3d(DBC[0], DBC[1], DBC[2]), // bbox min
                                    Eigen::Vector3d(DBC[3], DBC[4], DBC[5]), // bbox max
                                    Eigen::Vector3d(DBC[6], DBC[7], DBC[8]), // linear velocity
                                    Eigen::Vector3d(DBC[9], DBC[10], DBC[11]) * M_PI / 180, // angular velocity
                                    timeRange));
                        }
                        else if (extra == "NBC") {
                            double NBC[9];
                            for (int i = 0; i < 9; ++i) {
                                ss_shapes >> NBC[i];
                            }
                            // Optional start and stop time of NBC
                            std::array<double, 2> timeRange;
                            getOptional(ss_shapes, timeRange[0], 0.0);
                            getOptional(ss_shapes, timeRange[1], std::numeric_limits<double>::infinity());

                            inputShapeNBC.emplace_back(shapeI,
                                InputNBC(
                                    Eigen::Vector3d(NBC[0], NBC[1], NBC[2]), // bbox min
                                    Eigen::Vector3d(NBC[3], NBC[4], NBC[5]), // bbox max
                                    Eigen::Vector3d(NBC[6], NBC[7], NBC[8]), // force
                                    timeRange));
                        }
                        else if (extra == "initVel") {
                            ss_shapes >> initLV[0] >> initLV[1] >> initLV[2] >> initAV[0] >> initAV[1] >> initAV[2];
                        }
                        else if (extra == "meshSeq") {
                            std::string meshSeqFolderPath;
                            ss_shapes >> meshSeqFolderPath;
                            meshSeqFolderPath = resolvePath(meshSeqFolderPath, p_filePath);
                            inputShapeMeshSeqFolderPath.emplace_back(shapeI, meshSeqFolderPath);
                        }
                        else if (extra == "\\") {
                            std::getline(file, line_shapes);
                            ss_shapes = std::stringstream(line_shapes);
                        }
                        else if (extra[0] == '#') {
                            while (!ss_shapes.eof()) {
                                ss_shapes >> extra;
                                if (extra == "\\") {
                                    std::getline(file, line_shapes);
                                    ss_shapes = std::stringstream(line_shapes);
                                    break;
                                }
                            }
                        }
                        else {
                            spdlog::error("Uknown keyword in shape line {:d}: {}", shapeI, extra);
                        }

                        if (!ss_shapes.eof() && ss_shapes.fail()) {
                            spdlog::error("Failed to parse {} section of shape line {:d} ({})", extra, shapeI, ss_shapes.str());
                            break;
                        }
                    }
                    inputShapeMaterials.push_back(Eigen::Vector3d(density, E, nu));
                    inputShapeLVels.push_back(Eigen::Vector3d(vx, vy, vz));
                    inputShapeAVels.push_back(Eigen::Vector3d(avx, avy, avz));
                    inputShapeInitVels.push_back({ Eigen::Vector3d(initLV[0], initLV[1], initLV[2]), Eigen::Vector3d(initAV[0], initAV[1], initAV[2]) * M_PI / 180 });
                }
            }
            else if (token == "shapeMatrix") {
                std::string type;
                ss >> type;
                int count[3];
                ss >> count[0] >> count[1] >> count[2];
                double posX, posY, posZ;
                ss >> posX >> posY >> posZ;
                if (posX != posX) {
                    posX = 0.0;
                }
                if (posY != posY) {
                    posY = 0.0;
                }
                if (posZ != posZ) {
                    posZ = 0.0;
                }

                std::string line_shapes;
                std::getline(file, line_shapes);
                std::stringstream ss_shapes(line_shapes);

                std::string path;
                ss_shapes >> path;
                path = resolvePath(path, p_filePath);

                double x, y, z;
                ss_shapes >> x >> y >> z;
                double translationStep[3] = { x, y, z };

                ss_shapes >> x >> y >> z;
                Eigen::Matrix3d rotMat = (Eigen::AngleAxisd(x / 180.0 * M_PI, Eigen::Vector3d::UnitX())
                    * Eigen::AngleAxisd(y / 180.0 * M_PI, Eigen::Vector3d::UnitY())
                    * Eigen::AngleAxisd(z / 180.0 * M_PI, Eigen::Vector3d::UnitZ()))
                                             .toRotationMatrix();

                ss_shapes >> x >> y >> z;
                Eigen::Vector3d scale(x, y, z);

                Eigen::Vector3d material;
                material.setConstant(NAN);
                std::string extra;
                ss_shapes >> extra;
                if (extra == "material") {
                    ss_shapes >> material[0] >> material[1] >> material[2];
                }

                for (int xi = 0; xi < count[0]; ++xi) {
                    for (int yi = 0; yi < count[1]; ++yi) {
                        for (int zi = 0; zi < count[2]; ++zi) {
                            inputShapePaths.emplace_back(path);
                            inputShapeTranslates.emplace_back(Eigen::Vector3d(posX + translationStep[0] * xi,
                                posY + translationStep[1] * yi, posZ + translationStep[2] * zi));
                            inputShapeRotates.emplace_back(rotMat);
                            inputShapeScales.emplace_back(scale);
                            inputShapeMaterials.emplace_back(material);
                        }
                    }
                }
            }
            else if (token == "tol") {
                int amt;
                ss >> amt;
                assert(amt >= 0);
                tol.resize(amt);
                for (auto& tolI : tol) {
                    file >> tolI;
                }
            }
            else if (token == "noActiveSetConvergence") {
                useActiveSetConvergence = false;
            }
            else if (token == "useActiveSetConvergence") {
                useActiveSetConvergence = true;
            }
            else if (token == "warmStart") {
                ss >> warmStart;
                assert(warmStart >= 0);
            }
            else if (token == "view") {
                std::string type;
                ss >> type;
                if (type == "orthographic") {
                    orthographic = true;
                }
                else if (type == "perspective") {
                    orthographic = false;
                }
                else {
                    orthographic = false;
                    spdlog::info("use default perspective view");
                }
            }
            else if (token == "zoom") {
                ss >> zoom;
                assert(zoom > 0.0);
            }
            else if (token == "cameraTracking") {
                cameraTracking = true;
            }
            else if (token == "playBackSpeed") {
                ss >> playBackSpeed;
                if (playBackSpeed <= 0.0) {
                    playBackSpeed = 1.0;
                }
            }

            else if (token == "ground") {
                double groundFriction, groundY;
                ss >> groundFriction >> groundY;
                assert(groundFriction >= 0.0);
                collisionObjects.push_back(std::make_shared<HalfSpace<DIM>>(groundY, groundFriction));
            }
            else if (token == "halfSpace") {
                Eigen::Matrix<double, DIM, 1> origin, normal;
                ss >> origin[0] >> origin[1];
                if constexpr (DIM == 3) {
                    ss >> origin[2];
                }
                ss >> normal[0] >> normal[1];
                if constexpr (DIM == 3) {
                    ss >> normal[2];
                }
                normal.normalize();
                double stiffness, friction;
                ss >> stiffness >> friction;
                assert(friction >= 0.0);
                collisionObjects.push_back(std::make_shared<HalfSpace<DIM>>(
                    origin, normal, Eigen::Matrix<double, DIM, 1>::Zero(), friction));
            }
            else if (token == "meshCO") {
                std::string meshCOFilePath;
                ss >> meshCOFilePath;
                meshCOFilePath = resolvePath(meshCOFilePath, p_filePath);

                Eigen::Matrix<double, DIM, 1> origin;
                ss >> origin[0] >> origin[1];
                if constexpr (DIM == 3) {
                    ss >> origin[2];
                }

                double scale, stiffness, friction;
                ss >> scale >> stiffness >> friction;

                Eigen::Matrix3d rotMat = Eigen::Matrix3d::Identity();
                std::string extraSetting;
                ss >> extraSetting;
                if (extraSetting == "rotate") {
                    double x, y, z;
                    ss >> x >> y >> z;
                    rotMat = (Eigen::AngleAxisd(x / 180.0 * M_PI, Eigen::Vector3d::UnitX())
                        * Eigen::AngleAxisd(y / 180.0 * M_PI, Eigen::Vector3d::UnitY())
                        * Eigen::AngleAxisd(z / 180.0 * M_PI, Eigen::Vector3d::UnitZ()))
                                 .toRotationMatrix();
                }

                meshCollisionObjects.emplace_back(std::make_shared<MeshCO<DIM>>(meshCOFilePath.c_str(), origin, rotMat, scale, friction));
            }
            else if (token == "selfCollisionOn") {
                isSelfCollision = true;
            }
            else if (token == "selfCollisionOff") {
                isSelfCollision = false;
            }
            else if (token == "selfFric") {
                ss >> selfFric;
                if (selfFric < 0.0) {
                    selfFric = 0.0;
                }
            }

            else if (token == "constraintSolver") {
                std::string type;
                ss >> type;
                if (type == "IP") {
                    type = "interiorPoint";
                }
                constraintSolverType = getConstraintSolverTypeByStr(type);
                spdlog::info("Using collision constraint solver: {:s}",
                    getStrByConstraintSolverType(constraintSolverType));
            }
            else if (token == "constraintType") {
                std::string type;
                ss >> type;
                constraintType = getConstraintTypeByStr(type);
            }
            else if (token == "QPSolver") {
                std::string type;
                ss >> type;
                qpSolverType = getQPSolverTypeByStr(type);
            }
            else if (token == "appendStr") {
                ss >> appendStr;
            }

            else if (token == "restart") {
                restart = true;
                ss >> statusPath;
                statusPath = resolvePath(statusPath, p_filePath);
            }

            else if (token == "disableCout") {
                disableCout = true;
            }

            else if (token == "rotateModel") {
                assert(DIM == 3);
                ss >> rotAxis[0] >> rotAxis[1] >> rotAxis[2] >> rotDeg;
            }

            else if (token == "handleRatio") {
                ss >> handleRatio;
                assert((handleRatio > 0) && (handleRatio < 0.5));
            }

            else if (token == "tuning") {
                int amt;
                ss >> amt;
                assert(amt >= 0);
                tuning.resize(amt);
                for (auto& tuneI : tuning) {
                    file >> tuneI;
                }
            }
            else if (token == "dHat") {
                ss >> tuning[1];
                tuning[2] = tuning[1];
            }
            else if (token == "epsv") {
                ss >> tuning[4];
                tuning[5] = tuning[4];
            }
            else if (token == "fricIterAmt") {
                ss >> fricIterAmt;
            }
            else if (token == "useAbsParameters") {
                useAbsParameters = true;
            }
            else if (token == "kappaMinMultiplier" || token == "minBarrierStiffnessScale") {
                ss >> kappaMinMultiplier;
            }

            else if (token == "constraintOffset") {
                ss >> constraintOffset;
                spdlog::info("Using collision constraint offset: {:g}", constraintOffset);
            }
            else if (token == "CCDMethod" || token == "ccdMethod") {
                std::string type;
                ss >> type;
                this->ccdMethod = this->getCCDMethodTypeByStr(type);
            }
            else if (token == "CCDTolerance" || token == "ccdTolerance") {
                ss >> this->ccdTolerance;
            }
            else if (token == "section") {
                std::string section;
                bool sectionMatches = false;
                while (!ss.eof()) {
                    ss >> section;
                    if (section == "end") {
                        sectionMatches = true;
                        break;
                    }

                    if (section == "IP") {
                        section = "interiorPoint";
                    }

                    if (section == getStrByConstraintSolverType(constraintSolverType)) {
                        sectionMatches = true;
                        break;
                    }
                }

                if (!sectionMatches) {
                    spdlog::warn("Skipping section: {:s}", section);
                    while (std::getline(file, line)) {
                        ss = std::stringstream(line);
                        ss >> token;
                        if (token == "section") {
                            ss >> token;
                            if (token == "end") {
                                break;
                            }
                        }
                    }
                }
            }
            else if (token != "" && token[0] != '#') {
                spdlog::error("Unrecognized configuration option: {}", token);
            }
        }

        file.close();

        isConstrained = isSelfCollision || !(collisionObjects.empty() && meshCollisionObjects.empty()) || constraintSolverType == CST_IP;
        if (dampingRatio > 0) {
            dampingStiff = dampingRatio * std::pow(dt, 3) * 3 / 4;
        }
        return 0;
    }
    else {
        spdlog::error("Unable to open input file: {:s}", p_filePath);
        return -1;
    }
}

void Config::backUpConfig(const std::string& p_filePath)
{
    std::ifstream src(filePath, std::ios::binary);
    std::ofstream dst(p_filePath, std::ios::binary);
    dst << src.rdbuf();
}

void Config::appendInfoStr(std::string& inputStr) const
{
    std::string shapeName;
    int shapeI = 0;
    for (const auto& inputShapePathI : inputShapePaths) {
        shapeName += fs::path(inputShapePathI).stem().string() + "-";
        if (++shapeI >= 3) {
            if (shapeI > 3) {
                shapeName += std::to_string(inputShapePaths.size() - 3) + "others_";
            }
            break;
        }
    }
    shapeName.pop_back();

    inputStr += (shapeName + "_" + AnimScripter<DIM>::getStrByAnimScriptType(animScriptType) + "_" + getStrByEnergyType(energyType) + "_" + getStrByTimeIntegrationType(timeIntegrationType) + "_" + (isConstrained ? getStrByConstraintSolverType(constraintSolverType) : "") + "_");

    time_t rawTime = std::time(NULL);
    char buf[BUFSIZ];
    std::strftime(buf, sizeof(buf), "%Y%m%d%H%M%S", std::localtime(&rawTime));
    inputStr += buf;

    if (appendStr.length()) {
        inputStr += "_" + appendStr;
    }
}

EnergyType Config::getEnergyTypeByStr(const std::string& str)
{
    for (int i = 0; i < energyTypeStrs.size(); i++) {
        if (str == energyTypeStrs[i]) {
            return EnergyType(i);
        }
    }
    spdlog::warn("Unknown energy type: {}; using default energy type: NH", str);
    return ET_NH;
}
std::string Config::getStrByEnergyType(EnergyType energyType)
{
    assert(energyType < energyTypeStrs.size());
    return energyTypeStrs[energyType];
}
TimeIntegrationType Config::getTimeIntegrationTypeByStr(const std::string& str)
{
    for (int i = 0; i < timeIntegrationTypeStrs.size(); i++) {
        if (str == timeIntegrationTypeStrs[i]) {
            return TimeIntegrationType(i);
        }
    }
    spdlog::warn("Unknown time integration: {}; using default time integration type: BE", str);
    return TIT_BE;
}
std::string Config::getStrByTimeIntegrationType(TimeIntegrationType timeIntegrationType)
{
    assert(timeIntegrationType < timeIntegrationTypeStrs.size());
    return timeIntegrationTypeStrs[timeIntegrationType];
}
LinSysSolverType Config::getLinSysSolverTypeByStr(const std::string& solver)
{
#ifdef IPC_WITH_CHOLMOD
    if (solver == "cholmod" || solver == "CHOLMOD") {
        return LinSysSolverType::CHOLMOD;
    }
#endif
#ifdef IPC_WITH_AMGCL
    if (solver == "amgcl" || solver == "AMGCL") {
        return LinSysSolverType::AMGCL;
    }
#endif
    if (solver == "eigen" || solver == "EIGEN" || solver == "Eigen") {
        return LinSysSolverType::EIGEN;
    }
    spdlog::warn("Uknown linear system solver: {}; using default linear system solver: Eigen", solver);
    return LinSysSolverType::EIGEN;
}
ConstraintSolverType Config::getConstraintSolverTypeByStr(const std::string& str)
{
    for (int i = 0; i < constraintSolverTypeStrs.size(); i++) {
        if (str == constraintSolverTypeStrs[i]) {
            return ConstraintSolverType(i);
        }
    }
    spdlog::warn("Uknown CCD method: {}; using default constraint solver type: IP", str);
    return CST_IP;
}
std::string Config::getStrByConstraintSolverType(ConstraintSolverType constraintSolverType)
{
    assert(constraintSolverType < constraintSolverTypeStrs.size());
    return constraintSolverTypeStrs[constraintSolverType];
}
ccd::CCDMethod Config::getCCDMethodTypeByStr(const std::string& str)
{
    for (int i = 0; i < ccd::CCDMethod::NUM_CCD_METHODS; i++) {
        if (str == ccd::method_names[i]) {
            return ccd::CCDMethod(i);
        }
    }
    spdlog::warn("Uknown CCD method: {}; using default CCD method: {}",
        str, ccd::method_names[ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER]);
    return ccd::CCDMethod::FLOATING_POINT_ROOT_FINDER;
}
CollisionConstraintType Config::getConstraintTypeByStr(const std::string& str)
{
    for (int i = 0; i < constraintTypeStrs.size(); i++) {
        if (str == constraintTypeStrs[i]) {
            spdlog::info("Using collision constraint type: {:s}", constraintTypeStrs[i]);
            return CollisionConstraintType(i);
        }
    }
    spdlog::warn("Uknown collision constraint type: {:s}; using default collision constraint type: volume", str);
    return CollisionConstraintType::VOLUME;
}
QPSolverType Config::getQPSolverTypeByStr(const std::string& str)
{
    for (int i = 0; i < QPSolverTypeStrs.size(); i++) {
        if (str == QPSolverTypeStrs[i]) {
            spdlog::info("Using QP solver: {:s}", QPSolverTypeStrs[i]);
            return QPSolverType(i);
        }
    }
    spdlog::warn("Uknown QP solver: {:s}; using default QP solver: {:s}", str, QPSolverTypeStrs[0]);
    return QPSolverType(0);
}

} // namespace IPC
