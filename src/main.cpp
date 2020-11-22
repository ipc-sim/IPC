#include "Types.hpp"
#include "IglUtils.hpp"
#include "Config.hpp"
#include "Optimizer.hpp"
#include "NeoHookeanEnergy.hpp"
#include "FixedCoRotEnergy.hpp"
#include "GIF.hpp"
#include "Timer.hpp"
#include "getRSS.hpp"

#include <igl/readOBJ.h>
#ifdef USE_OPENGL
#include <igl/opengl/glfw/Viewer.h>
#include <igl/png/writePNG.h>
#endif

#include <sys/stat.h> // for mkdir

#include <fstream>
#include <string>
#include <ctime>

#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>

// optimization/simulation
IPC::Config config;
std::vector<const IPC::Mesh<DIM>*> triSoup;
int vertAmt_input;
IPC::Optimizer<DIM>* optimizer;
std::vector<IPC::Energy<DIM>*> energyTerms;
std::vector<double> energyParams;
bool savePNG = false;
bool offlineMode = false;
bool autoSwitch = false;
bool contactHandling = false;
bool optimization_on = false;
bool show_lines = true;
int iterNum = 0;
int converged = 0;
bool outerLoopFinished = false;
int numOfCCDFail = 0;

std::ofstream logFile;
std::string outputFolderPath = "output/";

// visualization
#ifdef USE_OPENGL
igl::opengl::glfw::Viewer viewer;
#endif
const int channel_initial = 0;
const int channel_result = 1;
int viewChannel = channel_result;
Eigen::MatrixXi SF;
std::vector<int> sTri2Tet;
bool viewUV = true; // view UV or 3D model
bool showSeam = false;
Eigen::MatrixXd seamColor, floorColor;
bool showBoundary = false;
int showDistortion = 2; // 0: don't show; 1: energy value; 2: other scalar field;
int showDistortion_init = showDistortion;
Eigen::MatrixXd faceColors_default;
bool showTexture = true; // show checkerboard
bool isLighting = true;
bool showFixedVerts = true;

std::vector<bool> isSurfNode;
std::vector<int> tetIndToSurf;
std::vector<int> surfIndToTet;
Eigen::MatrixXd V_surf;
Eigen::MatrixXi F_surf;

std::vector<bool> isCENode;
std::vector<int> tetIndToCE;
std::vector<int> CEIndToTet;
Eigen::MatrixXd V_CE;
Eigen::MatrixXi F_CE;

std::vector<int> compVAccSize;
std::vector<int> compFAccSize;

bool saveInfo_postDraw = false;
std::string infoName = "";

GifWriter GIFWriter;
uint32_t GIFDelay = 10; //*10ms
int GIFStep = 1;
double GIFScale = 0.6;

// timer
double secPast = 0.0;
time_t lastStart_world;
Timer timer, timer_step, timer_temp3, timer_mt;

void saveInfo(bool writePNG = true, bool writeGIF = true, int writeMesh = 1, double save_dt = 1e-2);
void saveScreenshot(const std::string& filePath, double scale = 1.0, bool writeGIF = false, bool writePNG = true);
void saveInfoForPresent(const std::string fileName = "info.txt", double save_dt = 1e-2);

void proceedOptimization(int proceedNum = 1)
{
    for (int proceedI = 0; (proceedI < proceedNum) && (!converged); proceedI++) {
        infoName = std::to_string(iterNum);
        saveInfo(false, false, 3);
        if (!offlineMode) {
            // PNG and gif output only works under online rendering mode
            saveInfo(false, true, false); // save gif
            if (savePNG) {
                saveScreenshot(outputFolderPath + infoName + ".png",
                    1.0, false, true);
            }
        }
        showDistortion = showDistortion_init;

        if (iterNum < config.tol.size()) {
            optimizer->setRelGL2Tol(config.tol[iterNum]);
        }
        else {
            if (config.tol.size()) {
                optimizer->setRelGL2Tol(config.tol.back());
            }
            else {
                optimizer->setRelGL2Tol();
            }
        }

        converged = optimizer->solve(1);
        if (converged == 2) {
            showDistortion = 0;
            converged = 0;
            logFile << "!!! maxIter reached for timeStep" << iterNum << std::endl;
        }
        iterNum = optimizer->getIterNum();

        saveInfoForPresent("info" + std::to_string(iterNum) + ".txt");

#ifdef FIRST_TIME_STEP
        converged = true;
#endif
    }
}

void updateViewerData_distortion(void)
{
#ifdef USE_OPENGL
    Eigen::MatrixXd color_distortionVis;

    switch (showDistortion) {
    case 2: { // show other triangle-based scalar fields
        Eigen::VectorXd l2StretchPerElem;
        Eigen::VectorXd faceWeight;
#if (DIM == 2)
        optimizer->getFaceFieldForVis(faceWeight);
#else
        Eigen::VectorXd faceWeight_tet;
        optimizer->getFaceFieldForVis(faceWeight_tet);
        faceWeight.conservativeResize(SF.rows());
        for (int sfI = 0; sfI < SF.rows(); sfI++) {
            if (sTri2Tet[sfI] < 0) {
                faceWeight[sfI] = 0.0; // surface mesh
            }
            else {
                faceWeight[sfI] = faceWeight_tet[sTri2Tet[sfI]];
            }
        }
        faceWeight *= -1.0;
#endif
        igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, faceWeight, true, color_distortionVis);
        color_distortionVis.array() += 0.2;
        color_distortionVis = color_distortionVis.array().min(1.0);
        break;
    }

    case 0: {
#if (DIM == 2)
        color_distortionVis = Eigen::MatrixXd::Ones(triSoup[viewChannel]->F.rows(), 3);
#else
        color_distortionVis = Eigen::MatrixXd::Ones(SF.rows(), 3);
#endif
        color_distortionVis.col(2).setZero();
        break;
    }

    case -1: {
        assert(faceColors_default.rows() == ((DIM == 3) ? SF.rows() : triSoup[viewChannel]->F.rows()));
        color_distortionVis = faceColors_default;
        break;
    }

    default:
        assert(0 && "unknown distortion visualization option!");
        break;
    }

    if (showSeam) {
        color_distortionVis.conservativeResize(color_distortionVis.rows() + seamColor.rows(), 3);
        color_distortionVis.bottomRows(seamColor.rows()) = seamColor;
    }

    if (floorColor.rows() > 0) {
        color_distortionVis.conservativeResize(color_distortionVis.rows() + floorColor.rows(), 3);
        color_distortionVis.bottomRows(floorColor.rows()) = floorColor;
    }

    viewer.data().set_colors(color_distortionVis);
#endif
}

void updateViewerData(void)
{
#ifdef USE_OPENGL
    Eigen::MatrixXd UV_vis = triSoup[viewChannel]->V;
    Eigen::MatrixXi F_vis = ((DIM == 2) ? triSoup[viewChannel]->F : SF);
    if (viewUV) {
        if constexpr (DIM == 2) {
            UV_vis.conservativeResize(UV_vis.rows(), 3);
            UV_vis.rightCols(1) = Eigen::VectorXd::Zero(UV_vis.rows());
        }
        viewer.core().align_camera_center((config.cameraTracking ? triSoup[viewChannel]->V : triSoup[viewChannel]->V_rest), F_vis);
        floorColor.conservativeResize(0, 3);
        for (const auto& coI : config.collisionObjects) {
            coI->draw(UV_vis, F_vis, floorColor);
        }
        for (const auto& coI : config.meshCollisionObjects) {
            coI->draw(UV_vis, F_vis, floorColor);
        }

        if ((UV_vis.rows() != viewer.data().V.rows()) || (F_vis.rows() != viewer.data().F.rows())) {
            viewer.data().clear();
        }
        viewer.data().set_mesh(UV_vis, F_vis);

        viewer.data().show_texture = false;
#if (DIM == 2)
        viewer.core().lighting_factor = 0.0;
#else
        if (isLighting) {
            viewer.core().lighting_factor = 0.6;
        }
        else {
            viewer.core().lighting_factor = 0.0;
        }
#endif

        viewer.data().set_points(Eigen::MatrixXd::Zero(0, 3), Eigen::RowVector3d(0.0, 0.0, 0.0));
        if (showFixedVerts) {
            for (const auto& fixedVI : triSoup[viewChannel]->fixedVert) {
                viewer.data().add_points(UV_vis.row(fixedVI), Eigen::RowVector3d(0.0, 0.0, 0.0));
            }
        }

        // draw codimensional segment collision objects
        viewer.data().clear_edges();
        for (int ceI = 0; ceI < triSoup[viewChannel]->CE.rows(); ++ceI) {
            viewer.data().add_edges(triSoup[viewChannel]->V.row(triSoup[viewChannel]->CE(ceI, 0)),
                triSoup[viewChannel]->V.row(triSoup[viewChannel]->CE(ceI, 1)),
                Eigen::RowVector3d(0.0, 0.0, 0.0));
        }
    }
    else {
        Eigen::MatrixXd V_vis = triSoup[viewChannel]->V_rest;
        viewer.core().align_camera_center(V_vis, F_vis);

        if ((V_vis.rows() != viewer.data().V.rows()) || (UV_vis.rows() != viewer.data().V_uv.rows()) || (F_vis.rows() != viewer.data().F.rows())) {
            viewer.data().clear();
        }
        viewer.data().set_mesh(V_vis, F_vis);

        if (showTexture) {
            viewer.data().set_uv(UV_vis);
            viewer.data().show_texture = true;
        }
        else {
            viewer.data().show_texture = false;
        }

        if (isLighting) {
            viewer.core().lighting_factor = 1.0;
        }
        else {
            viewer.core().lighting_factor = 0.0;
        }

        viewer.data().set_points(Eigen::MatrixXd::Zero(0, 3), Eigen::RowVector3d(0.0, 0.0, 0.0));
        if (showFixedVerts) {
            for (const auto& fixedVI : triSoup[viewChannel]->fixedVert) {
                viewer.data().add_points(V_vis.row(fixedVI), Eigen::RowVector3d(0.0, 0.0, 0.0));
            }
        }
    }
    updateViewerData_distortion();

    viewer.data().compute_normals();
#endif
}

void saveScreenshot(const std::string& filePath, double scale, bool writeGIF, bool writePNG)
{
#ifdef USE_OPENGL
    if (offlineMode) {
        return;
    }

    if (writeGIF) {
        scale = GIFScale;
    }
    viewer.data().point_size = 5 * scale;

    int width = static_cast<int>(scale * (viewer.core().viewport[2] - viewer.core().viewport[0]));
    int height = static_cast<int>(scale * (viewer.core().viewport[3] - viewer.core().viewport[1]));

    // Allocate temporary buffers for image
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R(width, height);
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> G(width, height);
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> B(width, height);
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A(width, height);

    // Draw the scene in the buffers
    viewer.core().draw_buffer(viewer.data(), false, R, G, B, A);

    if (writePNG) {
        // Save it to a PNG
        igl::png::writePNG(R, G, B, A, filePath);
    }

    if (writeGIF && (iterNum % GIFStep == 0)) {
        std::vector<uint8_t> img(width * height * 4);
        for (int rowI = 0; rowI < width; rowI++) {
            for (int colI = 0; colI < height; colI++) {
                int indStart = (rowI + (height - 1 - colI) * width) * 4;
                img[indStart] = R(rowI, colI);
                img[indStart + 1] = G(rowI, colI);
                img[indStart + 2] = B(rowI, colI);
                img[indStart + 3] = A(rowI, colI);
            }
        }
        GifWriteFrame(&GIFWriter, img.data(), width, height, GIFDelay);
    }

    viewer.data().point_size = 5;
#endif
}

void saveInfo(bool writePNG, bool writeGIF, int writeMesh, double save_dt)
{
    if (config.playBackSpeed < 1.0) {
        save_dt = std::max(config.dt, std::min(save_dt, config.playBackSpeed * GIFDelay / 100.0));
    }

    saveScreenshot(outputFolderPath + infoName + ".png", 0.5, writeGIF, writePNG);

    // Limit the number of files to 1 / save_dt per simulation second
    static int previousSaveIter = std::numeric_limits<int>::min() / 2;
    int itersBetweenSaves = round(save_dt / config.dt);
    if (iterNum != previousSaveIter && iterNum - previousSaveIter < itersBetweenSaves && !converged) {
        return;
    }
    previousSaveIter = iterNum;

    switch (writeMesh) {
    case 1:
        // triSoup[channel_result]->save(outputFolderPath + infoName + "_triSoup.obj");
        triSoup[channel_result]->saveAsMesh(outputFolderPath + infoName + "_mesh" + ((DIM == 2) ? ".obj" : ".msh"), false, SF);
        break;

    case 3: // save status
        optimizer->saveStatus();
        break;

    default:
        break;
    }
}

void saveInfoForPresent(const std::string fileName, double save_dt)
{
    if (config.playBackSpeed < 1.0) {
        save_dt = std::max(config.dt, std::min(save_dt, config.playBackSpeed * GIFDelay / 100.0));
    }

    // Limit the number of files to 1 / save_dt per simulation second
    static int previousSaveIter = std::numeric_limits<int>::min() / 2;
    int itersBetweenSaves = round(save_dt / config.dt);
    if (iterNum != previousSaveIter && iterNum - previousSaveIter < itersBetweenSaves && !converged) {
        return;
    }
    previousSaveIter = iterNum;

    std::ofstream file;
    file.open(outputFolderPath + fileName);
    assert(file.is_open());

    file << vertAmt_input << " " << triSoup[channel_initial]->F.rows() << std::endl;

    file << iterNum << " " << optimizer->getInnerIterAmt() << " " << double(optimizer->getInnerIterAmt()) / iterNum << std::endl;

    file << "current " << getCurrentRSS() / 1024.0 << " KB, peak " << getPeakRSS() / 1024.0 << " KB" << std::endl;

    timer.print(file);
    timer_step.print(file);
    timer_temp3.print(file);
    timer_mt.print(file);

    optimizer->outputCollStats(file);

    file.close();
}

void saveStats(void)
{
    FILE* out = fopen("output/resultsStats.txt", "a+");
    assert(out);

    fprintf(out, "%s %lf %lf\n", outputFolderPath.c_str(),
        double(optimizer->getInnerIterAmt()) / iterNum,
        timer.timing(1) / iterNum);

    fclose(out);
}

void toggleOptimization(void)
{
    optimization_on = !optimization_on;
    if (optimization_on) {
        if (converged) {
            optimization_on = false;
            spdlog::info("optimization converged.");
        }
        else {
#ifdef USE_OPENGL
            if (iterNum == 0) {
                GifBegin(&GIFWriter, (outputFolderPath + "anim.gif").c_str(),
                    GIFScale * (viewer.core().viewport[2] - viewer.core().viewport[0]),
                    GIFScale * (viewer.core().viewport[3] - viewer.core().viewport[1]), GIFDelay);

                saveScreenshot(outputFolderPath + "0.png", 0.5, true);
            }
            spdlog::info("start/resume optimization, press again to pause.");
            viewer.core().is_animating = true;
#endif
            time(&lastStart_world);
        }
    }
    else {
        spdlog::info("pause optimization, press again to resume.");
#ifdef USE_OPENGL
        viewer.core().is_animating = false;
#endif
        spdlog::info("World Time:\nTime past: {:g} s.", secPast);
        secPast += difftime(time(NULL), lastStart_world);
    }
}

#ifdef USE_OPENGL
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    if ((key >= '0') && (key <= '9')) {
        int changeToChannel = key - '0';
        if ((changeToChannel < triSoup.size()) && (viewChannel != changeToChannel)) {
            viewChannel = changeToChannel;
        }
    }
    else {
        switch (key) {
        case ' ': {
            proceedOptimization();
            viewChannel = channel_result;
            break;
        }

        case '/': {
            toggleOptimization();
            break;
        }

        case 'l':
        case 'L': {
            show_lines = !show_lines;
            viewer.data().show_lines = show_lines;
            break;
        }

        case 'u':
        case 'U': {
            viewUV = !viewUV;
            break;
        }

        case 'd':
        case 'D': {
            showDistortion++;
            if (showDistortion > 2) {
                showDistortion = 0;
            }
            break;
        }

        case 'o':
        case 'O': {
            infoName = std::to_string(iterNum);
            saveInfo(true, false, true);
            break;
        }

        case 'a':
        case 'A': {
            optimizer->checkGradient();
            break;
        }

        default:
            break;
        }
    }

    updateViewerData();

    return false;
}
#endif

#ifdef USE_OPENGL
bool postDrawFunc(igl::opengl::glfw::Viewer& viewer)
#else
bool postDrawFunc()
#endif
{
    if (autoSwitch && (iterNum == 0)) {
        toggleOptimization();
    }

    if (saveInfo_postDraw) {
        saveInfo_postDraw = false;
        saveInfo(outerLoopFinished, true, outerLoopFinished);
        // Note that the content saved in the screenshots are depends on where updateViewerData() is called
        if (outerLoopFinished) {
            //            triSoup[channel_result]->saveAsMesh(outputFolderPath + infoName + "_mesh_01UV.obj", true);
        }
    }

    if (outerLoopFinished) { //var name change!!!
        GifEnd(&GIFWriter);
        saveInfoForPresent();
        saveStats();
        if (autoSwitch) {
            exit(0);
        }
        else {
#ifdef USE_OPENGL
            viewer.core().is_animating = false;
#endif
            outerLoopFinished = false;
        }
    }

    return false;
}

#ifdef USE_OPENGL
void converge_preDrawFunc(igl::opengl::glfw::Viewer& viewer)
#else
void converge_preDrawFunc()
#endif
{
    infoName = "finalResult";

    secPast += difftime(time(NULL), lastStart_world);
    updateViewerData();

    optimization_on = false;
#ifdef USE_OPENGL
    viewer.core().is_animating = false;
#endif
    spdlog::info("optimization converged, with {:d} inner iterations in {:g}s.", optimizer->getInnerIterAmt(), secPast);
    logFile << "optimization converged, with " << optimizer->getInnerIterAmt() << " inner iterations in " << secPast << "s." << std::endl;
    spdlog::critical("simulation finished");
    outerLoopFinished = true;
}

#ifdef USE_OPENGL
bool preDrawFunc(igl::opengl::glfw::Viewer& viewer)
#else
bool preDrawFunc()
#endif
{
    static bool initViewerData = true;
    if (initViewerData) {
        updateViewerData();
        initViewerData = false;
    }

    if (optimization_on) {
        if (offlineMode) {
            while (!converged) {
                proceedOptimization();
            }
        }
        else {
            proceedOptimization();
        }

        //        viewChannel = channel_result;
        updateViewerData();

        if (converged) {
            saveInfo_postDraw = true;

#ifdef USE_OPENGL
            converge_preDrawFunc(viewer);
#else
            converge_preDrawFunc();
#endif
        }
    }
    return false;
}

struct CLIArgs {
    int progMode;
    std::string inputFileName;
    std::string folderTail;
    int logLevel = 0; // trace (all logs)
};

void init_cli_app(CLI::App& app, CLIArgs& args)
{
    app.add_set("progMode", args.progMode, { 0, 10, 11, 100 },
        "program mode 0=optimization, 10=auto-switch, "
        "11=auto-switch (save png), 100=offline",
        true);
    app.add_option("inputFileName", args.inputFileName, "input file name", true);
    app.add_option(
        "folderTail", args.folderTail, "output folder name tail/suffix", true);
    app.add_set("--logLevel,--log", args.logLevel, { 0, 1, 2, 3, 4, 5, 6 },
        "set log level (0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical,"
        " 6=off)",
        true);
}

int main(int argc, char* argv[])
{
#ifdef USE_TBB
#ifdef TBB_NUM_THREADS
    tbb::task_scheduler_init init(TBB_NUM_THREADS);
#endif
#endif

    CLI::App app{ "IPC" };
    CLIArgs args;
    init_cli_app(app, args);
    try {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    spdlog::set_level(static_cast<spdlog::level::level_enum>(args.logLevel));
    if (args.logLevel == 6) {
        std::cout.setstate(std::ios_base::failbit);
    }

    switch (args.progMode) {
    case 0:
        // optimization mode
        spdlog::info("Optimization mode");
        break;

    case 10:
        // autoswitch optimization mode
        autoSwitch = true;
        spdlog::info("Auto-switch optimization mode");
        break;

    case 11:
        // autoswitch optimization mode
        autoSwitch = true;
        savePNG = true;
        spdlog::info("Auto-switch optimization mode");
        break;

    case 100:
        // offline optimization mode
        autoSwitch = true;
        offlineMode = true;
        spdlog::info("Offline optimization mode without visualization");
        break;

    default:
        spdlog::error("No progMode {:d}", args.progMode);
        return 1;
    }

    // Optimization mode
    std::string meshFilePath;
    if (args.inputFileName.at(0) == '/') {
        spdlog::info("The input script file name is gloabl mesh file path.");
        meshFilePath = args.inputFileName;
        args.inputFileName = args.inputFileName.substr(args.inputFileName.find_last_of('/') + 1);
    }
    else {
        meshFilePath = args.inputFileName;
    }
    std::string meshName = args.inputFileName.substr(0, args.inputFileName.find_last_of('.'));
    // Load mesh
    Eigen::MatrixXd V, UV, N;
    Eigen::MatrixXi F, FUV, FN, E;
    const std::string suffix = meshFilePath.substr(meshFilePath.find_last_of('.'));
    bool loadSucceed = false;
    std::vector<std::vector<int>> borderVerts_primitive;
    std::vector<int> componentNodeRange, componentSFRange, componentCERange, componentCoDim;
    std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>> componentMaterial, componentLVels, componentAVels;
    std::vector<std::pair<Eigen::Vector3i, std::array<Eigen::Vector3d, 2>>> componentInitVels;
    std::vector<std::pair<std::vector<int>, std::array<Eigen::Vector3d, 2>>> DBCInfo;
    std::map<int, Eigen::Matrix<double, 1, DIM>> NBCInfo;
    std::vector<std::pair<int, std::string>> meshSeqFolderPath;
    if (suffix == ".txt") {
        loadSucceed = !config.loadFromFile(meshFilePath);
        if (loadSucceed) {
            assert(DIM == 3);

            Eigen::MatrixXd newV;
            Eigen::MatrixXi newF, newSF, newE;

            V.resize(0, 3);
            F.resize(0, 4);
            SF.resize(0, 3);
            E.resize(0, 2);
            componentNodeRange.emplace_back(0);
            componentSFRange.emplace_back(0);
            componentCERange.emplace_back(0);

            int DBCI = 0, NBCI = 0;
            for (int i = 0; i < (int)config.inputShapePaths.size(); ++i) {
                componentCoDim.emplace_back(3);

                const auto& inputShapePath = config.inputShapePaths[i];
                const auto& translate = config.inputShapeTranslates[i];
                const auto& rotate = config.inputShapeRotates[i];
                const auto& scale = config.inputShapeScales[i];
                int suffixI = inputShapePath.find_last_of('.');
                if (suffixI == std::string::npos) {
                    IPC::IglUtils::readNodeEle(inputShapePath, newV, newF, newSF);
                }
                else {
                    const std::string meshFileSuffix = inputShapePath.substr(suffixI);
                    if (meshFileSuffix == ".msh") {
                        if (!IPC::IglUtils::readTetMesh(inputShapePath, newV, newF, newSF)) {
                            spdlog::error("Unable to read input msh file: {:s}", inputShapePath);
                            exit(1);
                        }
                        newE.resize(0, 2);

                        if (config.inputShapeMaterials[i][0] > 0.0 && config.inputShapeMaterials[i][1] > 0.0 && config.inputShapeMaterials[i][2] > 0.0) {
                            Eigen::Vector3i startToEnd;
                            startToEnd[0] = i;
                            startToEnd[1] = F.rows();
                            startToEnd[2] = F.rows() + newF.rows();
                            componentMaterial.emplace_back(startToEnd, config.inputShapeMaterials[i]);
                        }
                    }
                    else if (meshFileSuffix == ".ele") {
                        std::string inputPathPre = inputShapePath.substr(0, suffixI);
                        IPC::IglUtils::readNodeEle(inputPathPre, newV, newF, newSF);
                        newE.resize(0, 2);

                        if (config.inputShapeMaterials[i][0] > 0.0 && config.inputShapeMaterials[i][1] > 0.0 && config.inputShapeMaterials[i][2] > 0.0) {
                            Eigen::Vector3i startToEnd;
                            startToEnd[0] = i;
                            startToEnd[1] = F.rows();
                            startToEnd[2] = F.rows() + newF.rows();
                            componentMaterial.emplace_back(startToEnd, config.inputShapeMaterials[i]);
                        }
                    }
                    else if (meshFileSuffix == ".obj") {
                        // for kinematic object
                        componentCoDim.back() = 2;
                        if (!igl::readOBJ(inputShapePath, newV, newSF)) {
                            spdlog::error("Unable to read input obj file: {:s}", inputShapePath);
                            exit(-1);
                        }
                        newF.resize(0, 4);
                        newE.resize(0, 2);
                    }
                    else if (meshFileSuffix == ".seg") {
                        // for kinematic object
                        componentCoDim.back() = 1;
                        if (!IPC::IglUtils::readSEG(inputShapePath, newV, newE)) {
                            Eigen::MatrixXi tempF;
                            if (!igl::readOBJ(inputShapePath.substr(0, suffixI) + ".obj", newV, tempF)) {
                                spdlog::error("Unable to read input seg or obj file: {:s}", inputShapePath);
                                exit(1);
                            }

                            std::set<std::pair<int, int>> edgesSet;
                            for (int sfI = 0; sfI < tempF.rows(); ++sfI) {
                                auto finder = edgesSet.find(std::pair<int, int>(tempF(sfI, 1), tempF(sfI, 0)));
                                if (finder == edgesSet.end()) {
                                    edgesSet.insert(std::pair<int, int>(tempF(sfI, 0), tempF(sfI, 1)));
                                }

                                finder = edgesSet.find(std::pair<int, int>(tempF(sfI, 2), tempF(sfI, 1)));
                                if (finder == edgesSet.end()) {
                                    edgesSet.insert(std::pair<int, int>(tempF(sfI, 1), tempF(sfI, 2)));
                                }

                                finder = edgesSet.find(std::pair<int, int>(tempF(sfI, 0), tempF(sfI, 2)));
                                if (finder == edgesSet.end()) {
                                    edgesSet.insert(std::pair<int, int>(tempF(sfI, 2), tempF(sfI, 0)));
                                }
                            }

                            newE.resize(edgesSet.size(), 2);
                            int newEI = 0;
                            for (const auto& eI : edgesSet) {
                                newE(newEI, 0) = eI.first;
                                newE(newEI, 1) = eI.second;
                                ++newEI;
                            }
                        }
                        newF.resize(0, 4);
                        newSF.resize(0, 3);
                    }
                    else if (meshFileSuffix == ".pt") {
                        // for kinematic object
                        componentCoDim.back() = 0;
                        Eigen::MatrixXi temp;
                        if (!igl::readOBJ(inputShapePath, newV, temp)) {
                            if (!igl::readOBJ(inputShapePath.substr(0, suffixI) + ".obj", newV, temp)) {
                                spdlog::error("Unable to read input pt or obj file: {:s}", inputShapePath);
                                exit(1);
                            }
                        }
                        newF.resize(0, 4);
                        newSF.resize(0, 3);
                        newE.resize(0, 2);
                    }
                    else {
                        spdlog::error("Unsupported tet mesh file format: {:s}", meshFileSuffix);
                        exit(1);
                    }

                    if (!std::isnan(config.inputShapeLVels[i][0]) && !std::isnan(config.inputShapeLVels[i][1]) && !std::isnan(config.inputShapeLVels[i][2])) {
                        Eigen::Vector3i startToEnd;
                        startToEnd[0] = i;
                        startToEnd[1] = V.rows();
                        startToEnd[2] = V.rows() + newV.rows();
                        componentLVels.emplace_back(startToEnd, config.inputShapeLVels[i]);
                    }
                    if (!std::isnan(config.inputShapeAVels[i][0]) && !std::isnan(config.inputShapeAVels[i][1]) && !std::isnan(config.inputShapeAVels[i][2])) {
                        Eigen::Vector3i startToEnd;
                        startToEnd[0] = i;
                        startToEnd[1] = V.rows();
                        startToEnd[2] = V.rows() + newV.rows();
                        componentAVels.emplace_back(startToEnd, config.inputShapeAVels[i]);
                    }
                    if (!std::isnan(config.inputShapeInitVels[i][0][0]) && !std::isnan(config.inputShapeInitVels[i][0][1]) && !std::isnan(config.inputShapeInitVels[i][0][2]) && !std::isnan(config.inputShapeInitVels[i][1][0]) && !std::isnan(config.inputShapeInitVels[i][1][1]) && !std::isnan(config.inputShapeInitVels[i][1][2])) {
                        Eigen::Vector3i startToEnd;
                        startToEnd[0] = i;
                        startToEnd[1] = V.rows();
                        startToEnd[2] = V.rows() + newV.rows();
                        componentInitVels.emplace_back(startToEnd, config.inputShapeInitVels[i]);
                    }

                    while (DBCI < config.inputShapeDBC.size() && config.inputShapeDBC[DBCI].first == i) {
                        // vertex selection
                        std::vector<int> selectedVerts;
                        IPC::IglUtils::Init_Dirichlet(newV, config.inputShapeDBC[DBCI].second[0],
                            config.inputShapeDBC[DBCI].second[1], selectedVerts);
                        for (auto& i : selectedVerts) {
                            i += V.rows();
                        }
                        if (selectedVerts.size()) {
                            DBCInfo.emplace_back(std::pair<std::vector<int>, std::array<Eigen::Vector3d, 2>>(selectedVerts,
                                { config.inputShapeDBC[DBCI].second[2], config.inputShapeDBC[DBCI].second[3] }));
                        }
                        ++DBCI;
                    }
                    while (NBCI < config.inputShapeNBC.size() && config.inputShapeNBC[NBCI].first == i) {
                        // vertex selection
                        std::vector<int> selectedVerts;
                        IPC::IglUtils::Init_Dirichlet(newV, config.inputShapeNBC[NBCI].second[0],
                            config.inputShapeNBC[NBCI].second[1], selectedVerts);
                        for (const auto& vI : selectedVerts) {
                            NBCInfo[vI + V.rows()] = config.inputShapeNBC[NBCI].second[2].template segment<DIM>(0).transpose();
                        }
                        ++NBCI;
                    }
                }

                int existVrt = (int)V.rows();
                for (int i = 0; i < (int)newV.rows(); ++i) {
                    Eigen::Vector3d p = newV.row(i).transpose();
                    newV.row(i) = (rotate * p.cwiseProduct(scale) + translate).transpose();
                }

                if (newF.rows()) {
                    newF.array() += existVrt;
                }
                if (newSF.rows()) {
                    newSF.array() += existVrt;
                }
                if (newE.rows()) {
                    newE.array() += existVrt;
                }

                auto append = [](auto& dst, const auto& src) {
                    int dstRows = dst.rows();
                    int srcRows = src.rows();
                    int srcCols = src.cols();
                    dst.conservativeResize(dstRows + srcRows, Eigen::NoChange);
                    dst.block(dstRows, 0, srcRows, srcCols) = src;
                };
                append(V, newV);
                if (newF.rows()) {
                    append(F, newF);
                }
                if (newSF.rows()) {
                    append(SF, newSF);
                }
                if (newE.rows()) {
                    append(E, newE);
                }

                componentNodeRange.emplace_back(V.rows());
                componentSFRange.emplace_back(SF.rows());
                componentCERange.emplace_back(E.rows());

                compVAccSize.emplace_back(V.rows());
                compFAccSize.emplace_back(F.rows());
            }

            UV = V.leftCols(DIM);
            if (config.rotDeg != 0.0) {
                const Eigen::Matrix3d rotMtr = Eigen::AngleAxis<double>(config.rotDeg / 180.0 * M_PI,
                    config.rotAxis)
                                                   .toRotationMatrix();
                Eigen::Matrix<double, 1, DIM> center = (UV.colwise().maxCoeff() - UV.colwise().minCoeff()) / 2.0;
#ifdef USE_TBB
                tbb::parallel_for(0, (int)UV.rows(), 1, [&](int vI)
#else
                for (int vI = 0; vI < UV.rows(); ++vI)
#endif
                    {
                        if constexpr (DIM == 3) {
                            UV.row(vI) = (rotMtr * (UV.row(vI) - center).transpose()).transpose() + center;
                        }
                        else {
                            Eigen::Vector3d pos;
                            pos.head(2) = (UV.row(vI) - center).transpose();
                            pos[2] = 0.0;
                            UV.row(vI) = (rotMtr * pos).transpose().head(2) + center;
                        }
                    }
#ifdef USE_TBB
                );
#endif
            }
            if (config.size > 0.0) {
                V *= config.size / (UV.colwise().maxCoeff() - UV.colwise().minCoeff()).maxCoeff();
                UV *= config.size / (UV.colwise().maxCoeff() - UV.colwise().minCoeff()).maxCoeff();
                V.leftCols(DIM).rowwise() -= UV.colwise().minCoeff();
                UV.rowwise() -= UV.colwise().minCoeff();
            }

#ifdef FIRST_TIME_STEP

            Eigen::RowVector3d maxCoords = V.colwise().maxCoeff();
            Eigen::RowVector3d minCoords = V.colwise().minCoeff();
            Eigen::RowVector3d coordsRange = maxCoords - minCoords;
#ifdef USE_TBB
            tbb::parallel_for(0, (int)V.rows(), 1, [&](int vI)
#else
            for (int vI = 0; vI < V.rows(); ++vI)
#endif
                {
                    // nonuniform scale
                    UV(vI, 0) *= 1.0 + 0.1 * (UV(vI, 0) - minCoords[0]) / coordsRange[0];

                    double yRatio = (UV(vI, 1) - minCoords[1]) / coordsRange[1];
                    UV(vI, 1) += yRatio * 0.1 * UV(vI, 2);

                    UV(vI, 2) *= 1.0 - 0.1 * (UV(vI, 2) - minCoords[2]) / coordsRange[2];
                }
#ifdef USE_TBB
            );
#endif

#endif // FIRST_TIME_STEP

            // nonuniform scale
            //                UV.col(0) *= 1.1;
            //                UV.col(1) *= 1.2;
            //                UV.col(2) *= 1.3;
            // shear
            //                UV.col(0) += 0.1 * UV.col(1);

            IPC::IglUtils::findBorderVerts(V, borderVerts_primitive, config.handleRatio);

            IPC::IglUtils::buildSTri2Tet(F, SF, sTri2Tet);
        }
    }
    else {
        spdlog::error("unkown script file format!");
        return -1;
    }
    if (!loadSucceed) {
        spdlog::error("failed to load script/mesh!");
        return -1;
    }
    vertAmt_input = V.rows();

    // construct mesh data structure
    IPC::Mesh<DIM>* temp = new IPC::Mesh<DIM>(V, F, SF, E, UV,
        componentNodeRange, componentSFRange, componentCERange, componentCoDim,
        componentMaterial, componentLVels, componentAVels, componentInitVels, DBCInfo, NBCInfo,
        config.inputShapeMeshSeqFolderPath,
        config.YM, config.PR, config.rho);
    // primitive test cases
    if ((suffix == ".txt") || (suffix == ".primitive")) {
        temp->borderVerts_primitive = borderVerts_primitive;
    }
    triSoup.emplace_back(temp);

    {
        // for output surface mesh
        isSurfNode.resize(0);
        isSurfNode.resize(temp->V.rows(), false);
        for (int tI = 0; tI < SF.rows(); ++tI) {
            isSurfNode[SF(tI, 0)] = true;
            isSurfNode[SF(tI, 1)] = true;
            isSurfNode[SF(tI, 2)] = true;
        }

        tetIndToSurf.resize(0);
        tetIndToSurf.resize(temp->V.rows(), -1);
        surfIndToTet.resize(0);
        surfIndToTet.resize(temp->V.rows(), -1);
        int sVI = 0;
        for (int vI = 0; vI < isSurfNode.size(); ++vI) {
            if (isSurfNode[vI]) {
                tetIndToSurf[vI] = sVI;
                surfIndToTet[sVI] = vI;
                ++sVI;
            }
        }

        V_surf.resize(sVI, 3);
        F_surf.resize(SF.rows(), 3);
        for (int tI = 0; tI < SF.rows(); ++tI) {
            F_surf(tI, 0) = tetIndToSurf[SF(tI, 0)];
            F_surf(tI, 1) = tetIndToSurf[SF(tI, 1)];
            F_surf(tI, 2) = tetIndToSurf[SF(tI, 2)];
        }
    }

    {
        // for output codimensional segment mesh
        isCENode.resize(0);
        isCENode.resize(temp->V.rows(), false);
        for (int ceI = 0; ceI < temp->CE.rows(); ++ceI) {
            isCENode[temp->CE(ceI, 0)] = true;
            isCENode[temp->CE(ceI, 1)] = true;
        }

        tetIndToCE.resize(0);
        tetIndToCE.resize(temp->V.rows(), -1);
        CEIndToTet.resize(0);
        CEIndToTet.resize(temp->V.rows(), -1);
        int ceVI = 0;
        for (int vI = 0; vI < isCENode.size(); ++vI) {
            if (isCENode[vI]) {
                tetIndToCE[vI] = ceVI;
                CEIndToTet[ceVI] = vI;
                ++ceVI;
            }
        }

        V_CE.resize(ceVI, 3);
        F_CE.resize(temp->CE.rows(), 2);
        for (int ceI = 0; ceI < temp->CE.rows(); ++ceI) {
            F_CE(ceI, 0) = tetIndToCE[temp->CE(ceI, 0)];
            F_CE(ceI, 1) = tetIndToCE[temp->CE(ceI, 1)];
        }
    }

    // Method Type
    std::string startDS = "Sim";

    // create output folder
    mkdir(outputFolderPath.c_str(), 0777);
    if ((suffix == ".txt") || (suffix == ".primitive")) {
        config.appendInfoStr(outputFolderPath);
        outputFolderPath += args.folderTail;
    }
    else {
        outputFolderPath += meshName + "_" + startDS + args.folderTail;
    }
    mkdir(outputFolderPath.c_str(), 0777);
    config.backUpConfig(outputFolderPath + "/config.txt");
    for (int coI = 0; coI < config.collisionObjects.size(); ++coI) {
        config.collisionObjects[coI]->saveMesh(outputFolderPath + "/ACO" + std::to_string(coI) + "_0.obj");
    }
    for (int mcoI = 0; mcoI < config.meshCollisionObjects.size(); ++mcoI) {
        config.meshCollisionObjects[mcoI]->saveMesh(outputFolderPath + "/MCO" + std::to_string(mcoI) + "_0.obj");
    }

    // create log file
    outputFolderPath += '/';
    logFile.open(outputFolderPath + "log.txt");
    if (!logFile.is_open()) {
        spdlog::error("failed to create log file, please ensure output directory is created successfully!");
        return -1;
    }

    spdlog::critical("output path: {:s}", outputFolderPath);

    // setup timer
    timer.new_activity("descent");

    timer_step.new_activity("matrixComputation");
    timer_step.new_activity("matrixAssembly");
    timer_step.new_activity("symbolicFactorization");
    timer_step.new_activity("numericalFactorization");
    timer_step.new_activity("backSolve");
    timer_step.new_activity("lineSearch_other");
    timer_step.new_activity("modifyGrad"); //previously boundarySplit
    timer_step.new_activity("modifySearchDir"); //previously interiorSplit
    timer_step.new_activity("updateHistory"); //previously cornerMerge
    timer_step.new_activity("lineSearch_eVal");
    timer_step.new_activity("fullyImplicit_eComp");
    timer_step.new_activity("solve_extraComp");
    timer_step.new_activity("compGrad");
    timer_step.new_activity("CCD");
    timer_step.new_activity("computeConstraintSets");

    timer_temp3.new_activity("init");
    timer_temp3.new_activity("initPrimal");
    timer_temp3.new_activity("initDual");
    timer_temp3.new_activity("initWeights");
    timer_temp3.new_activity("initCons");
    timer_temp3.new_activity("subdSolve");
    timer_temp3.new_activity("consSolve");
    timer_temp3.new_activity("build hash CCS");
    timer_temp3.new_activity("PT CCS");
    timer_temp3.new_activity("EE CCS");
    timer_temp3.new_activity("merge CCS");
    timer_temp3.new_activity("build hash CCD");
    timer_temp3.new_activity("PPET CCD");
    timer_temp3.new_activity("EE CCD");

    timer_mt.new_activity("query EE CCS");
    timer_mt.new_activity("type EE CCS");
    timer_mt.new_activity("dist EE CCS");
    timer_mt.new_activity("query EE CCD");
    timer_mt.new_activity("dist EE CCD");
    timer_mt.new_activity("EV EE CCD");
    timer_mt.new_activity("getE EE CCS");
    timer_mt.new_activity("getE EE CCD");
    timer_mt.new_activity("matSp graph");
    timer_mt.new_activity("matSp allocate");
    timer_mt.new_activity("H_elasticity");
    timer_mt.new_activity("massM");
    timer_mt.new_activity("H_IPC_ACO");
    timer_mt.new_activity("H_IPC_MCO");
    timer_mt.new_activity("H_IPC_Self");
    timer_mt.new_activity("hashP"); // 15
    timer_mt.new_activity("hashE");
    timer_mt.new_activity("hashT");
    timer_mt.new_activity("constructHash");
    timer_mt.new_activity("QE_locate");
    timer_mt.new_activity("QE_loop");
    timer_mt.new_activity("QE_collect");
    timer_mt.new_activity("paraEE CCS");
    timer_mt.new_activity("looping EE CCS");

    // * Our approach
    energyParams.emplace_back(1.0);
    switch (config.energyType) {
    case IPC::ET_NH:
        energyTerms.emplace_back(new IPC::NeoHookeanEnergy<DIM>());
        break;

    case IPC::ET_FCR:
        energyTerms.emplace_back(new IPC::FixedCoRotEnergy<DIM>());
        break;
    }
    //        energyTerms.back()->checkEnergyVal(*triSoup[0]);
    //        energyTerms.back()->checkGradient(*triSoup[0]);
    //        energyTerms.back()->checkHessian(*triSoup[0], true);

    optimizer = new IPC::Optimizer<DIM>(*triSoup[0], energyTerms, energyParams, false, Eigen::MatrixXd(), Eigen::MatrixXi(), Eigen::VectorXi(), config);
    optimizer->setTime(config.duration, config.dt);

    optimizer->precompute();
    double delay_10ms = std::min(10.0, optimizer->getDt() * 100.0);
    GIFStep = static_cast<int>(std::ceil(3.0 / delay_10ms));
    GIFDelay = static_cast<int>(delay_10ms * GIFStep); // always about 3x10ms, around 33FPS
    // set slower motion playback
    GIFStep = static_cast<int>(GIFStep * config.playBackSpeed);
    if (GIFStep < 1) {
        GIFStep = 1;
        logFile << config.playBackSpeed << "X playback can't be reached with 33FPS at dt=" << optimizer->getDt() << std::endl;
        logFile << "using playBackSpeed at " << 1.0 / static_cast<int>(std::ceil(3.0 / delay_10ms)) << "X" << std::endl;
    }
    triSoup.emplace_back(&optimizer->getResult());

    if (config.disableCout) {
        std::cout << "cout will be disabled from now on..." << std::endl;
        spdlog::set_level(spdlog::level::off);
        std::cout.setstate(std::ios_base::failbit);
    }

    if (offlineMode) {
        while (true) {
#ifdef USE_OPENGL
            preDrawFunc(viewer);
            postDrawFunc(viewer);
#else
            preDrawFunc();
            postDrawFunc();
#endif
        }
    }
    else {
#ifdef USE_OPENGL
        // Setup viewer and launch
        viewer.core().background_color << 1.0f, 1.0f, 1.0f, 0.0f;
        viewer.callback_key_down = &key_down;
        viewer.callback_pre_draw = &preDrawFunc;
        viewer.callback_post_draw = &postDrawFunc;
        viewer.data().show_lines = true;
        viewer.core().orthographic = config.orthographic;
        viewer.core().camera_zoom *= config.zoom;
        viewer.core().animation_max_fps = 60.0;
        viewer.data().point_size = 5;
        viewer.data().show_overlay = true;
#if (DIM == 3)
        if (!config.orthographic) {
            viewer.core().trackball_angle = Eigen::Quaternionf(Eigen::AngleAxisf(M_PI_4 / 2.0, Eigen::Vector3f::UnitX()));
        }
#endif
        viewer.launch();
#else
        spdlog::critical("Only offline mode is supported when OpenGL is disabled. See --help for more details.");
#endif
    }

    // Before exit
    logFile.close();
    for (auto& eI : energyTerms) {
        delete eI;
    }
    delete optimizer;
    delete triSoup[0];
}
