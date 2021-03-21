cmake_minimum_required(VERSION 3.9)
cmake_policy(SET CMP0025 NEW)
# use clang for compilation
set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_CXX_FLAGS "-std=c++14 -fopenmp")# -D_GLIBCXX_USE_CXX11_ABI=0")

# the following 2 are for fixing undefined symbol caused by CMake unable to find
# correct ABI implementation. See the source of this fix:
# https://gitlab.kitware.com/cmake/cmake/-/issues/18275
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")

project(IPC)

message("Your C++ compiler supports these C++ features:")
message("You compiler: ${CMAKE_CXX_COMPILER}")
foreach(i ${CMAKE_CXX_COMPILE_FEATURES})
  message("${i}")
endforeach()

################################################################################

# libigl
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" ON)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK              "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM               "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG               "Use PNG"                      ON)
option(LIBIGL_WITH_PYTHON            "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            OFF)
option(LIBIGL_QUIET_VIEWER           "Suppress viewer messages"     ON)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)

# project-options
option(IPC_WITH_EXACT_CCD "Use exact CCD as a verification"                     ON)
option(IPC_WITH_TESTS     "Build unit tests using Catch2"                       OFF)
option(IPC_WITH_GUROBI    "Enable Gurobi QP solver for SQP method"             OFF)
option(IPC_WITH_OSQP_MKL  "Use MKL Pardiso linear solver in OSQP if available"  OFF)
option(IPC_WITH_OPENGL    "Use libgil viewer for preview and output"            OFF)
set(IPC_LINSYSSOLVER "CHOLMOD" CACHE STRING "Linear system solver to use (options: CHOLMOD, AMGCL, EIGEN)")
set_property(CACHE IPC_LINSYSSOLVER PROPERTY STRINGS "CHOLMOD" "AMGCL" "EIGEN")
option(IPC_EXIT_UPON_QP_FAIL  "Exit if the QP/SQP intersects or blows-up" OFF)
option(IPC_WITH_MAIN                 "Build binary with library"     ON)
option(IPC_WITH_AMGCL                  "Build with AMGCL"              OFF)
option(IPC_WITH_TBB                    "Build with TBB"                OFF)
option(IPC_WITH_PREBUILT_EXT "For externals, use prebuilt version" OFF)
option(IPC_WITH_RATIONAL_CCD "With rational CCD if turned on when not using prebuilt externals" OFF)

if(NOT IPC_WITH_OPENGL)
  set(LIBIGL_WITH_VIEWER OFF CACHE BOOL "" FORCE)
  set(LIBIGL_WITH_OPENGL_GLFW_IMGUI OFF CACHE BOOL "" FORCE)
  set(LIBIGL_WITH_OPENGL_GLFW OFF CACHE BOOL "" FORCE)
  set(LIBIGL_WITH_OPENGL OFF CACHE BOOL "" FORCE)
  set(LIBIGL_WITH_PNG OFF CACHE BOOL "" FORCE)
endif()

# debug-options

################################################################################
# IPC Library
################################################################################
message("Creating library")

# project source files
file(GLOB SRCFILES
  "ipc/*.cpp"
  "ipc/Energy/*.cpp"
  "ipc/Energy/Physics_Elasticity/*.cpp"
  "ipc/Energy/Geometry_Isometry/*.cpp"
  "ipc/Energy/Connectivity/*.cpp"
  "ipc/Energy/Collision/*.cpp"
  "ipc/CollisionObject/*.cpp"
  "ipc/CCD/*.cpp"
  "ipc/Utils/*.cpp"
  "ipc/LinSysSolver/EigenLibSolver.cpp"
  "ipc/LinSysSolver/CHOLMODSolver.cpp"
#  "src/LinSysSolver/AMGCLSolver.cpp"
  "ipc/TimeStepper/*.cpp"
)

FOREACH(item ${SRCFILES})
    IF(${item} MATCHES "main.cpp")
        LIST(REMOVE_ITEM SRCFILES ${item})
    ENDIF(${item} MATCHES "main.cpp")
ENDFOREACH(item)

add_library(${PROJECT_NAME}_dev ${SRCFILES})
if(IPC_WITH_MAIN)
  add_executable(${PROJECT_NAME}_bin "ipc/main.cpp")
  # Link against the library
  target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_dev)
endif(IPC_WITH_MAIN)

target_include_directories(${PROJECT_NAME}_dev PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#target_include_directories(${PROJECT_NAME}_dev PUBLIC
#  "src/ipc"
#  "src/ipc/Energy"
#  "src/ipc/Energy/Physics_Elasticity"
#  "src/ipc/Energy/Geometry_Isometry"
#  "src/ipc/Energy/Connectivity"
#  "src/ipc/Energy/Collision"
#  "src/ipc/CollisionObject/"
#  "src/ipc/CCD/"
#  "src/ipc/Utils"
#  "src/ipc/Utils/SVD"
#  "src/ipc/LinSysSolver"
#  "src/ipc/TimeStepper"
#)

#if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#  target_link_libraries(${PROJECT_NAME}_dev PUBLIC stdc++fs)
#else()
#  target_link_libraries(${PROJECT_NAME}_dev PUBLIC c++fs)
#endif()

# Link with eigen when using prebuilt. If not with prebuilt externals, library will
# be linked to igl::core in IPCDependencies.cmake which links to Eigen
if(IPC_WITH_PREBUILT_EXT)
  message(STATUS "Linking eigen")
  find_path(EIGEN_INCLUDE_DIR REQUIRED)
  target_include_directories(${PROJECT_NAME}_dev PUBLIC ${EIGEN_INCLUDE_DIR})
endif()

################################################################################

### Configuration
if(IPC_WITH_PREBUILT_EXT)
  message(STATUS "Using prebuilt external libs")
else()
  message(STATUS "Using external downloads")
  set(IPC_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# include(Warnings)
# The following file handles situations for both with and without prebuilt externals
include(${PROJECT_NAME}Dependencies)

################################################################################
# Etienne Vouga's CTCD Library
################################################################################
add_library(EVCTCD "ipc/CCD/EVCTCD/CTCD.cpp")
target_include_directories(EVCTCD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(NOT IPC_WITH_PREBUILT_EXT)
    # Link against Eigen in libigl
    target_link_libraries(EVCTCD PUBLIC igl::core)
else()
    # find and link libigl outside which has Eigen
    target_include_directories(EVCTCD PUBLIC ${EIGEN_INCLUDE_DIR})
    message(STATUS "Test eigen")
    foreach(subdir ${SUBDIR})
      message(STATUS "${subdir}")
    endforeach()
endif()
# Turn of floating point contraction for CCD robustness
#target_compile_options(EVCTCD PUBLIC "-ffp-contract=off")
target_link_libraries(${PROJECT_NAME}_dev PUBLIC EVCTCD)

################################################################################
# Packages that can potentially be found by invoking ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Find${PACKAGENAME}.cmake
# using find_package. This also means by invocation of find_package, a corresponding
# FindPACKAGENAME.cmake is needed (at least in module mode)
################################################################################
# SuiteSparse
set(SUITESPARSE_INCLUDE_DIR_HINTS $ENV{SUITESPARSE_INC})
set(SUITESPARSE_LIBRARY_DIR_HINTS $ENV{SUITESPARSE_LIB})
find_package(SuiteSparse REQUIRED)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${SUITESPARSE_LIBRARIES})
target_include_directories(${PROJECT_NAME}_dev PUBLIC ${SUITESPARSE_INCLUDE_DIRS})

# Set version of SVD to use
message(STATUS "Using scalar implmentation of SVD")
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC -DUSE_SCALAR_IMPLEMENTATION)

# Figure out SSE level support
message(STATUS "Seaching for SSE...")
find_package(SSE)

# Figure out AVX level support
message(STATUS "Searching for AVX...")
find_package(AVX)

# Figure out FMA level support
message(STATUS "Searching for FMA...")
find_package(FMA)

# Add SSE, AVX, and FMA flags to compiler flags
string(REPLACE " " ";" SIMD_FLAGS "${SSE_FLAGS} ${AVX_FLAGS} ${FMA_FLAGS}")
target_compile_options(${PROJECT_NAME}_dev PUBLIC ${SIMD_FLAGS})

# Add -pthread to compilation and linking
find_package(Threads)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${CMAKE_THREAD_LIBS_INIT})

# code only support 3D
message(STATUS "Building for three dimensions")
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC DIM=3)

# Ensure that PRId64 is defined
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC -D__STDC_FORMAT_MACROS)

if(IPC_WITH_RATIONAL_CCD)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC WITH_RATIONAL_CCD)
endif()

#TODO: Choosing linear solver should not be at compile time but rather at
# Set the selected linear system solver
if(IPC_LINSYSSOLVER STREQUAL "CHOLMOD")
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC LINSYSSOLVER_USE_CHOLMOD)
elseif(IPC_LINSYSSOLVER STREQUAL "AMGCL")
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC LINSYSSOLVER_USE_AMGCL)
elseif(NOT (IPC_LINSYSSOLVER STREQUAL "EIGEN"))
  message(WARNING "Invalid linear system solver: ${IPC_LINSYSSOLVER}. Falling back to EIGEN.")
endif()

if(IPC_EXIT_UPON_QP_FAIL)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC EXIT_UPON_QP_FAIL)
endif()

set(LIB_INSTALL_DIR "libs")
set(INCLUDE_INSTALL_DIR "include")
set(BINARY_INSTALL_DIR "bin")

install(DIRECTORY "ipc"
        DESTINATION ${INCLUDE_INSTALL_DIR}
        FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp")
install(TARGETS "EVCTCD" "${PROJECT_NAME}_dev"
        LIBRARY DESTINATION ${LIB_INSTALL_DIR})
if(IPC_WITH_MAIN)
  install(TARGETS "${PROJECT_NAME}_bin"
          RUNTIME DESTINATION ${BIN_INSTALL_DIR})
endif()
