cmake_minimum_required(VERSION 3.9)
project(IPC)

################################################################################

# Detects whether this is a top-level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(IPC_TOPLEVEL_PROJECT ON)
else()
    set(IPC_TOPLEVEL_PROJECT OFF)
endif()

# Project options
option(IPC_WITH_OPENGL "Use libgil viewer for preview and output" ON)
option(IPC_WITH_AMGCL "Enable AMGCL linear system solver (select with IPC_LINSYSSOLVER=\"AMGCL\")" ON)
option(IPC_WITH_CHOLMOD "Enable CHOLMOD linear system solver (select with IPC_LINSYSSOLVER=\"CHOLMOD\")" ON)
set(IPC_LINSYSSOLVER "CHOLMOD" CACHE STRING "Linear system solver to use (options: CHOLMOD, AMGCL, EIGEN)")
set_property(CACHE IPC_LINSYSSOLVER PROPERTY STRINGS "CHOLMOD" "AMGCL" "EIGEN")
## SQP and QP options
option(IPC_WITH_GUROBI       "Enable Gurobi QP solver for SQP method"             OFF)
option(IPC_WITH_OSQP_MKL     "Use MKL Pardiso linear solver in OSQP if available"  ON)
option(IPC_EXIT_UPON_QP_FAIL "Exit if the QP/SQP intersects or blows-up"          OFF)
## CCD options
option(IPC_WITH_EXACT_CCD "Use exact CCD as a verification"                OFF)
option(IPC_WITH_FPRP      "Use exact floating-point CCD as a verification" OFF)
mark_as_advanced(IPC_WITH_FPRP)
## Top level options
option(IPC_WITH_TESTS "Build unit tests using Catch2" ${IPC_TOPLEVEL_PROJECT})

# libigl
# WARNING: LIBIGL_USE_STATIC_LIBRARY=ON does not work with C++ 17
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK              "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM               "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    ON)
option(LIBIGL_WITH_PNG               "Use PNG"                      ON)
option(LIBIGL_WITH_PYTHON            "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            ON)
option(LIBIGL_QUIET_VIEWER           "Suppress viewer messages"     ON)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)

if(NOT IPC_WITH_OPENGL)
  set(LIBIGL_WITH_VIEWER OFF CACHE BOOL "" FORCE)
  set(LIBIGL_WITH_OPENGL_GLFW_IMGUI OFF CACHE BOOL "" FORCE)
  set(LIBIGL_WITH_OPENGL_GLFW OFF CACHE BOOL "" FORCE)
  set(LIBIGL_WITH_OPENGL OFF CACHE BOOL "" FORCE)
  set(LIBIGL_WITH_PNG OFF CACHE BOOL "" FORCE)
endif()

################################################################################

### Configuration
set(IPC_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# include(Warnings)
include(${PROJECT_NAME}Dependencies)

################################################################################
# IPC Library
################################################################################

# project source files
file(GLOB SRCFILES
  "src/*.cpp"
  "src/Energy/*.cpp"
  "src/Energy/Physics_Elasticity/*.cpp"
  "src/Energy/Geometry_Isometry/*.cpp"
  "src/Energy/Connectivity/*.cpp"
  "src/Energy/Collision/*.cpp"
  "src/CollisionObject/*.cpp"
  "src/CCD/*.cpp"
  "src/Utils/*.cpp"
  "src/LinSysSolver/EigenLibSolver.cpp"
  "src/LinSysSolver/CHOLMODSolver.cpp"
  "src/LinSysSolver/AMGCLSolver.cpp"
  "src/TimeStepper/*.cpp"
)

FOREACH(item ${SRCFILES})
    IF(${item} MATCHES "main.cpp")
        LIST(REMOVE_ITEM SRCFILES ${item})
    ENDIF(${item} MATCHES "main.cpp")
ENDFOREACH(item)

add_library(${PROJECT_NAME}_dev ${SRCFILES})

target_include_directories(${PROJECT_NAME}_dev PUBLIC
  "src"
  "src/Energy"
  "src/Energy/Physics_Elasticity"
  "src/Energy/Geometry_Isometry"
  "src/Energy/Connectivity"
  "src/Energy/Collision"
  "src/CollisionObject/"
  "src/CCD/"
  "src/Utils"
  "src/Utils/SVD"
  "src/LinSysSolver"
  "src/TimeStepper"
)

################################################################################
# Simulations Executable
################################################################################

add_executable(${PROJECT_NAME}_bin "src/main.cpp")
# Link against the library
target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_dev)
# CLI11
target_link_libraries(${PROJECT_NAME}_bin PUBLIC CLI11::CLI11)

################################################################################
# Required Libraries
################################################################################

# SuiteSparse
if(IPC_WITH_CHOLMOD)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${SUITESPARSE_LIBRARIES})
  target_include_directories(${PROJECT_NAME}_dev PUBLIC ${SUITESPARSE_INCLUDE_DIRS})
endif()

# OSQP
target_link_libraries(${PROJECT_NAME}_dev PUBLIC osqp::osqp)

# libigl
target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::core igl::triangle igl::tetgen)
if(IPC_WITH_OPENGL)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::opengl_glfw igl::opengl_glfw_imgui igl::png)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_OPENGL)
endif()

if(LIBIGL_WITH_PREDICATES)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::predicates)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_PREDICATES)
endif()

# CCD
target_link_libraries(${PROJECT_NAME}_dev PUBLIC CCDWrapper)
if(IPC_WITH_EXACT_CCD)
  message(WARNING "The exact CCD methods RootParity and BSC miss some collisions! Please use TightInclusion instead.")
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_EXACT_CCD)
endif()
if(IPC_WITH_FPRP)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_FPRP_CCD)
endif()

# tbb
target_link_libraries(${PROJECT_NAME}_dev PUBLIC TBB::tbb)

# Logger
target_link_libraries(${PROJECT_NAME}_dev PUBLIC spdlog::spdlog)

# AMGCL
if(IPC_WITH_AMGCL)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC amgcl::amgcl)
endif()

# target_link_libraries(${PROJECT_NAME}_dev PUBLIC pardiso600-GNU720-X86-64 gfortran)
# find_package(Pardiso)
# if(PARDISO_FOUND)
#     find_package(LAPACK)
#     if(LAPACK_FOUND)
#         target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${LAPACK_LIBRARIES})
#     else()
#         message(FATAL_ERROR “unable to find lapack”)
#     endif()

#     # find_package(OpenMP)
#     # if( OpenMP_CXX_FOUND )
#     #   target_link_libraries(polyfem PUBLIC ${OpenMP_CXX_LIBRARIES})
#     #   # target_compile_definitions(polyfem PUBLIC ${OpenMP_CXX_FLAGS})
#     # else()
#     #   message(FATAL_ERROR “unable to find omp”)
#     # endif()

#     target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${PARDISO_LIBRARIES})
# else()
#     message(WARNING “Pardiso not found, solver will not be available.“)
# endif()
# target_link_libraries(${PROJECT_NAME}_dev PUBLIC gfortran)

if(IPC_WITH_GUROBI)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC EigenGurobi::EigenGurobi)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_GUROBI)
endif()

# Extra warnings
# target_link_libraries(${PROJECT_NAME}_dev PRIVATE warnings::all)

# MshIO
target_link_libraries(${PROJECT_NAME}_dev PUBLIC mshio::mshio)

# Filesystem
target_link_libraries(${PROJECT_NAME}_dev PUBLIC ghc::filesystem)

################################################################################
# Sub-Projects
################################################################################

add_subdirectory(src/Projects)

################################################################################
# Code Format Tool
################################################################################

find_program(CLANG_FORMAT_INSTALLED "clang-format")
if(CLANG_FORMAT_INSTALLED)
  message(STATUS "Building code format tool")
  add_subdirectory(Format)
else()
  message(STATUS "clang-format not found: skipping Format subdirectory")
endif()


################################################################################
# Compiler options
################################################################################

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

# Use C++17
target_compile_features(${PROJECT_NAME}_dev PUBLIC cxx_std_17)

# code only support 3D
message(STATUS "Building for three dimensions")
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC DIM=3)

# Ensure that PRId64 is defined
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC -D__STDC_FORMAT_MACROS)

# Unit tests
if(IPC_WITH_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()

# MKL Pardiso OSQP linear system solver
if(IPC_WITH_OSQP_MKL)
  find_library(MKL_RT_LIBRARY mkl_rt
               PATHS $ENV{MKLROOT}/lib
                     $ENV{MKLROOT}/lib/intel64
                     $ENV{INTEL}/mkl/lib/intel64
                     /opt/intel/mkl/lib)
  if(MKL_RT_LIBRARY)
    message(STATUS "MKL found: enabling OSQP MKL Pardiso linear system solver")
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC OSQP_USE_MKL_PARDISO)
  else()
    message(WARNING "MKL not found: disabling OSQP MKL Pardiso linear system solver")
  endif()
endif()

# Set the selected linear system solver
if(IPC_LINSYSSOLVER STREQUAL "CHOLMOD")
  if(IPC_WITH_CHOLMOD)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC LINSYSSOLVER_USE_CHOLMOD)
  else()
    message(FATAL_ERROR "CHOLMOD selected as IPC_LINSYSSOLVER, but CHOLMOD is disabled (IPC_WITH_CHOLMOD=OFF)!")
  endif()
elseif(IPC_LINSYSSOLVER STREQUAL "AMGCL")
  if(IPC_WITH_AMGCL)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC LINSYSSOLVER_USE_AMGCL)
  else()
    message(FATAL_ERROR "AMGCL selected as IPC_LINSYSSOLVER, but AMGCL is disabled (IPC_WITH_AMGCL=OFF)!")
  endif()
elseif(NOT (IPC_LINSYSSOLVER STREQUAL "EIGEN"))
  message(WARNING "Invalid linear system solver: ${IPC_LINSYSSOLVER}. Falling back to EIGEN.")
endif()

if(IPC_WITH_AMGCL)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_AMGCL)
endif()

if(IPC_WITH_CHOLMOD)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_CHOLMOD)
endif()

# Quiet the viewer
if(LIBIGL_QUIET_VIEWER)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC IGL_VIEWER_VIEWER_QUIET)
endif()

if(IPC_EXIT_UPON_QP_FAIL)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC EXIT_UPON_QP_FAIL)
endif()
