# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(IPCDownloadExternal)

################################################################################
# Required libraries
################################################################################

message(STATUS "Downloading externals")

# libigl
if(IPC_WITH_PREBUILT_EXT)
  message(STATUS "Building with prebuilt")
  find_library(PREDICATES_LIBRARIES REQUIRED)
  find_path(PREDICATES_INCLUDE_DIR REQUIRED)
  find_path(IGL_INCLUDE_DIR REQUIRED)
  target_include_directories(${PROJECT_NAME}_dev SYSTEM INTERFACE ${PREDICATES_INCLUDE_DIR})
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${PREDICATES_LIBRARIES})
  target_include_directories(${PROJECT_NAME}_dev SYSTEM INTERFACE ${IGL_INCLUDE_DIR})
elseif(NOT TARGET igl)
  download_libigl()
  add_subdirectory(${IPC_EXTERNAL}/libigl EXCLUDE_FROM_ALL)
  # libigl
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::predicates igl::core)
  if(LIBIGL_WITH_TRIANGLE)
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::triangle)
  endif(LIBIGL_WITH_TRIANGLE)
  if(LIBIGL_WITH_TETGEN)
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::tetgen)
  endif(LIBIGL_WITH_TETGEN)
else()
  message(STATUS "Seems like igl is ready")
endif()
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_PREDICATES)

# exact-ccd
if(IPC_WITH_EXACT_CCD)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_EXACT_CCD)
  if(IPC_WITH_PREBUILT_EXT)
    find_path(EXACT_CCD_INCLUDE_DIR REQUIRED)
    find_library(EXACT_CCD_LIBRARIES REQUIRED)
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${EXACT_CCD_LIBRARIES})
    target_include_directories(${PROJECT_NAME}_dev PUBLIC ${EXACT_CCD_INCLUDE_DIR})
  elseif(IPC_WITH_EXACT_CCD AND NOT TARGET exact-ccd::exact-ccd)
    download_exact_ccd()
    add_subdirectory(${IPC_EXTERNAL}/exact-ccd EXCLUDE_FROM_ALL)
    add_library(exact-ccd::exact-ccd ALIAS exact-ccd)
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC exact-ccd::exact-ccd)
  endif()
endif(IPC_WITH_EXACT_CCD)

# spdlog
if(NOT IPC_WITH_PREBUILT_EXT AND NOT TARGET spdlog::spdlog)
  add_library(spdlog INTERFACE)
  add_library(spdlog::spdlog ALIAS spdlog)
  download_spdlog()
  target_include_directories(spdlog SYSTEM INTERFACE ${IPC_EXTERNAL}/spdlog/include)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC spdlog::spdlog)
endif()
if(IPC_WITH_PREBUILT_EXT)
  find_path(SPDLOG_INCLUDE_DIR REQUIRED)
  target_include_directories(${PROJECT_NAME}_dev PUBLIC ${SPDLOG_INCLUDE_DIR})
endif()

################################################################################
# Optional Libraries
################################################################################
if(IPC_WITH_OPENGL AND NOT IPC_WITH_PREBUILT_EXT)
  message(FATAL_ERROR "We are not supporting OpenGL for now")
  #target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::opengl_glfw igl::opengl_glfw_imgui igl::png)
  #target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_OPENGL)
endif()

# tbb
if(IPC_WITH_TBB AND NOT IPC_WITH_PREBUILT_EXT)
  download_tbb()
  set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
  set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  add_subdirectory(${IPC_EXTERNAL}/tbb EXCLUDE_FROM_ALL)
  add_library(TBB::tbb ALIAS tbb_static)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_TBB)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC TBB::tbb)
endif()

# AMGCL
if(IPC_WITH_AMGCL AND NOT IPC_WITH_PREBUILT_EXT)
  download_amgcl()
  set(Boost_USE_MULTITHREADED TRUE)
  add_subdirectory(${IPC_EXTERNAL}/amgcl EXCLUDE_FROM_ALL)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC amgcl::amgcl)
endif()

if(IPC_WITH_GUROBI AND NOT IPC_WITH_PREBUILT_EXT)
  download_eigen_gurobi()
  add_subdirectory(${IPC_EXTERNAL}/eigen-gurobi EXCLUDE_FROM_ALL)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC EigenGurobi::EigenGurobi)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_GUROBI)
endif()

# Rational CCD
if(IPC_WITH_RATIONAL_CCD AND NOT IPC_WITH_PREBUILT_EXT)
  download_rational_ccd()
  add_subdirectory(${IPC_EXTERNAL}/rational_ccd)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC RationalCCD)
endif()

################################################################################
# Parameter parsing library for main
################################################################################
if(IPC_WITH_MAIN)
  message(STATUS "Linking main with parameter parsing libs")
  if(IPC_WITH_GFLAGS)
    message(STATUS "Linking main with GFLAGS")
    find_library(GFLAGS_LIBRARIES REQUIRED)
    find_path(GFLAGS_INCLUDE_DIR REQUIRED)
    target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${GFLAGS_LIBRARIES})
    target_include_directories(${PROJECT_NAME}_bin PUBLIC ${GFLAGS_INCLUDE_DIR})
  else()
    message(STATUS "Linking main with CLI11")
    download_cli11()
    add_subdirectory(${IPC_EXTERNAL}/cli11)
    target_link_libraries(${PROJECT_NAME}_bin PUBLIC CLI11::CLI11)
  endif()
endif(IPC_WITH_MAIN)
