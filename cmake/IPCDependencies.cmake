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
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::predicates igl::core igl::triangle igl::tetgen)
else()
  message(STATUS "Seems like igl is ready")
endif()
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_PREDICATES)

# TBB
if(IPC_WITH_PREBUILT_EXT)
  message(STATUS "TBB is skipped, we won't use TBB when using prebuilt externals")
elseif(IPC_WITH_TBB AND NOT TARGET TBB::tbb)
  download_tbb()
  set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
  set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  add_subdirectory(${IPC_EXTERNAL}/tbb EXCLUDE_FROM_ALL)
  add_library(TBB::tbb ALIAS tbb_static)
endif()

# exact-ccd

elseif(IPC_WITH_EXACT_CCD AND NOT TARGET exact-ccd::exact-ccd)
  download_exact_ccd()
  add_subdirectory(${IPC_EXTERNAL}/exact-ccd EXCLUDE_FROM_ALL)
  add_library(exact-ccd::exact-ccd ALIAS exact-ccd)
endif()

# spdlog
if(NOT IPC_WITH_PREBUILT_EXT AND NOT TARGET spdlog::spdlog)
  add_library(spdlog INTERFACE)
  add_library(spdlog::spdlog ALIAS spdlog)
  download_spdlog()
  target_include_directories(spdlog SYSTEM INTERFACE ${IPC_EXTERNAL}/spdlog/include)
endif()
if(IPC_WITH_PREBUILT_EXT)
  find_path(SPDLOG REQUIRED)
  target_include_directories(spdlog SYSTEM INTERFACE ${SPDLOG_INCLUDE_DIR})
endif()

# CLI11
if(NOT IPC_WITH_GFLAGS AND NOT TARGET CLI11::CLI11)
  download_cli11()
  add_subdirectory(${IPC_EXTERNAL}/cli11)
endif()

# eigen-gurobi
if(IPC_WITH_GUROBI AND NOT TARGET EigenGurobi::EigenGurobi)
  download_eigen_gurobi()
  add_subdirectory(${IPC_EXTERNAL}/eigen-gurobi EXCLUDE_FROM_ALL)
  add_library(EigenGurobi::EigenGurobi ALIAS EigenGurobi)
endif()

# Rational CCD
if(IPC_WITH_RATIONAL_CCD)
  download_rational_ccd()
  add_subdirectory(${IPC_EXTERNAL}/rational_ccd)
endif()


################################################################################
# Required Libraries
################################################################################

if(IPC_WITH_OPENGL)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::opengl_glfw igl::opengl_glfw_imgui igl::png)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_OPENGL)
endif()

if(IPC_WITH_RATIONAL_CCD)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC RationalCCD)
endif()

# tbb
if(IPC_WITH_TBB)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_TBB)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC TBB::tbb)
endif()

# exact-ccd
if(IPC_WITH_EXACT_CCD)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC exact-ccd::exact-ccd)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_EXACT_CCD)
endif()

# Logger
if(IPC_WITH_PREBUILT_EXT)
  find_package(spdlog REQUIRED)
endif()
target_link_libraries(${PROJECT_NAME}_dev PUBLIC spdlog::spdlog)

# AMGCL
if(IPC_WITH_AMGCL)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC amgcl::amgcl)
endif()

if(IPC_WITH_GUROBI)
  target_link_libraries(${PROJECT_NAME}_dev PUBLIC EigenGurobi::EigenGurobi)
  target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_GUROBI)
endif()
