cmake_minimum_required(VERSION 3.9)

################################################################################
# project-options
option(BUILD_DIAGNOSTIC_PROJECT      "Build the Diagnostic sub-project"     OFF)
option(BUILD_MESH_PROCESSING_PROJECT "Build the MeshProcessing sub-project" ON)
################################################################################

if(BUILD_DIAGNOSTIC_PROJECT)
  add_subdirectory(Diagnostic)
endif()
if(BUILD_MESH_PROCESSING_PROJECT)
  add_subdirectory(MeshProcessing)
endif()
