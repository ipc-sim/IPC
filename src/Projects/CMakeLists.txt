cmake_minimum_required(VERSION 3.9)

################################################################################
# project-options
option(IPC_BUILD_DIAGNOSTIC_PROJECT      "Build the Diagnostic sub-project"     ON)
option(IPC_BUILD_MESH_PROCESSING_PROJECT "Build the MeshProcessing sub-project" ON)
################################################################################

if(IPC_BUILD_DIAGNOSTIC_PROJECT)
  add_subdirectory(Diagnostic)
endif()
if(IPC_BUILD_MESH_PROCESSING_PROJECT)
  add_subdirectory(MeshProcessing)
endif()
