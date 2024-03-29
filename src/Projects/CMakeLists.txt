cmake_minimum_required(VERSION 3.9)

if(IPC_BUILD_DIAGNOSTIC_PROJECT)
  add_subdirectory(Diagnostic)
endif()

if(IPC_BUILD_MESH_PROCESSING_PROJECT)
  add_subdirectory(MeshProcessing)
endif()
