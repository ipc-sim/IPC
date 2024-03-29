cmake_minimum_required(VERSION 3.9)

add_executable(meshprocessing
  MeshProcessing.cpp
)

set_target_properties(meshprocessing PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
################################################################################
# Required Libraries
################################################################################
target_link_libraries(meshprocessing PUBLIC ${PROJECT_NAME}_dev)

# libigl
target_link_libraries(meshprocessing PUBLIC igl::core igl_copyleft::tetgen)

# tbb
target_link_libraries(meshprocessing PUBLIC TBB::tbb)
