cmake_minimum_required(VERSION 3.9)

add_executable(diagnostic
  Diagnostic.cpp
)

set_target_properties(diagnostic PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
################################################################################
# Required Libraries
################################################################################
target_link_libraries(diagnostic PUBLIC ${PROJECT_NAME}_dev)

# Cholmod solver
if(IPC_WITH_CHOLMOD)
  target_link_libraries(diagnostic PUBLIC ${CHOLMOD_LIBRARIES})
endif()

# OSQP
target_link_libraries(diagnostic PUBLIC osqp::osqp)

# libigl
target_link_libraries(diagnostic PUBLIC igl::core)

# tbb
target_link_libraries(diagnostic PUBLIC TBB::tbb)

# AMGCL
if(IPC_WITH_AMGCL)
  target_link_libraries(diagnostic PUBLIC amgcl::amgcl)
endif()
