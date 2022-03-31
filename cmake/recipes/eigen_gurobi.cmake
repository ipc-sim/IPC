# Eigen Gurobi Wrapper
if(TARGET EigenGurobi::EigenGurobi)
    return()
endif()

message(STATUS "Third-party: creating target 'EigenGurobi::EigenGurobi'")

include(FetchContent)
FetchContent_Declare(
    eigen_gurobi
    GIT_REPOSITORY https://github.com/zfergus/eigen-gurobi.git
    GIT_TAG 51b1aacb3c5733555d09fe362887d618ee97826d
)
FetchContent_MakeAvailable(eigen_gurobi)

add_library(EigenGurobi::EigenGurobi ALIAS EigenGurobi)
