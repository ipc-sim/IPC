# Eigen Gurobi Wrapper
if(TARGET EigenGurobi::EigenGurobi)
    return()
endif()

message(STATUS "Third-party: creating target 'EigenGurobi::EigenGurobi'")

include(CPM)
CPMAddPackage("gh:zfergus/eigen-gurobi#51b1aacb3c5733555d09fe362887d618ee97826d")

add_library(EigenGurobi::EigenGurobi ALIAS EigenGurobi)
