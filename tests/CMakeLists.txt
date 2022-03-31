cmake_minimum_required(VERSION 3.1)
################################################################################

add_executable(${PROJECT_NAME}_tests
    main.cpp

    # Test SQP Constraints
    Collisions/CollisionConstraintTests.cpp
)

################################################################################
# Required Libraries
################################################################################

# IPC
target_link_libraries(${PROJECT_NAME}_tests PUBLIC ${PROJECT_NAME}_dev)

# Catch2
include(catch2)
target_link_libraries(${PROJECT_NAME}_tests PUBLIC Catch2::Catch2)

# finite-diff
include(finite_diff)
target_link_libraries(${PROJECT_NAME}_tests PUBLIC finitediff::finitediff)

################################################################################
# Register tests
################################################################################

set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)

# Include Catch2 and provide function `catch_discover_tests` to register tests.
FetchContent_GetProperties(catch2)
include("${catch2_SOURCE_DIR}/contrib/Catch.cmake")

catch_discover_tests(${PROJECT_NAME}_tests)
