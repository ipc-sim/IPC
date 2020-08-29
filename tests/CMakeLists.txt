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
target_link_libraries(${PROJECT_NAME}_tests PUBLIC Catch2::Catch2)

# finite-diff
target_link_libraries(${PROJECT_NAME}_tests PUBLIC FiniteDiff)

################################################################################
# Add Tests
################################################################################

foreach(source IN ITEMS ${test_sources})
   source_group("tests" FILES "${source}")
endforeach()

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(${PROJECT_NAME}_tests)
