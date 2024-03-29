# spdlog (https://github.com/gabime/spdlog)
# License: MIT
if(TARGET spdlog::spdlog)
    return()
endif()

message(STATUS "Third-party: creating target 'spdlog::spdlog'")

option(SPDLOG_INSTALL "Generate the install target" ON)
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME "spdlog")

include(CPM)
CPMAddPackage("gh:gabime/spdlog@1.9.2")

set_target_properties(spdlog PROPERTIES POSITION_INDEPENDENT_CODE ON)

set_target_properties(spdlog PROPERTIES FOLDER external)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR
   "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(spdlog PRIVATE
        "-Wno-sign-conversion"
    )
endif()
