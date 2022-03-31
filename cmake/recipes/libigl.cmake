#
# Copyright 2020 Adobe. All rights reserved.
# This file is licensed to you under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may obtain a copy
# of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under
# the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR REPRESENTATIONS
# OF ANY KIND, either express or implied. See the License for the specific language
# governing permissions and limitations under the License.
#
if(TARGET igl::core)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::core'")

# libigl
# WARNING: LIBIGL_USE_STATIC_LIBRARY=ON does not work with C++ 17
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK              "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM               "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    ON)
option(LIBIGL_WITH_PNG               "Use PNG"                      ON)
option(LIBIGL_WITH_PYTHON            "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            ON)
option(LIBIGL_QUIET_VIEWER           "Suppress viewer messages"     ON)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)

set(LIBIGL_WITH_VIEWER ${IPC_WITH_OPENGL} CACHE BOOL "" FORCE)
set(LIBIGL_WITH_OPENGL_GLFW_IMGUI ${IPC_WITH_OPENGL} CACHE BOOL "" FORCE)
set(LIBIGL_WITH_OPENGL_GLFW ${IPC_WITH_OPENGL} CACHE BOOL "" FORCE)
set(LIBIGL_WITH_OPENGL ${IPC_WITH_OPENGL} CACHE BOOL "" FORCE)
set(LIBIGL_WITH_PNG ${IPC_WITH_OPENGL} CACHE BOOL "" FORCE)

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.3.0
    GIT_SHALLOW TRUE
)
FetchContent_GetProperties(libigl)
if(libigl_POPULATED)
    return()
endif()
FetchContent_Populate(libigl)

include(eigen)

set(LIBIGL_WITH_PREDICATES ON CACHE BOOL "Use exact predicates" FORCE)

list(APPEND CMAKE_MODULE_PATH ${libigl_SOURCE_DIR}/cmake)
include(${libigl_SOURCE_DIR}/cmake/libigl.cmake ${libigl_BINARY_DIR})

# Install rules
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME libigl)
set_target_properties(igl PROPERTIES EXPORT_NAME core)
install(DIRECTORY ${libigl_SOURCE_DIR}/include/igl DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS igl igl_common EXPORT Libigl_Targets)
install(EXPORT Libigl_Targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/igl NAMESPACE igl::)
