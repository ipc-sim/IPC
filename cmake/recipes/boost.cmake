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

if(TARGET Boost::boost)
    return()
endif()

message(STATUS "Third-party: creating targets 'Boost::boost'...")

include(FetchContent)
FetchContent_Declare(
    boost-cmake
    GIT_REPOSITORY https://github.com/Orphis/boost-cmake.git
    GIT_TAG 7f97a08b64bd5d2e53e932ddf80c40544cf45edf
)

set(PREVIOUS_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(OLD_CMAKE_POSITION_INDEPENDENT_CODE ${CMAKE_POSITION_INDEPENDENT_CODE})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# This guy will download boost using FetchContent
FetchContent_GetProperties(boost-cmake)
if(NOT boost-cmake_POPULATED)
    FetchContent_Populate(boost-cmake)
    # File lcid.cpp from Boost_locale.cpp doesn't compile on MSVC, so we exclude them from the default
    # targets being built by the project (only targets explicitly used by other targets will be built).
    add_subdirectory(${boost-cmake_SOURCE_DIR} ${boost-cmake_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ${OLD_CMAKE_POSITION_INDEPENDENT_CODE})
set(CMAKE_CXX_FLAGS "${PREVIOUS_CMAKE_CXX_FLAGS}")
