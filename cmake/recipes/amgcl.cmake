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

if(TARGET amgcl::amgcl)
    return()
endif()

message(STATUS "Third-party: creating target 'amgcl::amgcl'")

include(FetchContent)
FetchContent_Declare(
    amgcl
    GIT_REPOSITORY https://github.com/ddemidov/amgcl.git
    GIT_TAG 1.4.2
    GIT_SHALLOW TRUE
)

function(amgcl_import_target)
    macro(ignore_package NAME VERSION_NUM)
        include(CMakePackageConfigHelpers)
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${NAME}/${NAME}Config.cmake "")
        write_basic_package_version_file(
            ${CMAKE_CURRENT_BINARY_DIR}/${NAME}/${NAME}ConfigVersion.cmake
            VERSION ${VERSION_NUM}
            COMPATIBILITY AnyNewerVersion
        )
        set(${NAME}_DIR ${CMAKE_CURRENT_BINARY_DIR}/${NAME} CACHE PATH "")
        set(${NAME}_ROOT ${CMAKE_CURRENT_BINARY_DIR}/${NAME} CACHE PATH "")
    endmacro()

    include(boost)

    ignore_package(Boost 1.71.0)
    set(Boost_INCLUDE_DIRS "")
    set(Boost_LIBRARIES "")

    # Prefer Config mode before Module mode to prevent lib from loading its own FindXXX.cmake
    set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)

    # Ready to include third-party lib
    FetchContent_MakeAvailable(amgcl)

    target_link_libraries(amgcl INTERFACE Boost::boost)
endfunction()

amgcl_import_target()
