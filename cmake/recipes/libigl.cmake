# libigl (https://github.com/libigl/libigl)
# License: MPL-2.0
if(TARGET igl::core)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::core'")


set(LIBIGL_IMGUI  ${IPC_WITH_OPENGL} CACHE BOOL "" FORCE)
set(LIBIGL_GLFW   ${IPC_WITH_OPENGL} CACHE BOOL "" FORCE)
set(LIBIGL_OPENGL ${IPC_WITH_OPENGL} CACHE BOOL "" FORCE)
set(LIBIGL_STB    ${IPC_WITH_OPENGL} CACHE BOOL "" FORCE)
set(LIBIGL_PREDICATES ON CACHE BOOL "Use exact predicates" FORCE)
set(LIBIGL_COPYLEFT_TETGEN ${IPC_BUILD_MESH_PROCESSING_PROJECT} CACHE BOOL "Use Tetgen" FORCE)

include(eigen)

include(CPM)
CPMAddPackage("gh:libigl/libigl#36930e5d19bed6c7703ec194fed4fa36df945284")