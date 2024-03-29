# OSQP library

if(TARGET osqp::osqp)
    return()
endif()

message(STATUS "Third-party: creating target 'osqp::osqp'")

# Make sure the right types are used
set(DFLOAT OFF CACHE BOOL "Use float numbers instead of doubles"   FORCE)
set(DLONG  OFF CACHE BOOL "Use long integers (64bit) for indexing" FORCE)

include(CPM)
CPMAddPackage("gh:oxfordcontrol/osqp@0.4.1")

if(UNIX AND NOT APPLE)
  set_target_properties(osqpstatic PROPERTIES INTERFACE_LINK_LIBRARIES ${CMAKE_DL_LIBS})
endif()

add_library(osqp::osqp ALIAS osqpstatic)
