# MshIO (https://github.com/qnzhou/MshIO)
# License: Apache-2.0

if(TARGET mshio)
    return()
endif()

message(STATUS "Third-party: creating target 'mshio'")


include(CPM)
CPMAddPackage("gh:qnzhou/MshIO#29d0263b45bbbb2931ecbe892d0d7f0f3a493d0c")