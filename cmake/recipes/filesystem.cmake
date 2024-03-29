# GHC Filesystem
if(TARGET ghc::filesystem)
    return()
endif()

message(STATUS "Third-party: creating target 'ghc::filesystem'")

include(CPM)
CPMAddPackage("gh:gulrak/filesystem@1.5.10")

add_library(ghc::filesystem ALIAS ghc_filesystem)