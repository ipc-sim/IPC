# GHC Filesystem
if(TARGET ghc::filesystem)
    return()
endif()

message(STATUS "Third-party: creating target 'ghc::filesystem'")

include(FetchContent)
FetchContent_Declare(
    filesystem
    GIT_REPOSITORY https://github.com/gulrak/filesystem.git
    GIT_TAG v1.5.10
)
FetchContent_MakeAvailable(filesystem)

add_library(ghc::filesystem ALIAS ghc_filesystem)
