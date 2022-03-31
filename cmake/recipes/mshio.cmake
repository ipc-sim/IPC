# MshIO
# License: Apache-2.0

if(TARGET mshio)
    return()
endif()

message(STATUS "Third-party: creating target 'mshio'")


include(FetchContent)
FetchContent_Declare(
    mshio
    GIT_REPOSITORY https://github.com/qnzhou/MshIO.git
    GIT_TAG a82fb59e25b5e4f586b29ec85a29c3bbc8117307
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(mshio)
