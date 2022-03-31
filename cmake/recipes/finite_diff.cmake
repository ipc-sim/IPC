if(TARGET finitediff::finitediff)
    return()
endif()

message(STATUS "Third-party: creating target 'finitediff::finitediff'")

include(FetchContent)
FetchContent_Declare(
    finite-diff
    GIT_REPOSITORY https://github.com/zfergus/finite-diff.git
    GIT_TAG c5961a6014b725ac00cfc2df9ecba048766e9b46
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(finite-diff)
