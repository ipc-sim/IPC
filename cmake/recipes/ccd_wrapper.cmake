# CCD Wrapper (includes Rational CCD)

if(TARGET ccd_wrapper::ccd_wrapper)
    return()
endif()

message(STATUS "Third-party: creating target 'ccd_wrapper::ccd_wrapper'")

option(CCD_WRAPPER_WITH_FPRF                              "Enable floating-point root finder method"       ON)
set(CCD_WRAPPER_WITH_RP  ${IPC_WITH_EXACT_CCD} CACHE BOOL "Enable root parity method"                   FORCE)
set(CCD_WRAPPER_WITH_BSC ${IPC_WITH_EXACT_CCD} CACHE BOOL "Enable Bernstein sign classification method" FORCE)
option(CCD_WRAPPER_WITH_TIGHT_INCLUSION                   "Enable Tight Inclusion method"                  ON)
option(TIGHT_INCLUSION_WITH_NO_ZERO_TOI                   "Enable refinement if CCD produces a zero ToI"   ON)

include(FetchContent)
FetchContent_Declare(
    ccd_wrapper
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/CCD-Wrapper.git
    GIT_TAG 23907dadf3e1eef606e38450ada4aa4f96fd9f71
)
FetchContent_MakeAvailable(ccd_wrapper)
