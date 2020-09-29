include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
  set(IPC_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
  set(IPC_EXTRA_OPTIONS "")
endif()

function(custom_download_project name)
  download_project(
    PROJ         ${name}
    SOURCE_DIR   ${IPC_EXTERNAL}/${name}
    DOWNLOAD_DIR ${IPC_EXTERNAL}/.cache/${name}
    QUIET
    ${IPC_EXTRA_OPTIONS}
    ${ARGN}
  )
endfunction()

################################################################################

# OSQP
function(download_osqp)
  custom_download_project(osqp
    GIT_REPOSITORY https://github.com/oxfordcontrol/osqp.git
    # GIT_TAG        657a3b117320c4f8ecb57e27011e75e63f61ce4d
    GIT_TAG        v0.4.1
  )
endfunction()

# libigl
function(download_libigl)
  custom_download_project(libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG        efee81b7dbc81ec87adaca1197b47f4faab961d3
  )
endfunction()

# TBB
function(download_tbb)
  custom_download_project(tbb
    # GIT_REPOSITORY https://github.com/intel/tbb.git
    # GIT_TAG        2018_U5
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG        344fa84f34089681732a54f5def93a30a3056ab9
  )
endfunction()

# Logger
function(download_spdlog)
    custom_download_project(spdlog
       GIT_REPOSITORY https://github.com/gabime/spdlog.git
       GIT_TAG        v1.4.2
    )
endfunction()

# AMGCL
function(download_amgcl)
    custom_download_project(amgcl
       GIT_REPOSITORY https://github.com/ddemidov/amgcl.git
       GIT_TAG        461a66ce6d197a3816218bf94ffd114a367c1ef1
    )
endfunction()

# Catch2
function(download_catch2)
    custom_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.10.2
    )
endfunction()

# finite-diff
function(download_finite_diff)
    custom_download_project(finite-diff
        GIT_REPOSITORY https://github.com/zfergus/finite-diff.git
        GIT_TAG        dac4e629defb0edef050e56050eefb00574ecbaa
    )
endfunction()

# CLI11
function(download_cli11)
    custom_download_project(cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG        b979d3a37098da7ad6291d9bd3c4fdec3a705043
    )
endfunction()

# Eigen Gurobi Wrapper
function(download_eigen_gurobi)
  custom_download_project(eigen-gurobi
    GIT_REPOSITORY https://github.com/zfergus/eigen-gurobi.git
    GIT_TAG        51b1aacb3c5733555d09fe362887d618ee97826d
  )
endfunction()

# CCD Wrapper (includes Rational CCD)
function(download_ccd_wrapper)
  custom_download_project(ccd-wrapper
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/CCD-Wrapper.git
    GIT_TAG        ae980163f2507f2a678e08bfeb36f774536ea958
  )
endfunction()
