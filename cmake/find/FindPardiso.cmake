# Find Pardiso library
# --------------------
#
# Defines the following variables:
#   PARDISO_LIBRARIES    Path to the Pardiso libraries to link with
#   PARDISO_INCLUDE_DIR  Path to the Pardiso include directory
#
################################################################################

set(PARDISO_SEARCH_PATHS
		${PARDISO_INSTALL_PREFIX}
		"$ENV{PARDISO_INSTALL_PREFIX}"
		"/usr/local/"
		"$ENV{PROGRAMFILES}/PARDISO"
		"$ENV{PROGRAMW6432}/PARDISO"
		"$ENV{HOME}/.local/"
		"$ENV{HOME}/.pardiso/"
		"/usr/local/lib/"
)

find_library(PARDISO_LIBRARIES
		pardiso500-GNU461-X86-64
		pardiso500-GNU472-X86-64
		pardiso500-GNU481-X86-64
		pardiso500-INTEL1301-X86-64
		pardiso500-MACOS-X86-64
		pardiso500-MPI-GNU450-X86-64
		pardiso500-MPI-GNU461-X86-64
		pardiso500-MPI-GNU463-X86-64
		pardiso500-MPI-GNU472-X86-64
		pardiso500-MPI-INTEL1301-X86-64
		pardiso500-WIN-X86-64
		pardiso600-GNU720-X86-64
		pardiso600-GNU800-X86-64
		pardiso600-MACOS-X86-64
	PATHS ${PARDISO_SEARCH_PATHS}
	PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARDISO DEFAULT_MSG PARDISO_LIBRARIES)

mark_as_advanced(PARDISO_LIBRARIES)
