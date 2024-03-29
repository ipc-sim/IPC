################################################################################
# See comments and discussions here:
# http://stackoverflow.com/questions/5088460/flags-to-enable-thorough-and-verbose-g-warnings
################################################################################

if(TARGET IPC::warnings)
  return()
endif()

set(IPC_FLAGS
  # -Wall
  # -Wextra
  # -pedantic

  # -Wconversion
  #-Wunsafe-loop-optimizations # broken with C++11 loops
  # -Wunused

  # -Wno-long-long
  # -Wpointer-arith
  # -Wformat=2
  # -Wuninitialized
  # -Wcast-qual
  # -Wmissing-noreturn
  # -Wmissing-format-attribute
  # -Wredundant-decls

  # -Werror=implicit
  # -Werror=nonnull
  # -Werror=init-self
  # -Werror=main
  # -Werror=missing-braces
  # -Werror=sequence-point
  # -Werror=return-type
  # -Werror=trigraphs
  # -Werror=array-bounds
  # -Werror=write-strings
  # -Werror=address
  # -Werror=int-to-pointer-cast
  # -Werror=pointer-to-int-cast

  # -Wno-unused-variable
  # -Wunused-but-set-variable
  # -Wno-unused-parameter

  #-Weffc++
  # -Wno-old-style-cast
  # -Wno-sign-conversion
  #-Wsign-conversion

  # -Wshadow

  # -Wstrict-null-sentinel
  # -Woverloaded-virtual
  # -Wsign-promo
  # -Wstack-protector
  # -Wstrict-aliasing
  # -Wstrict-aliasing=2

  # Warn whenever a switch statement has an index of enumerated type and
  # lacks a case for one or more of the named codes of that enumeration.
  # -Wswitch
  # This is annoying if all cases are already covered.
  # -Wswitch-default
  # This is annoying if there is a default that covers the rest.
  # -Wswitch-enum
  # -Wswitch-unreachable
  # -Wcovered-switch-default # Annoying warnings from nlohmann::json

  # -Wcast-align
  # -Wdisabled-optimization
  #-Winline # produces warning on default implicit destructor
  # -Winvalid-pch
  # -Wmissing-include-dirs
  # -Wpacked
  # -Wno-padded
  # -Wstrict-overflow
  # -Wstrict-overflow=2

  # -Wctor-dtor-privacy
  # -Wlogical-op
  # -Wnoexcept
  # -Woverloaded-virtual
  # -Wundef

  # -Wnon-virtual-dtor
  # -Wdelete-non-virtual-dtor
  # -Werror=non-virtual-dtor
  # -Werror=delete-non-virtual-dtor

  # -Wno-sign-compare

  ###########
  # GCC 6.1 #
  ###########

  # -Wnull-dereference
  # -fdelete-null-pointer-checks
  # -Wduplicated-cond
  # -Wmisleading-indentation

  #-Weverything

  ###########################
  # Enabled by -Weverything #
  ###########################

  #-Wdocumentation
  #-Wdocumentation-unknown-command
  #-Wfloat-equal

  #-Wglobal-constructors
  #-Wexit-time-destructors
  #-Wmissing-variable-declarations
  #-Wextra-semi
  #-Wweak-vtables
  #-Wno-source-uses-openmp
  #-Wdeprecated
  #-Wnewline-eof
  #-Wmissing-prototypes

  #-Wno-c++98-compat
  #-Wno-c++98-compat-pedantic

  ###########################
  # Need to check if those are still valid today
  ###########################

  #-Wimplicit-atomic-properties
  #-Wmissing-declarations
  #-Wmissing-prototypes
  #-Wstrict-selector-match
  #-Wundeclared-selector
  #-Wunreachable-code

  # Not a warning, but enable link-time-optimization
  # TODO: Check out modern CMake version of setting this flag
  # https://cmake.org/cmake/help/latest/module/CheckIPOSupported.html
  #-flto

  # Gives meaningful stack traces
  # -fno-omit-frame-pointer
  # -fno-optimize-sibling-calls

  ####################
  # Disabled warning #
  ####################
  -Wno-unused-command-line-argument
)

# Flags above don't make sense for MSVC
if(MSVC)
  set(IPC_FLAGS)
endif()

include(CheckCXXCompilerFlag)

add_library(IPC_warnings INTERFACE)
add_library(IPC::warnings ALIAS IPC_warnings)

foreach(FLAG IN ITEMS ${IPC_FLAGS})
  string(REPLACE "=" "-" FLAG_VAR "${FLAG}")
  if(NOT DEFINED IS_SUPPORTED_${FLAG_VAR})
    check_cxx_compiler_flag("${FLAG}" IS_SUPPORTED_${FLAG_VAR})
  endif()
  if(IS_SUPPORTED_${FLAG_VAR})
    target_compile_options(IPC_warnings INTERFACE ${FLAG})
  endif()
endforeach()
