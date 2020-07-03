set(GIT_PATH "NOTSET" CACHE STRING "Absolute path to the git executable")
set(PYTHON_NAME "python" CACHE STRING "Name of the python executable, eg. python2")
set(PYTHON_PATH "NOTSET" CACHE STRING "Absolute path to the python executable")
set(CLANGFORMAT_PATH "NOTSET" CACHE STRING "Absolute path to the clang-format executable")
set(CLANGFORMAT_STYLE "file" CACHE STRING "Paremter pass to clang-format -style=<here>")
set(IGNORE_LIST "" CACHE STRING "Semi colon separeted list of directories to ignore")

if("${GIT_PATH}" STREQUAL "NOTSET")
	find_program(FIND_GIT git)
	if("${FIND_GIT}" STREQUAL "FIND_GIT-NOTFOUND")
		message(FATAL_ERROR "Could not find 'git' please set GIT_PATH:STRING")
	else()
		set(GIT_PATH ${FIND_GIT})
		message(STATUS "Found: ${GIT_PATH}")
	endif()
else()
	if(NOT EXISTS ${GIT_PATH})
		message(WARNING "Could not find git: ${GIT_PATH}")
	else()
		message(STATUS "Found: ${GIT_PATH}")
	endif()
endif()

if("${PYTHON_PATH}" STREQUAL "NOTSET")
	find_program(FIND_PYTHON ${PYTHON_NAME})
	if("${FIND_PYTHON}" STREQUAL "FIND_PYTHON-NOTFOUND")
		message(FATAL_ERROR "Could not find '${PYTHON_NAME}' please set PYTHON_PATH:STRING")
	else()
		set(PYTHON_PATH ${FIND_PYTHON})
		message(STATUS "Found: ${PYTHON_PATH}")
	endif()
else()
	if(NOT EXISTS ${PYTHON_PATH})
		message(WARNING "Could not find python: ${PYTHON_PATH}")
	else()
		message(STATUS "Found: ${PYTHON_PATH}")
	endif()
endif()

if("${CLANGFORMAT_PATH}" STREQUAL "NOTSET")
	find_program(FIND_CLANGFORMAT clang-format)
	if("${FIND_CLANGFORMAT}" STREQUAL "FIND_CLANGFORMAT-NOTFOUND")
		message(FATAL_ERROR "Could not find 'clang-format' please set CLANGFORMAT_PATH:STRING")
	else()
		set(CLANGFORMAT_PATH ${FIND_CLANGFORMAT})
		message(STATUS "Found: ${CLANGFORMAT_PATH}")
	endif()
else()
	if(NOT EXISTS ${CLANGFORMAT_PATH})
		message(WARNING "Could not find clang-format: ${CLANGFORMAT_PATH}")
	else()
		message(STATUS "Found: ${CLANGFORMAT_PATH}")
	endif()
endif()

set(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/git-cmake-format.py)
string(REGEX REPLACE "\\\\" "/" GIT_PATH ${GIT_PATH})
string(REGEX REPLACE "\\\\" "/" PYTHON_PATH ${PYTHON_PATH})
string(REGEX REPLACE "\\\\" "/" CLANGFORMAT_PATH ${CLANGFORMAT_PATH})

execute_process(COMMAND git rev-parse --show-toplevel
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_ROOT
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT GIT_ROOT)
  message(WARNING "Not in a git repository")
else()
  message(STATUS "Found git root: ${GIT_ROOT}")
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/git-pre-commit-hook
    ${GIT_ROOT}/.git/hooks/pre-commit)
  if(EXISTS ${GIT_ROOT}/Working)
    execute_process(COMMAND git rev-parse --show-toplevel
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Working
      OUTPUT_VARIABLE GIT_WORKING_ROOT
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/git-pre-commit-hook
      ${GIT_ROOT}/Working/.git/hooks/pre-commit)
  endif()
endif()

add_custom_target(format ALL
	${PYTHON_PATH} ${SCRIPT} --cmake ${GIT_PATH} ${CLANGFORMAT_PATH} -style=${CLANGFORMAT_STYLE} -ignore=${IGNORE_LIST}
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
