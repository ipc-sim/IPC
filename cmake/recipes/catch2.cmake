# Catch2 (https://github.com/catchorg/Catch2)
# License: BSL-1.0
if(TARGET Catch2::Catch2)
    return()
endif()

message(STATUS "Third-party: creating target 'Catch2::Catch2'")

option(CATCH_CONFIG_CPP17_STRING_VIEW "Enable support for std::string_view" ON)
option(CATCH_INSTALL_DOCS "Install documentation alongside library" OFF)
option(CATCH_INSTALL_EXTRAS "Install extras alongside library" OFF)

include(CPM)
CPMAddPackage("gh:catchorg/Catch2@2.13.6")