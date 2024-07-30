cmake_minimum_required(VERSION 3.9.0) # for IPO
project(KATSS
	VERSION 0.9.0
	DESCRIPTION "K-mer counting function utilities."
	LANGUAGES C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED TRUE)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_is_supported OUTPUT ipo_error)
if(ipo_is_supported)
	message(STATUS "IPO / LTO enabled")
else()
	message(STATUS "IPO / LTO not supported: <${ipo_error}>")
endif()

# Set cmake variables 
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH "Installation directory for libraries")
set(INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include/KATSS" CACHE PATH "Installation directory for headers")

option(VERBOSE "Output error/warning messages to stderr" ON)
option(SKIP_INSTALL_ALL "Don't install any targets" OFF)
option(SKIP_INSTALL_LIBRARIES "Don't install shared and static libraries" OFF)
option(SKIP_INSTALL_STATIC "Don't install static library" OFF)
option(SKIP_INSTALL_SHARED "Don't install shared library" OFF)
option(SKIP_INSTALL_HEADER "Don't install header files" OFF)

# Find zlib library
find_package(ZLIB REQUIRED)

add_subdirectory(source)
