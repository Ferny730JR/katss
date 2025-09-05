cmake_minimum_required(VERSION 3.10.0)
project(KATSS
	VERSION 0.9.0
	DESCRIPTION "K-mer counting function utilities."
	LANGUAGES C)

# Default to a release build.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "No build type selected; defaulting to Release")
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()

set(CMAKE_C_STANDARD 11)

# Check if C11 threads is supported, otherwise find any thread functionality
include(CheckIncludeFile)
CHECK_INCLUDE_FILE(threads.h HAVE_C11_THREADS)
if(NOT HAVE_C11_THREADS)
	if(WIN32)
		set(C11_THREADS_DEFINE "__STDC_NO_THREADS__=1")
	endif()
	set(THREAD_LIB tinycthread)
else()
	find_package(Threads REQUIRED)
	set(THREAD_LIB Threads::Threads)
endif()

# Find if OpenMP is supported
find_package(OpenMP)

# Check for IPO optimizations
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

set(CMAKE_C_FLAGS_DEBUG "-O0 -ggdb3")
add_subdirectory(source)

# Install rules
if(NOT SKIP_INSTALL_ALL)
	if(NOT SKIP_INSTALL_BINARIES)
		install(TARGETS ikke
		        RUNTIME DESTINATION "${INSTALL_BIN_DIR}")
	endif()

	if(NOT SKIP_INSTALL_LIBRARIES)
		set_target_properties(kkctr_static PROPERTIES
			VERSION ${PROJECT_VERSION}
			OUTPUT_NAME katss)
		install(TARGETS kkctr_static
			LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
			ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")
	endif()

	if(NOT SKIP_INSTALL_HEADER)
		install(FILES ${KKCTR_PUBLIC_HEADERS}
			DESTINATION "${INSTALL_INC_DIR}/katss")
	endif()
endif()
