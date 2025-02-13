cmake_minimum_required(VERSION 3.9.0) # 3.9 for IPO
project(SeqFile
	VERSION 0.1.0
	DESCRIPTION "An efficient fasta/fastq processing library"
	LANGUAGES C)

# Default to a release build.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "No build type selected; defaulting to Release")
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()

set(CMAKE_C_STANDARD 11)

include(CheckIncludeFile)
CHECK_INCLUDE_FILE(threads.h HAVE_C11_THREADS)
if(NOT HAVE_C11_THREADS)
	set(THREAD_LIB tinycthread)
else()
	find_package(Threads REQUIRED)
	set(THREAD_LIB Threads::Threads)
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_is_supported OUTPUT ipo_error)
if(ipo_is_supported)
	message(STATUS "IPO / LTO enabled")
else()
	message(STATUS "IPO / LTO not supported: <${ipo_error}>")
endif()

# Set cmake variables 
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH "Installation directory for libraries")
set(INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Installation directory for headers")

# SeqFile user configuration files
option(SEQF_BUILD_SHARED "Build shared SeqFile library" ON)
option(SEQF_BUILD_STATIC "Build static SeqFile library" ON)
option(SEQF_BUILD_TESTS "Build SeqFile unit tests" ON)
option(SEQF_SKIP_INSTALL_ALL "Don't install any targets" OFF)
option(SEQF_SKIP_INSTALL_LIBRARIES "Don't install shared and static libraries" OFF)
option(SEQF_SKIP_INSTALL_STATIC "Don't install static library" OFF)
option(SEQF_SKIP_INSTALL_SHARED "Don't install shared library" OFF)
option(SEQF_SKIP_INSTALL_HEADER "Don't install header files" OFF)

# Find compression library
message(CHECK_START "Finding ISA-L")
find_library(COMPRESSION_LIB
	NAMES isal
	PATHS "/usr/lib"
	      "/usr/lib64"
	      "/usr/local/lib"
	      "/usr/local/lib64"
	      "$ENV{HOME}/.local/lib"
	      "$ENV{HOME}/.local/lib64")
find_path(COMPRESSION_INC
	NAMES igzip_lib.h
	PATHS "$ENV{HOME}/.local/include"
	      "/usr/include"
	      "/usr/local/include")
if(COMPRESSION_LIB AND COMPRESSION_INC)
	set(SEQF_HAS_ISA_L TRUE)
	message(CHECK_PASS "found: " ${COMPRESSION_LIB})
else()
	set(SEQF_HAS_ISA_L FALSE)
	message(CHECK_FAIL "not found, attempting ZLIB...")
	find_package(ZLIB REQUIRED)
	set(COMPRESSION_LIB ZLIB::ZLIB)
	set(COMPRESSION_INC "")
endif()

set(CMAKE_C_FLAGS_DEBUG "-O0 -ggdb3")

# Begin building project
add_subdirectory(include)
add_subdirectory(source)
if(SEQF_BUILD_TESTS)
	add_subdirectory(test)
endif(SEQF_BUILD_TESTS)

# Install the library 
if(NOT SEQF_SKIP_INSTALL_LIBRARIES AND NOT SEQF_SKIP_INSTALL_ALL)
	if(NOT SEQF_SKIP_INSTALL_STATIC AND SEQF_BUILD_STATIC)
		install(TARGETS seqf_static LIBRARY DESTINATION ${INSTALL_LIB_DIR})
	endif()
	if(NOT SEQF_SKIP_INSTALL_SHARED AND SEQF_BUILD_SHARED)
		install(TARGETS seqf_shared LIBRARY DESTINATION ${INSTALL_LIB_DIR})
	endif()
endif()

if(NOT SEQF_SKIP_INSTALL_HEADER AND NOT SEQF_SKIP_INSTALL_ALL)
	install(FILES ${SEQF_PUBLIC_HEADERS} DESTINATION ${INSTALL_INC_DIR})
endif()
