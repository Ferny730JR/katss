cmake_minimum_required(VERSION 3.9) # for IPO / LTO
project(KatssKmerCounter
	VERSION 0.9.0
	DESCRIPTION "K-mer counting function utilities."
	LANGUAGES C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED TRUE)

# Set cmake variables 
set(KKCTR_INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" 
	CACHE PATH "Installation directory for libraries")
set(KKCTR_INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include/KATSS" 
	CACHE PATH "Installation directory for headers")

option(KKCTR_VERBOSE "Output error/warning messages to stderr" OFF)
option(KKCTR_SKIP_INSTALL_ALL "Don't install any targets" OFF)
option(KKCTR_SKIP_INSTALL_LIBRARIES "Don't install shared and static libraries" OFF)
option(KKCTR_SKIP_INSTALL_STATIC "Don't install static library" OFF)
option(KKCTR_SKIP_INSTALL_SHARED "Don't install shared library" OFF)
option(KKCTR_SKIP_INSTALL_HEADER "Don't install header files" OFF)

# Set public header/API for library
set(KKCTR_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(KKCTR_PUBLIC_HEADERS
	"${KKCTR_INCLUDE_DIR}/counter.h"
	"${KKCTR_INCLUDE_DIR}/hash_functions.h"
	"${KKCTR_INCLUDE_DIR}/enrichments.h"
	"${KKCTR_INCLUDE_DIR}/seqseq.h")

add_subdirectory(source)

# Install the library 
if(NOT KKCTR_SKIP_INSTALL_LIBRARIES AND NOT KKCTR_SKIP_INSTALL_ALL)
	if(NOT KKCTR_SKIP_INSTALL_STATIC)
		install(TARGETS kkctr_static LIBRARY DESTINATION ${KKCTR_INSTALL_LIB_DIR})
	endif()
	if(NOT KKCTR_SKIP_INSTALL_SHARED)
		install(TARGETS kkctr_shared LIBRARY DESTINATION ${KKCTR_INSTALL_LIB_DIR})
	endif()
endif()

if(NOT KKCTR_SKIP_INSTALL_HEADER AND NOT KKCTR_SKIP_INSTALL_ALL)
	install(FILES ${KKCTR_PUBLIC_HEADERS} DESTINATION ${KKCTR_INSTALL_INC_DIR})
endif()
