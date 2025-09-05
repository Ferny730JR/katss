set(STRUCTURE_SOURCE_FILES
	"bpp_tables.c"
	"structure.c")

add_library(katss_structure STATIC ${STRUCTURE_SOURCE_FILES})
target_include_directories(katss_structure PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(katss_structure PRIVATE 
	KATSS_MEMORYUTILS
	KATSS_STRING
	seqf_static
	${THREAD_LIB})
target_compile_definitions(katss_structure PRIVATE
	${C11_THREADS_DEFINE})

# Find ViennaRNA library
find_library(LIBRARY_VRNA 
	NAMES RNA
	HINTS "/usr/lib"
		  "/usr/lib64"
		  "/usr/local/lib"
		  "/usr/local/lib64"
		  "$ENV{HOME}/.local/lib"
		  "$ENV{HOME}/.local/lib64")
find_path(INCLUDE_DIR_VRNA
	NAMES ViennaRNA/fold.h
	PATHS "$ENV{HOME}/.local/include"
		  "/usr/include"
		  "/usr/local/include")

if(LIBRARY_VRNA AND INCLUDE_DIR_VRNA)
	message(STATUS "Found ViennaRNA: " ${LIBRARY_VRNA})
	target_link_libraries(katss_structure PRIVATE ${LIBRARY_VRNA})
	target_include_directories(katss_structure PRIVATE ${INCLUDE_DIR_VRNA})
else()
	message(SEND_ERROR "Couldn't find libRNA")
endif()

include(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(m log2 "" HAVE_LIB_M)
if(HAVE_LIB_M AND NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	set(EXTRA_LIBS ${EXTRA_LIBS} m)
endif()
target_link_libraries(katss_structure PRIVATE ${EXTRA_LIBS})

if(ipo_is_supported)
	set_property(TARGET katss_structure PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()
