# Set source files for library
set(KATSS_SOURCE_FILES
	"${CMAKE_CURRENT_SOURCE_DIR}/hash_functions.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/tables.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/seqseq.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/counter.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/recounter.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/uncounter.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/enrichments.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/ushuffle.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/katss_helpers.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/katss_count.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/katss_enrichment.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/katss_ikke.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/thread_safe_rand.c"
	)

# Check for math library, which is used by katss
include(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(m log2 "" HAVE_LIB_M)
if(HAVE_LIB_M AND NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	set(EXTRA_LIBS ${EXTRA_LIBS} m)
endif()

# Create katss kmer counting static library
add_library(kkctr_static STATIC ${KATSS_SOURCE_FILES})

target_include_directories(kkctr_static PUBLIC
	${KKCTR_INCLUDE_DIR})

target_link_libraries(kkctr_static PRIVATE
	${THREAD_LIB}
	KATSS_MEMORYUTILS
	T_TEST_LIB
	seqf_static
	${EXTRA_LIBS})

target_compile_definitions(kkctr_static PRIVATE
	KATSS_VERBOSE=$<BOOL:${KKCTR_VERBOSE}>
	${C11_THREADS_DEFINE})

if(ipo_is_supported)
	set_property(TARGET kkctr_static PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
	target_compile_options(kkctr_static PRIVATE "-O3")
endif(ipo_is_supported)

set_target_properties(kkctr_static PROPERTIES
	VERSION ${PROJECT_VERSION}
	OUTPUT_NAME kkctr)
