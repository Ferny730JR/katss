# Set source files for library
set(KATSS_SOURCE_FILES
	"${CMAKE_CURRENT_SOURCE_DIR}/hash_functions.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/tables.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/seqseq.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/counter.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/uncounter.c"
	"${CMAKE_CURRENT_SOURCE_DIR}/enrichments.c")

# Create objects for source files
add_library(OBJLIB OBJECT "${KATSS_SOURCE_FILES}")

# Shared libraries need PIC
#set_property(TARGET OBJLIB PROPERTY POSITION_INDEPENDENT_CODE 1)
target_compile_definitions(OBJLIB PRIVATE KATSS_VERBOSE=$<BOOL:${KKCTR_VERBOSE}>)

# Set optimization flags for release builds
if(ipo_is_supported)
	set_property(TARGET OBJLIB PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
	target_compile_options(OBJLIB PRIVATE "-O3")
endif()

# Specify include directories
target_link_libraries(OBJLIB PRIVATE tinycthread KATSS_MEMORYUTILS KATSS_RNAFILES)
target_include_directories(OBJLIB PRIVATE ${KKCTR_INCLUDE_DIR})

# Create static library
add_library(kkctr_static STATIC $<TARGET_OBJECTS:OBJLIB>)

if(ipo_is_supported)
	set_property(TARGET kkctr_static PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()

# Link libraries with katss
include(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(m log2 "" HAVE_LIB_M)
if(HAVE_LIB_M AND NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	set(EXTRA_LIBS ${EXTRA_LIBS} m)
endif()

#target_link_libraries(kkctr_shared # RNAFILES links tinycthread
#	ZLIB::ZLIB KATSS_MEMORYUTILS KATSS_RNAFILES ${EXTRA_LIBS})
target_link_libraries(kkctr_static KATSS_MEMORYUTILS KATSS_RNAFILES ${EXTRA_LIBS})

# Set public include directory
#target_include_directories(kkctr_shared PUBLIC ${KKCTR_INCLUDE_DIR})
target_include_directories(kkctr_static PUBLIC ${KKCTR_INCLUDE_DIR})

# Set library properties
#set_target_properties(kkctr_shared PROPERTIES
#	VERSION ${PROJECT_VERSION}
#	OUTPUT_NAME kkctr)
set_target_properties(kkctr_static PROPERTIES
	VERSION ${PROJECT_VERSION}
	OUTPUT_NAME kkctr)
