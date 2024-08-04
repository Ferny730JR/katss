# Memory utilities library
add_library(KATSS_MEMORYUTILS OBJECT "memory_utils.c")
target_include_directories(KATSS_MEMORYUTILS PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(KATSS_MEMORYUTILS PRIVATE KATSS_VERBOSE=$<BOOL:${VERBOSE}>)

# RNA file parsing file
add_library(KATSS_RNAFILES OBJECT "rnafiles.c")
target_include_directories(KATSS_RNAFILES PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(KATSS_RNAFILES PRIVATE ${COMPRESSION_INC})
target_link_libraries(KATSS_RNAFILES PRIVATE ${THREAD_LIB} ${COMPRESSION_LIB})
target_compile_definitions(KATSS_RNAFILES PRIVATE _HAS_ISA_L_=$<BOOL:${KATSS_HAS_ISA_L}>)

# KATSS string library
add_library(KATSS_STRING OBJECT "string_utils.c")
target_include_directories(KATSS_STRING PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Set optimizations
if(ipo_is_supported)
    set_property(TARGET KATSS_MEMORYUTILS PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set_property(TARGET KATSS_RNAFILES    PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set_property(TARGET KATSS_STRING      PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
    target_compile_options(KATSS_MEMORYUTILS PRIVATE "-O3")
    target_compile_options(KATSS_RNAFILES    PRIVATE "-O3")
    target_compile_options(KATSS_STRING      PRIVATE "-O3")
endif()
