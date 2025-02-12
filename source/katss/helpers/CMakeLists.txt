# Memory utilities library
add_library(KATSS_MEMORYUTILS OBJECT "memory_utils.c")
target_include_directories(KATSS_MEMORYUTILS PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(KATSS_MEMORYUTILS PRIVATE KATSS_VERBOSE=$<BOOL:${VERBOSE}>)

# KATSS string library
add_library(KATSS_STRING OBJECT "string_utils.c")
target_include_directories(KATSS_STRING PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Set optimizations
if(ipo_is_supported)
    set_property(TARGET KATSS_MEMORYUTILS PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set_property(TARGET KATSS_STRING      PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
    target_compile_options(KATSS_MEMORYUTILS PRIVATE "-O3")
    target_compile_options(KATSS_STRING      PRIVATE "-O3")
endif()
