add_executable(ikke ikke.c ikke_cmdl.c)
target_link_libraries(ikke PUBLIC
	kkctr_static
	KATSS_STRING)
if(ipo_is_supported)
    set_property(TARGET ikke PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
    target_compile_options(ikke "-O3")
endif()

add_executable(kstruct kstruct.c kstruct_cmdl.c)
target_link_libraries(kstruct PUBLIC
	$<$<BOOL:${OpenMP_FOUND}>:OpenMP::OpenMP_C>
	katss_structure
	KATSS_MEMORYUTILS
	KATSS_STRING
	stdc++)
if(ipo_is_supported)
    set_property(TARGET kstruct PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
    target_compile_options(kstruct "-O3")
endif()
