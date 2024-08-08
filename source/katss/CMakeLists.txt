# Don't install KmerCounter 
set(KKCTR_SKIP_INSTALL_ALL ON CACHE BOOL "Don't install any targets")
set(KKCTR_VERBOSE ON CACHE BOOL "Output error/warning messages to stderr")

# build modules
add_subdirectory(helpers)
add_subdirectory(KmerCounter)
add_subdirectory(Structure)
