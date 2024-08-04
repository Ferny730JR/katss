# Don't install and build tests for tinycthread
set(TINYCTHREAD_DISABLE_TESTS ON CACHE BOOL "Disable TinyCThread unit tests")
set(TINYCTHREAD_INSTALL OFF CACHE BOOL "Install a static library for TinyCThread")

# Build dependencies first, then katss, lastly binaries
if(NOT HAVE_C11_THREADS)
	add_subdirectory(tinycthread)
endif()
add_subdirectory(katss)
add_subdirectory(bin)
