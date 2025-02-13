# Don't install and build tests for tinycthread
set(TINYCTHREAD_DISABLE_TESTS ON CACHE BOOL "Disable TinyCThread unit tests")
set(TINYCTHREAD_INSTALL OFF CACHE BOOL "Install a static library for TinyCThread")

# Don't install or compile shared library for SeqFile
set(SEQF_SKIP_INSTALL_ALL ON CACHE BOOL "Don't install seqf")
set(SEQF_BUILD_SHARED OFF CACHE BOOL "Disable seqf shared library")
set(SEQF_BUILD_TESTS OFF CACHE BOOL "Don't build seqfile unit tests")

# Build dependencies first, then katss, lastly binaries
if(NOT HAVE_C11_THREADS)
	add_subdirectory(tinycthread)
endif()
add_subdirectory(seqfile)
add_subdirectory(katss)
add_subdirectory(bin)
