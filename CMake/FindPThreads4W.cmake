# FindPThreads4W.cmake — locate the pthreads4w (pthreads-win32) library
# and provide the PThreads4W::PThreads4W imported target.
#
# vcpkg's `pthreads` port ships a PThreads4WConfig.cmake, which is where
# this target historically came from (ci.yml / release.yml build with
# the vcpkg toolchain). The cvcpkg `pthreads4w` bundle ships only the
# headers and libraries (pthread.h, _ptw32.h, pthreadVC*.lib/.dll), so
# this module bridges the gap:
#
#   1. Try CONFIG mode first — if a PThreads4WConfig.cmake is on the
#      prefix path (vcpkg), defer to it. This keeps the vcpkg-based CI
#      workflows byte-for-byte identical in behavior.
#   2. Otherwise find the headers/library on CMAKE_PREFIX_PATH (cvcpkg
#      prefix) and synthesize the same imported target.

find_package(PThreads4W CONFIG QUIET)
if(TARGET PThreads4W::PThreads4W)
  set(PThreads4W_FOUND TRUE)
  return()
endif()

find_path(PTHREADS4W_INCLUDE_DIR
  NAMES pthread.h
  PATH_SUFFIXES include)

# MSVC C-cleanup ABI ("VC") variants, version 3 then 2.
find_library(PTHREADS4W_LIBRARY
  NAMES pthreadVC3 pthreadVC2
  PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PThreads4W
  REQUIRED_VARS PTHREADS4W_LIBRARY PTHREADS4W_INCLUDE_DIR)

if(PThreads4W_FOUND AND NOT TARGET PThreads4W::PThreads4W)
  add_library(PThreads4W::PThreads4W UNKNOWN IMPORTED)
  set_target_properties(PThreads4W::PThreads4W PROPERTIES
    IMPORTED_LOCATION "${PTHREADS4W_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PTHREADS4W_INCLUDE_DIR}")
endif()

mark_as_advanced(PTHREADS4W_INCLUDE_DIR PTHREADS4W_LIBRARY)
