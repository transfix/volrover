# VolumeRover CMakeLists.txt - Qt6 Version
# Based on original CMakeLists.txt, updated for Qt6
# $Id: CMakeLists.txt 1405 2010-02-26 21:16:13Z transfix $

cmake_minimum_required(VERSION 3.16)
project(VolumeRover VERSION 2.1.0 LANGUAGES C CXX)

# On Windows, ensure boost::asio's WinSock2.h is included before any
# transitive inclusion of windows.h.  Without these defines, MSVC errors
# with "WinSock.h has already been included" once boost::asio is pulled in
# (e.g. via boost::process or boost::interprocess from libcvc).
if(WIN32)
  add_compile_definitions(WIN32_LEAN_AND_MEAN _WIN32_WINNT=0x0601 NOMINMAX)
endif()

# Set C++ standard.  libcvc 3.1.0 requires C++20 (uses std::format,
# concepts, etc.) and is consumed via FetchContent below.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# AppleClang treats narrowing conversions as errors in C++11+.
# This legacy codebase has many benign narrowing conversions, so
# downgrade to warnings rather than fixing hundreds of call sites.
if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang")
  add_compile_options(-Wno-error=c++11-narrowing -Wno-c++11-narrowing)
endif()

# Enable automatic MOC, UIC, RCC for Qt
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Global build variables
option(USING_MAKE_ENV_VARS "If this is enabled, use the CPPFLAGS and LDFLAGS environment variables" OFF)
if(USING_MAKE_ENV_VARS)
  string(COMPARE NOTEQUAL "$ENV{CPPFLAGS}" "" CPPFLAGS)
  string(COMPARE NOTEQUAL "$ENV{LDFLAGS}" "" LDFLAGS)
  if(CPPFLAGS AND LDFLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} $ENV{CPPFLAGS}")
    set(CMAKE_LD_FLAGS "$ENV{LDFLAGS}")
  endif()
endif()

option(USING_MAKE_ENV_VARS_FOR_INSTALL "If this is enabled, use the LDFLAGS environment variable as a directory to look up libs for a bundle" ON)
mark_as_advanced(USING_MAKE_ENV_VARS_FOR_INSTALL)
if(USING_MAKE_ENV_VARS_FOR_INSTALL)
  set(ENV_LIB_DIR "$ENV{ENV}/lib")
endif()

# Set the macOS deployment target (update to modern version)
if(APPLE)
  option(USING_MACOSX_DEV_TARGET "If this is enabled, you need to specify the mac osx target for build" OFF)
  if(USING_MACOSX_DEV_TARGET)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures for macOS")
  endif()
endif()

if(NOT BUNDLE_LIB_SEARCH_DIRS)
  set(BUNDLE_LIB_SEARCH_DIRS "${ENV_LIB_DIR}")
endif()
message(STATUS "BUNDLE_LIB_SEARCH_DIRS: ${BUNDLE_LIB_SEARCH_DIRS}")

# CMake Modules
set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_SOURCE_DIR}/CMake/cuda
  ${CMAKE_MODULE_PATH})

# Build output 
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# Set the revision number from git (updated from SVN)
option(USING_GIT_VERSION "Use git for version info" ON)
set(CVC_GIT_REVISION "0")
if(USING_GIT_VERSION)
  find_package(Git)
  if(Git_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE CVC_GIT_REVISION
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    message(STATUS "Git revision: ${CVC_GIT_REVISION}")
  endif()
endif()

# The VolumeRover version number
set(VolumeRover_VERSION_MAJOR 2)
set(VolumeRover_VERSION_MINOR 1)
set(VolumeRover_VERSION_PATCH 0)

# For CI/dev builds, append git short hash to distinguish artifacts.
# Tagged release builds (when VOLROVER_RELEASE is ON) use the clean version.
option(VOLROVER_RELEASE "Build a release without git hash in version string" OFF)
if(NOT VOLROVER_RELEASE AND CVC_GIT_REVISION AND NOT CVC_GIT_REVISION STREQUAL "0")
  set(VolumeRover_VERSION ${VolumeRover_VERSION_MAJOR}.${VolumeRover_VERSION_MINOR}.${VolumeRover_VERSION_PATCH}.${CVC_GIT_REVISION})
else()
  set(VolumeRover_VERSION ${VolumeRover_VERSION_MAJOR}.${VolumeRover_VERSION_MINOR}.${VolumeRover_VERSION_PATCH})
endif()
message(STATUS "VolumeRover version: ${VolumeRover_VERSION}")

# macOS Bundle settings
if(APPLE)
  set(CMAKE_BUNDLE_NAME
    "VolumeRover-${VolumeRover_VERSION_MAJOR}.${VolumeRover_VERSION_MINOR}.${VolumeRover_VERSION_PATCH}")
  set(CMAKE_BUNDLE_LOCATION "${CMAKE_INSTALL_PREFIX}")
  # make sure CMAKE_INSTALL_PREFIX ends in /
  string(LENGTH "${CMAKE_INSTALL_PREFIX}" LEN)
  math(EXPR LEN "${LEN} -1")
  string(SUBSTRING "${CMAKE_INSTALL_PREFIX}" ${LEN} 1 ENDCH)
  if(NOT "${ENDCH}" STREQUAL "/")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/")
  endif()
  set(CMAKE_INSTALL_PREFIX 
    "${CMAKE_INSTALL_PREFIX}${CMAKE_BUNDLE_NAME}.app/Contents")
endif()

# Include macros for setting up library dependencies
# Use SetupQt6.cmake for Qt6 migration
option(USE_QT6 "Use Qt6 instead of Qt4" ON)
if(USE_QT6)
  include(SetupQt6)
else()
  include(SetupQt)
endif()

include(SetupFFTW)
include(SetupGSL)

# Setup GLEW (required for OpenGL extension support)
find_package(GLEW REQUIRED)
if(GLEW_FOUND)
  message(STATUS "Found GLEW ${GLEW_VERSION}")
  message(STATUS "GLEW includes: ${GLEW_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "GLEW not found. Install libglew-dev")
endif()

# Setup Boost (required)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
# Note: system is header-only since Boost 1.69 and should not be listed as a compiled COMPONENT
find_package(Boost 1.70.0 REQUIRED COMPONENTS thread date_time regex filesystem program_options)
if(Boost_FOUND)
  message(STATUS "Found Boost ${Boost_VERSION}")
  message(STATUS "Boost includes: ${Boost_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "Boost not found. Install libboost-all-dev or set BOOST_ROOT")
endif()

# log4cplus logging library (system package)
# Uses FindLog4cplus.cmake: tries config mode first, falls back to manual detection
find_package(Log4cplus REQUIRED)
message(STATUS "Found log4cplus")

# pthreads (needed by XmlRPC/VolumeViewer multi-tile server code)
if(WIN32)
  find_package(PThreads4W REQUIRED)
else()
  find_package(Threads REQUIRED)
endif()

option(DISABLE_CGAL "This disables CGAL even if it is available on the build system." OFF)
if(NOT DISABLE_CGAL)
  include(SetupCGAL)
endif()

# Possibly prevent crashing on macOS (still relevant)
if(APPLE)
  add_definitions(-D_GLIBCXX_FULLY_DYNAMIC_STRING)
endif()

# Prevent windows.h min/max macros from clobbering std::min/std::max
if(WIN32)
  add_definitions(-DNOMINMAX)
endif()

# Suppress libcvc's PascalCase compat typedefs (App, ThreadInfo, etc.).
# volrover provides its own subclasses in inc/cvc_compat.h that auto-supply
# the volrover-owned cvc::app context to the new (app& ctx, ...) ctors.
add_definitions(-DCVC_COMPAT_APP_DEFINED)

# ************* CPack configuration for VolumeRover *************
# Include the packaging system for distributing binaries

if(WIN32)
  set(CPACK_GENERATOR "NSIS;ZIP")
elseif(APPLE)
  set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)
  set(CPACK_GENERATOR "DragNDrop;ZIP")
else()
  set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)
  set(CPACK_GENERATOR "TGZ;DEB")
endif()

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Volumetric image visualization and processing")
set(CPACK_PACKAGE_VENDOR "Computational Visualization Center")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_VERSION_MAJOR "${VolumeRover_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${VolumeRover_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${VolumeRover_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION "${VolumeRover_VERSION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "VolumeRover ${VolumeRover_VERSION_MAJOR}.${VolumeRover_VERSION_MINOR}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "VolumeRover-${VolumeRover_VERSION}")
set(CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")

# Installers for 32- vs. 64-bit CMake
if(CMAKE_CL_64)
  set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")
  set(CPACK_NSIS_PACKAGE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} (Win64)")
  set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION} (Win64)")
else()
  set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
  set(CPACK_NSIS_PACKAGE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
  set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} VolumeRover")
  set(CPACK_NSIS_HELP_LINK "https://github.com/transfix/volrover")
  set(CPACK_NSIS_URL_INFO_ABOUT "https://github.com/transfix/volrover")
  set(CPACK_NSIS_CONTACT "transfix@ices.utexas.edu")
  set(CPACK_NSIS_MODIFY_PATH ON)
  set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}")
endif()

if(NOT DEFINED CPACK_SYSTEM_NAME)
  if("${CMAKE_SYSTEM_NAME}" STREQUAL "CYGWIN")
    set(CPACK_SYSTEM_NAME Cygwin)
  else()
    set(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  endif()
endif()

if(${CPACK_SYSTEM_NAME} MATCHES Windows)
  if(CMAKE_CL_64)
    set(CPACK_SYSTEM_NAME win64-x64)
  else()
    set(CPACK_SYSTEM_NAME win32-x86)
  endif()
endif()

if(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
  if(CYGWIN)
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_SOURCE_PACKAGE_FILE_NAME}")
  else()
    set(CPACK_PACKAGE_FILE_NAME 
      "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")
  endif()
endif()

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Joe Rivera")
set(CPACK_PACKAGE_CONTACT "transfix@ices.utexas.edu")
set(CPACK_STRIP_FILES ON)
set(CPACK_SOURCE_STRIP_FILES "")
set(CPACK_PACKAGE_EXECUTABLES "VolumeRover2" "VolumeRover")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/README.md")

# Install README into the package
install(FILES README.md DESTINATION .)

# Debian package dependencies
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libqt6core6, libqt6gui6, libqt6widgets6, libqt6opengl6, libboost-filesystem1.74.0, libboost-system1.74.0")

# dpkg-shlibdeps can't find bundled libraries (libcvc, libQGLViewer, gtest,
# gmock) in the build tree; disable automatic shlib dependency generation
# and rely on the explicit CPACK_DEBIAN_PACKAGE_DEPENDS list above.
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS OFF)

include(CPack)
include(InstallRequiredSystemLibraries)

# Installation location option
option(USING_STANDARD_INSTALL_LOCATION "If on, it will not put everything under $prefix/VolumeRover2" ON)
if(USING_STANDARD_INSTALL_LOCATION)
  set(ADDED_PREFIX ".")
else()
  set(ADDED_PREFIX "VolumeRover2")
endif()

# ************* libcvc dependency *************
# Resolution order, falling through on failure:
#   1. find_package(cvc CONFIG) against a pre-installed libcvc on
#      CMAKE_PREFIX_PATH (fastest, no download or build).
#   2. Download the prebuilt release tarball matching this platform /
#      build type from https://github.com/transfix/libcvc/releases and
#      find_package against the extracted tree.
#   3. FetchContent_Declare(GIT_REPOSITORY ...) and build from source.
#      Override with -DLIBCVC_FORCE_SOURCE=ON or
#      -DFETCHCONTENT_SOURCE_DIR_LIBCVC=/path/to/libcvc for local dev.
#
# volrover links cvc::xmlrpc (used by VolumeGridRover, CompServ, all the
# Segmentation/Seg* modules, SecStruct, GenSeg, VolumeRover2). libcvc
# release archives <= 3.1.1 do not ship cvc::xmlrpc, so the resolver
# verifies the target exists before accepting a prebuilt package and
# falls through to the source build otherwise. Once libcvc starts
# shipping xmlrpc, the prebuilt path activates automatically.
set(LIBCVC_TAG     "v3.2.1" CACHE STRING "libcvc release tag to consume")
set(LIBCVC_VERSION "3.2.1"  CACHE STRING "libcvc package version to look up")
set(LIBCVC_GIT_SHA "63a7306785359818459ce19eb1343c16a8298e49"
    CACHE STRING "libcvc commit SHA used when falling through to source build")
option(LIBCVC_FORCE_SOURCE "Skip prebuilt libcvc and always build from source" OFF)
set(LIBCVC_NEEDS_XMLRPC TRUE)

include(FetchContent)

# Disable libcvc sub-features we don't need during volrover build
set(CVC_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(CVC_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
set(CVC_ENABLE_COVERAGE OFF CACHE BOOL "" FORCE)
set(CVC_BUILD_VOLROVER3 OFF CACHE BOOL "" FORCE)
set(CVC_USING_XMLRPC ON CACHE BOOL "" FORCE)
set(CVC_USING_IMOD_MRC OFF CACHE BOOL "" FORCE)

function(_libcvc_check_usable out_var)
  set(${out_var} FALSE PARENT_SCOPE)
  if(NOT cvc_FOUND)
    return()
  endif()
  if(LIBCVC_NEEDS_XMLRPC AND NOT TARGET cvc::xmlrpc)
    message(STATUS "libcvc package found (${cvc_VERSION}) but cvc::xmlrpc is missing")
    return()
  endif()
  set(${out_var} TRUE PARENT_SCOPE)
endfunction()

set(_libcvc_usable      FALSE)
set(_libcvc_source_built FALSE)

# 1. Pre-installed package
if(NOT LIBCVC_FORCE_SOURCE)
  find_package(cvc ${LIBCVC_VERSION} CONFIG QUIET)
  _libcvc_check_usable(_libcvc_usable)
  if(_libcvc_usable)
    message(STATUS "Using pre-installed libcvc ${cvc_VERSION}")
  endif()
endif()

# 2. Prebuilt release tarball
if(NOT _libcvc_usable AND NOT LIBCVC_FORCE_SOURCE)
  if(WIN32)
    set(_libcvc_plat "windows-x86_64")
    set(_libcvc_ext  "zip")
  elseif(APPLE)
    set(_libcvc_plat "macos-arm64")
    set(_libcvc_ext  "zip")
  elseif(UNIX)
    set(_libcvc_plat "linux-x86_64")
    set(_libcvc_ext  "tar.gz")
  endif()
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(_libcvc_cfg "debug")
  else()
    set(_libcvc_cfg "release")
  endif()
  if(DEFINED _libcvc_plat)
    set(_libcvc_archive "libcvc-${LIBCVC_VERSION}-${_libcvc_plat}-${_libcvc_cfg}.${_libcvc_ext}")
    set(_libcvc_url     "https://github.com/transfix/libcvc/releases/download/${LIBCVC_TAG}/${_libcvc_archive}")
    set(_libcvc_stage   "${CMAKE_BINARY_DIR}/libcvc-prebuilt")
    set(_libcvc_cfgfile "${_libcvc_stage}/lib/cmake/cvc/cvcConfig.cmake")
    if(NOT EXISTS "${_libcvc_cfgfile}")
      message(STATUS "Downloading prebuilt libcvc: ${_libcvc_url}")
      file(MAKE_DIRECTORY "${_libcvc_stage}")
      file(DOWNLOAD "${_libcvc_url}" "${_libcvc_stage}/${_libcvc_archive}"
           STATUS _libcvc_dl TLS_VERIFY ON)
      list(GET _libcvc_dl 0 _libcvc_dl_code)
      if(_libcvc_dl_code EQUAL 0)
        execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf "${_libcvc_archive}"
          WORKING_DIRECTORY "${_libcvc_stage}"
          RESULT_VARIABLE _libcvc_extract)
        if(NOT _libcvc_extract EQUAL 0)
          message(STATUS "Failed to extract ${_libcvc_archive} (exit ${_libcvc_extract})")
        endif()
      else()
        message(STATUS "Prebuilt libcvc download failed: ${_libcvc_dl}")
      endif()
    endif()
    if(EXISTS "${_libcvc_cfgfile}")
      list(PREPEND CMAKE_PREFIX_PATH "${_libcvc_stage}")
      find_package(cvc ${LIBCVC_VERSION} CONFIG QUIET)
      _libcvc_check_usable(_libcvc_usable)
      if(_libcvc_usable)
        message(STATUS "Using prebuilt libcvc ${cvc_VERSION} from ${_libcvc_url}")
      endif()
    endif()
  endif()
endif()

# 3. Source build (fallback)
if(NOT _libcvc_usable)
  message(STATUS "Building libcvc ${LIBCVC_TAG} from source")
  FetchContent_Declare(
    libcvc
    GIT_REPOSITORY https://github.com/transfix/libcvc.git
    GIT_TAG        ${LIBCVC_GIT_SHA}
  )
  FetchContent_MakeAvailable(libcvc)
  set(_libcvc_source_built TRUE)
endif()

# On case-insensitive filesystems (macOS, Windows), libcvc's in-tree
# headers can be shadowed by volrover's own include directories named
# similarly. Only relevant for the source-build path; the imported
# cvc::cvc target from find_package carries an absolute include path
# that does not collide.
if(_libcvc_source_built)
  include_directories(BEFORE ${libcvc_SOURCE_DIR}/inc)
  include_directories(BEFORE ${libcvc_BINARY_DIR}/inc)
endif()

add_subdirectory(src)

# ************* Unit tests *************
option(VOLROVER_BUILD_TESTS "Build volrover unit tests" ON)
if(VOLROVER_BUILD_TESTS)
  enable_testing()

  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  add_subdirectory(tests)
endif()
