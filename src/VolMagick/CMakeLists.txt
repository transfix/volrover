# joe - 07/23/2011 - adding cvc
cmake_minimum_required(VERSION 3.16)

if(WIN32)
  add_definitions(-D__WINDOWS__)
endif()

set(INCLUDE_FILES 
  ../../inc/VolMagick/BoundingBox.h
  ../../inc/VolMagick/Dimension.h
  ../../inc/VolMagick/Exceptions.h
  ../../inc/VolMagick/VolMagick.h
  ../../inc/VolMagick/Voxels.h
  ../../inc/VolMagick/Volume.h
  ../../inc/VolMagick/VolumeFileInfo.h
  ../../inc/VolMagick/VolumeFile_IO.h
  ../../inc/VolMagick/Utility.h
  ../../inc/VolMagick/endians.h
  ../../inc/VolMagick/VolumeCache.h
  ../../inc/VolMagick/Null_IO.h
  ../../inc/VolMagick/RawIV_IO.h
  ../../inc/VolMagick/RawV_IO.h
  ../../inc/VolMagick/MRC_IO.h
  ../../inc/VolMagick/Spider_IO.h
  ../../inc/VolMagick/VTK_IO.h
  ../../inc/VolMagick/VoxelOperationStatusMessenger.h
  ../../inc/VolMagick/StdErrOpStatus.h
)

SET(SOURCE_FILES
  ContrastEnhancement.cpp
  VTK_IO.cpp
  GDTVFilter.cpp
  Null_IO.cpp
  BilateralFilter.cpp
  MRC_IO.cpp
  VolumeCache.cpp
  AnisotropicDiffusion.cpp
  Voxels.cpp
  Volume.cpp
  VolumeFileInfo.cpp
  VolumeFile_IO.cpp
  Utility.cpp
  RawV_IO.cpp
  Spider_IO.cpp
  RawIV_IO.cpp
  VoxelOperationStatusMessenger.cpp
  StdErrOpStatus.cpp
)

# the core CVC lib that contains CVC::App
SET(VOLMAGICK_LINK_LIBS CVC)

# The VolMagick version number.
SET(VOLMAGICK_VERSION_MAJOR 1)
SET(VOLMAGICK_VERSION_MINOR 0)
SET(VOLMAGICK_VERSION_PATCH ${CVC_SVN_REVISION})

# The full version string
SET(VOLMAGICK_VERSION ${VOLMAGICK_VERSION_MAJOR}.${VOLMAGICK_VERSION_MINOR}.${VOLMAGICK_VERSION_PATCH})
message(STATUS "VolMagick version: ${VOLMAGICK_VERSION}")

add_definitions(-DVOLMAGICK_VERSION_STRING=\"${VOLMAGICK_VERSION}\")

if(CVC_USING_HDF5)
  find_package(HDF5 COMPONENTS C CXX)
  if(HDF5_FOUND)
    message("HDF5 found! Enabling CVC file support")

    add_definitions(-DVOLMAGICK_USING_HDF5 ${HDF5_DEFINITIONS})
    set(INCLUDE_FILES
       ${INCLUDE_FILES}
       ../../inc/VolMagick/HDF5_IO.h
    )
    set(SOURCE_FILES
       ${SOURCE_FILES}
       HDF5_IO.cpp
    )

    # Resolve the HDF5 link libs. The upstream/cvcpkg HDF5 CONFIG package
    # (hdf5-config.cmake) sets HDF5_FOUND and defines imported targets, but leaves
    # HDF5_LIBRARIES EMPTY — so HDF5_IO.cpp's H5:: C++ API went unlinked (undefined
    # H5::DataSpace at final link). When empty, derive from the config's imported
    # targets: one C++ (hdf5_cpp) + one C (hdf5), preferring shared, probing
    # cvcpkg's non-namespaced names AND the namespaced forms (vcpkg/brew). The _cpp
    # target transitively pulls the C lib; we list both for a deterministic line.
    # The system FindHDF5-module path (apt/brew) sets HDF5_LIBRARIES and is used
    # as-is.
    set(VOLMAGICK_HDF5_LIBS ${HDF5_LIBRARIES})
    if(NOT VOLMAGICK_HDF5_LIBS)
      foreach(_cpp hdf5_cpp-shared hdf5::hdf5_cpp-shared hdf5_cpp-static hdf5::hdf5_cpp-static hdf5_cpp hdf5::hdf5_cpp)
        if(TARGET ${_cpp})
          list(APPEND VOLMAGICK_HDF5_LIBS ${_cpp})
          break()
        endif()
      endforeach()
      foreach(_c hdf5-shared hdf5::hdf5-shared hdf5-static hdf5::hdf5-static hdf5 hdf5::hdf5)
        if(TARGET ${_c})
          list(APPEND VOLMAGICK_HDF5_LIBS ${_c})
          break()
        endif()
      endforeach()
    endif()
    if(NOT VOLMAGICK_HDF5_LIBS) # last resort: raw component vars from FindHDF5 module
      set(VOLMAGICK_HDF5_LIBS ${HDF5_CXX_LIBRARIES} ${HDF5_C_LIBRARIES})
    endif()
    SET(VOLMAGICK_LINK_LIBS ${VOLMAGICK_LINK_LIBS} ${VOLMAGICK_HDF5_LIBS})

    # HDF5 include dirs for target_include_directories. CRITICAL: do NOT mix the
    # Debian/Ubuntu /usr/include/hdf5/serial headers with a deps/cvcpkg HDF5. Their
    # hsize_t differs — serial is `unsigned long long`, the deps build is
    # `uint64_t` (== unsigned long) — so if the system serial path shadows the deps
    # headers, HDF5_IO.cpp emits H5:: methods taking `unsigned long long*` that the
    # linked libhdf5_cpp (which exports `unsigned long*`) does not provide →
    # undefined references at the final link. Prefer find_package's include dirs;
    # if it gave none, pull them off the imported target; only fall back to the
    # system serial path for a bare apt build where neither is available.
    set(VOLMAGICK_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS} ${HDF5_C_INCLUDE_DIRS} ${HDF5_CXX_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR})
    if(NOT VOLMAGICK_HDF5_INCLUDE_DIRS)
      foreach(_t hdf5_cpp-shared hdf5::hdf5_cpp-shared hdf5-shared hdf5::hdf5-shared hdf5_cpp-static hdf5-static)
        if(TARGET ${_t})
          get_target_property(_inc ${_t} INTERFACE_INCLUDE_DIRECTORIES)
          if(_inc)
            list(APPEND VOLMAGICK_HDF5_INCLUDE_DIRS ${_inc})
            break()
          endif()
        endif()
      endforeach()
    endif()
    if(NOT VOLMAGICK_HDF5_INCLUDE_DIRS AND EXISTS "/usr/include/hdf5/serial")
      list(APPEND VOLMAGICK_HDF5_INCLUDE_DIRS "/usr/include/hdf5/serial")
    endif()
    if(VOLMAGICK_HDF5_INCLUDE_DIRS)
      list(REMOVE_DUPLICATES VOLMAGICK_HDF5_INCLUDE_DIRS)
    endif()
  else(HDF5_FOUND)
    message("HDF5 not found! CVC file support disabled")
  endif(HDF5_FOUND)
endif(CVC_USING_HDF5)

option(VOLMAGICK_USING_IMOD_MRC "Add support for using imod's MRC loading routines." OFF)
mark_as_advanced(VOLMAGICK_USING_IMOD_MRC)
if(VOLMAGICK_USING_IMOD_MRC)
  add_definitions(-DVOLMAGICK_USING_IMOD_MRC -DNOTIFFLIBS)
  file(GLOB IMOD_MRC_LIBIIMOD_HEADERS ../../inc/VolMagick/libiimod/*.h)
  set(INCLUDE_FILES
     ${INCLUDE_FILES}
     ${IMOD_MRC_LIBIIMOD_HEADERS}
     ../../inc/VolMagick/IMOD_MRC_IO.h
  )
  set(SOURCE_FILES
     ${SOURCE_FILES}
     IMOD_MRC_IO.cpp
     libiimod/b3dutil.c
     libiimod/diffusion.c
     libiimod/iilikemrc.c
     libiimod/iimage.c
     libiimod/iimrc.c
     libiimod/iitif.c
     libiimod/ilist.c
     libiimod/mrcfiles.c
     libiimod/mrcsec.c
     libiimod/mrcslice.c
     libiimod/plist.c
     libiimod/sliceproc.c
     libiimod/tiffstub.c
     libiimod/islice.c
     libiimod/parallelwrite.c
  )
endif(VOLMAGICK_USING_IMOD_MRC)

option(VOLMAGICK_USING_VOLMAGICK_INR "Add support for INR format.  Little used." OFF)
mark_as_advanced(VOLMAGICK_USING_VOLMAGICK_INR)
if(VOLMAGICK_USING_VOLMAGICK_INR)
  add_definitions(-DVOLMAGICK_USING_VOLMAGICK_INR)
  set(INCLUDE_FILES
     ${INCLUDE_FILES}
     ../../inc/VolMagick/INR_IO.h
  )
  set(SOURCE_FILES
     ${SOURCE_FILES}
     INR_IO.cpp
  ) 
endif(VOLMAGICK_USING_VOLMAGICK_INR)

# QtImage support disabled - Qt3 no longer supported
# option(VOLMAGICK_USING_QTIMAGE "Add support for loading volumes as images using Qt" OFF)
# mark_as_advanced(VOLMAGICK_USING_QTIMAGE)
if(FALSE) # Disabled - requires Qt3 which is no longer supported
  add_definitions(-DVOLMAGICK_USING_QTIMAGE)
  set(INCLUDE_FILES
    ${INCLUDE_FILES}
    ../../inc/VolMagick/QtImage_IO.h
  )
  set(SOURCE_FILES
    ${SOURCE_FILES}
    QtImage_IO.cpp
  )
  set(VOLMAGICK_LINK_LIBS
    ${VOLMAGICK_LINK_LIBS}
    ${QT_LIBRARIES}
  )
endif()

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${INCLUDE_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
)

ADD_LIBRARY (VolMagick STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

# Add HDF5 include directories if enabled
if(CVC_USING_HDF5 AND HDF5_FOUND)
  target_include_directories(VolMagick PUBLIC ${VOLMAGICK_HDF5_INCLUDE_DIRS})
endif()

set_target_properties(VolMagick PROPERTIES
  AUTOMOC OFF
  AUTOUIC OFF
  AUTORCC OFF
)

target_link_libraries(VolMagick
  PUBLIC
    Boost::thread
    Boost::date_time
    Boost::regex
    Boost::filesystem
    ${VOLMAGICK_LINK_LIBS}
)
