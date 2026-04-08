cmake_minimum_required(VERSION 3.16)

# Setup Qt6 or Qt4
if(USE_QT6)
  find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets Xml OpenGL OpenGLWidgets)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(QT_LIBRARIES Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Xml Qt6::OpenGL Qt6::OpenGLWidgets)
else()
  SetupQt()
endif()

# Ensure large file support
add_definitions(-D_FILE_OFFSET_BITS=64)

set(SOURCE_FILES
  VolumeGridRover.cpp
  PointClassFile.cpp
  ContourFile.cpp
  SDF2D.cpp
#  bspline_opt.cpp
#  bspline_fit.cpp
  sdf_opt.cpp
)

set(INCLUDE_FILES
  ../../inc/VolumeGridRover/VolumeGridRover.h
  ../../inc/VolumeGridRover/SurfRecon.h
  ../../inc/VolumeGridRover/PointClassFile.h
  ../../inc/VolumeGridRover/ContourFile.h
  ../../inc/VolumeGridRover/SDF2D.h
#  ../../inc/VolumeGridRover/bspline_opt.h
  ../../inc/VolumeGridRover/sdf_opt.h
)

file(GLOB UI_FORMS *.ui)

source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})
source_group("UI Files" FILES ${UI_FORMS})

include_directories(
  ../../inc
  ${CMAKE_BINARY_DIR}/src/VolumeGridRover  #for generated headers
)

OPTION(VOLUMEGRIDROVER_USING_MEDAX "Compile VolumeGridRover using the medial axis code" ON)
MARK_AS_ADVANCED(VOLUMEGRIDROVER_USING_MEDAX)
IF(VOLUMEGRIDROVER_USING_MEDAX)
  SET(INCLUDE_FILES
    ${INCLUDE_FILES}
    ../../inc/VolumeGridRover/medax.h
  )
  ADD_DEFINITIONS(-DUSING_VOLUMEGRIDROVER_MEDAX)

  OPTION(VOLUMEGRIDROVER_MEDAX_INSERT_EDGES "Use the 'insert edges' option of the medial axis code" ON)
  mark_as_advanced(VOLUMEGRIDROVER_MEDAX_INSERT_EDGES)
  IF(VOLUMEGRIDROVER_MEDAX_INSERT_EDGES)
    ADD_DEFINITIONS(-DMEDAX_INSERT_EDGES)
  ENDIF(VOLUMEGRIDROVER_MEDAX_INSERT_EDGES)

  OPTION(VOLUMEGRIDROVER_MEDAX_ALWAYS_INSERT "Use the 'always insert' option of the medial axis code" ON)
  mark_as_advanced(VOLUMEGRIDROVER_MEDAX_ALWAYS_INSERT)
  IF(VOLUMEGRIDROVER_MEDAX_ALWAYS_INSERT)
    ADD_DEFINITIONS(-DMEDAX_ALWAYS_INSERT)
  ENDIF(VOLUMEGRIDROVER_MEDAX_ALWAYS_INSERT)
ENDIF(VOLUMEGRIDROVER_USING_MEDAX)

OPTION(VOLUMEGRIDROVER_EM_CLUSTERING "Compile VolumeGridRover using Sangmin's EM clustering - causes some serious memory problems if enabled!!" OFF)
MARK_AS_ADVANCED(VOLUMEGRIDROVER_EM_CLUSTERING)
IF(VOLUMEGRIDROVER_EM_CLUSTERING)
  #Note: EM.cpp is kept around for reference.  It isn't actually used.

  SET(INCLUDE_FILES
    ${INCLUDE_FILES}
    ../../inc/VolumeGridRover/EM.h
  )

  ADD_DEFINITIONS(-DUSING_EM_CLUSTERING)
ENDIF(VOLUMEGRIDROVER_EM_CLUSTERING)

set(MOC_HEADERS
  ../../inc/VolumeGridRover/VolumeGridRover.h
  ../../inc/VolumeGridRover/bspline_opt.h
  ../../inc/VolumeGridRover/sdf_opt.h
)

# Create the library (Qt6 uses AUTOMOC/AUTOUIC, Qt4 uses manual wrapping)
if(USE_QT6)
  add_library(VolumeGridRover STATIC
    ${SOURCE_FILES}
    ${INCLUDE_FILES}
    ${UI_FORMS}
  )
elseif(QT4_FOUND)
  qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})
  qt4_wrap_ui(UI_FILES ${UI_FORMS})

  add_library(VolumeGridRover STATIC
    ${SOURCE_FILES}
    ${INCLUDE_FILES}
    ${MOC_SOURCES}
    ${UI_FILES}
  )
else()
  message(FATAL_ERROR "VolumeGridRover requires Qt4 or Qt6.")
endif()

SetupCGAL(VolumeGridRover)
SetupGSL(VolumeGridRover)

if(NOT BUILD_SEGMENTATION_LIB)
  message(SEND_ERROR "VolumeGridRover requires the Virus segmentation library!")
endif(NOT BUILD_SEGMENTATION_LIB)

SET(LINK_LIBS
  VolMagick
  ColorTable2
  GenSeg
  XmlRPC
  QGLViewer
)

#If we're using Tiling, the VolumeGridRover needs to link
#with Tiling and GeometryFileTypes.
# 11/3/11 edwardsj - converting to new contour tiler
IF(BUILD_TILING_LIB)
  SET(LINK_LIBS
    ${LINK_LIBS}
    ContourFilterLib
    ContourTiler
    ContourFilterLib
#    Tiling
#    GeometryFileTypes
  )
ENDIF(BUILD_TILING_LIB)

#If we're enabling 2D isocontouring, we need to link with
#libContour
option(VOLUMEGRIDROVER_ISOCONTOURING "Enable 2D isocontouring via the VolumeGridRover" ON)
mark_as_advanced(VOLUMEGRIDROVER_ISOCONTOURING)
if(VOLUMEGRIDROVER_ISOCONTOURING)
  add_definitions(-DVOLUMEGRIDROVER_ISOCONTOURING)  
  set(LINK_LIBS
    ${LINK_LIBS}
    Contour
  )
endif(VOLUMEGRIDROVER_ISOCONTOURING)

target_link_libraries(VolumeGridRover
  PRIVATE
    Boost::thread
    Boost::date_time
    Boost::regex
    Boost::filesystem
    ${LINK_LIBS}
)

option(BUILD_VOLUMEGRIDROVER_STANDALONE "Build standalone VolumeGridRover app." OFF)
if(BUILD_VOLUMEGRIDROVER_STANDALONE)
  set(STANDALONE_INCLUDE_FILES
    ../../inc/VolumeGridRover/VolumeGridRoverMainWindow.h
  )
  set(STANDALONE_SOURCE_FILES
    VolumeGridRoverMainWindow.cpp
    main.cpp
  )
  set(STANDALONE_UI_FORMS
    VolumeGridRoverMainWindowBase.ui
  )

  source_group("Source Files" FILES ${STANDALONE_SOURCE_FILES})
  source_group("UI Files" FILES ${STANDALONE_UI_FORMS})
  source_group("Header Files" FILES ${STANDALONE_INCLUDE_FILES})

  # Qt6 uses AUTOMOC/AUTOUIC, Qt4 uses manual wrapping
  if(USE_QT6)
    add_executable(GridRover 
       ${STANDALONE_INCLUDE_FILES}
       ${STANDALONE_SOURCE_FILES}
       ${STANDALONE_UI_FORMS}
    )
  elseif(QT4_FOUND)
    set(STANDALONE_MOC_HEADERS
      ../../inc/VolumeGridRover/VolumeGridRoverMainWindow.h
    )
    qt4_wrap_cpp(STANDALONE_MOC_SOURCES ${STANDALONE_MOC_HEADERS})

    add_executable(GridRover 
       ${STANDALONE_INCLUDE_FILES}
       ${STANDALONE_SOURCE_FILES}
       ${STANDALONE_MOC_SOURCES}
       ${STANDALONE_UI_FORMS}
    )
  else()
    message(FATAL_ERROR "GridRover requires Qt4 or Qt6.")
  endif()

  target_link_libraries(GridRover VolumeGridRover ${QT_LIBRARIES})
endif(BUILD_VOLUMEGRIDROVER_STANDALONE)
