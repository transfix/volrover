cmake_minimum_required(VERSION 3.16)

SetupQt()

# Prevent Qt's foreach macro from interfering with Boost.Foreach namespace
add_definitions(-DQT_NO_FOREACH)

if(ENABLE_STEREO_DISPLAY)
  add_definitions(-DENABLE_STEREO_DISPLAY)
endif()

set(ISOCONTOURING_LIB "")
option(VOLUMEROVER2_ISOCONTOURING_WITH_LBIE "Build VolumeRover2 with LBIE isocontouring support (can call FastContouring depending on method selected)" ON)
mark_as_advanced(VOLUMEROVER2_ISOCONTOURING_WITH_LBIE)
if(VOLUMEROVER2_ISOCONTOURING_WITH_LBIE)
  add_definitions(-DISOCONTOURING_WITH_LBIE)
  set(ISOCONTOURING_LIB LBIE)
endif(VOLUMEROVER2_ISOCONTOURING_WITH_LBIE)

option(VOLUMEROVER2_ISOCONTOURING_WITH_FASTCONTOURING "Build VolumeRover2 with direct FastContouring isocontouring support" ON)
mark_as_advanced(VOLUMEROVER2_ISOCONTOURING_WITH_FASTCONTOURING)
if(VOLUMEROVER2_ISOCONTOURING_WITH_FASTCONTOURING)
  add_definitions(-DISOCONTOURING_WITH_FASTCONTOURING)
  set(ISOCONTOURING_LIB FastContouring)
endif(VOLUMEROVER2_ISOCONTOURING_WITH_FASTCONTOURING)

# Enabling both at the same time isn't a problem now.  In fact, it's preferred!
# This is only an issue for using the deprecated VolumeViewerPage
#if(VOLUMEROVER2_ISOCONTOURING_WITH_LBIE AND VOLUMEROVER2_ISOCONTOURING_WITH_FASTCONTOURING)
#  message(SEND_ERROR "VolumeRover2: Cannot choose both LBIE and FastContouring based isocontouring!")
#endif(VOLUMEROVER2_ISOCONTOURING_WITH_LBIE AND VOLUMEROVER2_ISOCONTOURING_WITH_FASTCONTOURING)

option(VOLUMEROVER2_THREADED_ISOCONTOURING "Run interactive isocontouring routines in a separate thread from the main GUI" ON)
mark_as_advanced(VOLUMEROVER2_THREADED_ISOCONTOURING)
if(VOLUMEROVER2_THREADED_ISOCONTOURING)
  add_definitions(-DTHREADED_ISOCONTOURING)
endif(VOLUMEROVER2_THREADED_ISOCONTOURING)

add_definitions(-DVOLUMEROVER_VERSION_STRING=\"${VolumeRover_VERSION}\")

option(VOLUMEROVER2_FLOAT_COLORTABLE "Use float values for uploaded color table for greater precision.  Only works with CG and ARB volume rendering paths for now" OFF)
mark_as_advanced(VOLUMEROVER2_FLOAT_COLORTABLE)
if(VOLUMEROVER2_FLOAT_COLORTABLE)
  add_definitions(-DVOLUMEROVER2_FLOAT_COLORTABLE)
endif(VOLUMEROVER2_FLOAT_COLORTABLE)

#FILE(GLOB SOURCE_FILES *.cpp)
set(SOURCE_FILES
#  VolumeViewerPage.cpp
  main.cpp
  ${CMAKE_SOURCE_DIR}/src/volrover_app_instance.cpp
  CVCMainWindow.cpp
  Viewers.cpp
  VolumeViewer.cpp
  BoundingBoxModify.cpp
  ImportData.cpp
#  VolumeRoverMain.cpp
  DimensionModify.cpp
  GeometryInterface.cpp
  ContoursInterface.cpp
  VolumeMemoryInterface.cpp
  NewVolumeDialog.cpp
  VolumeInterface.cpp
  RemapVoxels.cpp
  ReconstructionDialog.cpp
  AnisotropicDiffusionDialog.cpp
  BilateralFilterDialog.cpp
  ContourTilerDialog.cpp
  ContrastEnhancementDialog.cpp
  GDTVFilterDialog.cpp
  PocketTunnelDialog.cpp
  CurationDialog.cpp
  TightCoconeDialog.cpp
  SuperSecondaryStructuresDialog.cpp
  SkeletonizationDialog.cpp
  LBIEQualityImprovementDialog.cpp
  MultiTileServerDialog.cpp
  VolumeViewer2.cpp
  Application.cpp
)
#FILE(GLOB INCLUDE_FILES ../../inc/VolumeRover2/*.[Hh])
set(INCLUDE_FILES
  ../../inc/VolumeRover2/DimensionModify.h
  ../../inc/VolumeRover2/VolMagickEventHandler.h
  ../../inc/VolumeRover2/DataWidget.h
#  ../../inc/VolumeRover2/VolumeViewerPage.h
  ../../inc/VolumeRover2/NewVolumeDialog.h
#  ../../inc/VolumeRover2/VolumeRoverMain.h
  ../../inc/VolumeRover2/ImportData.h
  ../../inc/VolumeRover2/BoundingBoxModify.h
  ../../inc/VolumeRover2/CVCMainWindow.h
  ../../inc/VolumeRover2/RemapVoxels.h
  ../../inc/VolumeRover2/VolumeViewer.h
  ../../inc/VolumeRover2/VolumeMemoryInterface.h
  ../../inc/VolumeRover2/Viewers.h
  ../../inc/VolumeRover2/VolumeInterface.h
  ../../inc/VolumeRover2/VolumeRoverExtension.h
  ../../inc/VolumeRover2/MultiTileServerDialog.h
  ../../inc/VolumeRover2/VolumeViewer2.h
  ../../inc/VolumeRover2/GeometryInterface.h
  ../../inc/VolumeRover2/ContoursInterface.h
  ../../inc/VolumeRover2/ReconstructionDialog.h
  ../../inc/VolumeRover2/AnisotropicDiffusionDialog.h
  ../../inc/VolumeRover2/BilateralFilterDialog.h
  ../../inc/VolumeRover2/ContourTilerDialog.h
  ../../inc/VolumeRover2/ContrastEnhancementDialog.h
  ../../inc/VolumeRover2/GDTVFilterDialog.h
  ../../inc/VolumeRover2/PocketTunnelDialog.h
  ../../inc/VolumeRover2/CurationDialog.h
  ../../inc/VolumeRover2/TightCoconeDialog.h
  ../../inc/VolumeRover2/SuperSecondaryStructuresDialog.h
  ../../inc/VolumeRover2/SkeletonizationDialog.h
  ../../inc/VolumeRover2/LBIEQualityImprovementDialog.h
)

# Qt6-compatible UI forms (Qt4 .ui files work with Qt6 AUTOUIC)
set(UI_FORMS
  HLSSurfaceDialog.ui
  AnisotropicDiffusionDialog.ui
  RemapVoxels.ui
  SecondaryStructureDialog.ui
  HistogramDialog.ui
  CurationDialog.ui
  ContrastEnhancementDialog.ui
  NewVolumeDialog.ui
  ImportData.ui
  TightCoconeDialog.ui
  SuperSecondaryStructuresDialog.ui
  DimensionModify.ui
  SegmentVirusMapDialog.ui
  UnknownData.ui
  BoundingBoxModify.ui
  BilateralFilterDialog.ui
  ContourTilerDialog.ui
  PocketTunnelDialog.ui
  LBIE_dialog.ui
  LBIE_qualityImprovement.ui
  ReconstructionDialog.ui
  CVCMainWidget.ui
  VolumeInterface.ui
  GeometryInterface.ui
  ContoursInterface.ui
  VolumeMemoryInterface.ui
  GDTVFilterDialog.ui
  SkeletonizationDialog.ui
  VolumeViewerPage.ui
  VolumeViewerPageManipulators.ui
  AddPropertyDialog.ui
  MainOptions.ui
  SaveFileDialog.ui
  SliceRenderingDialog.ui
  CloseDataDialog.ui
  SaveImageDialog.ui
  ViewGeometryOptionsDialog.ui
  SelectCurrentVolumeDialog.ui
  MultiTileServerDialog.ui
  generateRawvDialog.ui
  generateMeshDialog.ui
  mslevelsetDialog.ui
  mplevelsetDialog.ui
)

source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})

include_directories(
  ../../inc
  ${CMAKE_BINARY_DIR}/src/VolumeRover2  #for generated headers
)

set(MOC_HEADERS
#  ../../inc/VolumeRover2/VolumeRoverMain.h
  ../../inc/VolumeRover2/RemapVoxels.h
  ../../inc/VolumeRover2/ImportData.h
  ../../inc/VolumeRover2/NewVolumeDialog.h
  ../../inc/VolumeRover2/VolumeInterface.h
  ../../inc/VolumeRover2/GeometryInterface.h
  ../../inc/VolumeRover2/ContoursInterface.h
  ../../inc/VolumeRover2/VolumeMemoryInterface.h
  ../../inc/VolumeRover2/DimensionModify.h
#  ../../inc/VolumeRover2/VolumeViewerPage.h
  ../../inc/VolumeRover2/BoundingBoxModify.h
  ../../inc/VolumeRover2/CVCMainWindow.h
  ../../inc/VolumeRover2/Viewers.h
  ../../inc/VolumeRover2/VolumeViewer.h
  ../../inc/VolumeRover2/ReconstructionDialog.h
  ../../inc/VolumeRover2/AnisotropicDiffusionDialog.h
  ../../inc/VolumeRover2/BilateralFilterDialog.h
  ../../inc/VolumeRover2/ContourTilerDialog.h
  ../../inc/VolumeRover2/ContrastEnhancementDialog.h
  ../../inc/VolumeRover2/GDTVFilterDialog.h
  ../../inc/VolumeRover2/PocketTunnelDialog.h
  ../../inc/VolumeRover2/CurationDialog.h
  ../../inc/VolumeRover2/TightCoconeDialog.h
  ../../inc/VolumeRover2/SuperSecondaryStructuresDialog.h
  ../../inc/VolumeRover2/SkeletonizationDialog.h
  ../../inc/VolumeRover2/LBIEQualityImprovementDialog.h
  ../../inc/VolumeRover2/MultiTileServerDialog.h
  ../../inc/VolumeRover2/VolumeViewer2.h
)

# the core CVC app functionality is now provided by libcvc (cvc::cvc)
set(LIBS cvc::cvc)

#
# arand: need to fix the CGAL dependence so that CGAL is included once
#        but whenever a library using CGAL is used
#
# joe - 6/10/2011 - added SetupCGAL to make this easier.


if(BUILD_POCKETTUNNEL_LIB AND NOT DISABLE_CGAL)
  set(LIBS ${LIBS} PocketTunnel)
endif(BUILD_POCKETTUNNEL_LIB AND NOT DISABLE_CGAL)

if(BUILD_SWEETMESH_LIB)
  set(LIBS ${LIBS} SweetMesh)
endif(BUILD_SWEETMESH_LIB)

if(BUILD_SKELETONIZATION_LIB AND NOT DISABLE_CGAL)
  set(LIBS ${LIBS} Skeletonization)
endif(BUILD_SKELETONIZATION_LIB AND NOT DISABLE_CGAL)

if(BUILD_TIGHTCOCONE_LIB AND NOT DISABLE_CGAL)
  set(LIBS ${LIBS} TightCocone)
endif(BUILD_TIGHTCOCONE_LIB AND NOT DISABLE_CGAL)

if(BUILD_SUPERSECONDARYSTRUCTURES_LIB AND NOT DISABLE_CGAL)
  set(LIBS ${LIBS} SuperSecondaryStructures)
endif(BUILD_SUPERSECONDARYSTRUCTURES_LIB AND NOT DISABLE_CGAL)

if(BUILD_CURATION_LIB AND NOT DISABLE_CGAL)
  set(LIBS ${LIBS} Curation)
endif(BUILD_CURATION_LIB AND NOT DISABLE_CGAL)

if(BUILD_SECONDARYSTRUCTURES_LIB AND NOT DISABLE_CGAL)
 set(LIBS ${LIBS} SecondaryStructures)
 set(LIBS ${LIBS} Histogram)
endif(BUILD_SECONDARYSTRUCTURES_LIB AND NOT DISABLE_CGAL)

if(BUILD_RECONSTRUCTION_LIB)
  set(LIBS ${LIBS} Reconstruction)
endif(BUILD_RECONSTRUCTION_LIB)
	
if(BUILD_HLEVELSET_LIB)
  set(LIBS ${LIBS} HLevelSet)
endif(BUILD_HLEVELSET_LIB)

if(BUILD_SEGMENTATION_LIB AND NOT DISABLE_CGAL)
  set(LIBS ${LIBS} SegCapsid SegMonomer SegMed SegSubunit SecStruct)
endif(BUILD_SEGMENTATION_LIB AND NOT DISABLE_CGAL)

if(BUILD_TILING_LIB AND NOT DISABLE_CGAL)
  set(LIBS ${LIBS} ContourTiler CGAL::CGAL gsl gslcblas mpfr gmp ContourFilterLib)
endif(BUILD_TILING_LIB AND NOT DISABLE_CGAL)

if(BUILD_MMHLS_LIB)
  set(LIBS ${LIBS} ExtMCMesher generateRawV generateMesh)
endif(BUILD_MMHLS_LIB)

if(BUILD_MSLEVELSET_LIB)
  FIND_PACKAGE(CUDA)
  INCLUDE_DIRECTORIES(
  ${CUDA_INCLUDE_DIRS}
  ${CUDA_CUT_INCLUDE_DIR}
  )
  set(LIBS ${LIBS}  ${CUDA_LIBRARIES} ${CUDA_CUT_LIBRARIES} MSLevelSet)
  if(BUILD_HOSEGMENTATION_LIB)
  set(LIBS ${LIBS} HOSegmentation)
  	if(BUILD_MPSEGMENTATION_LIB)
	set(LIBS ${LIBS} MPSegmentation)
	endif(BUILD_MPSEGMENTATION_LIB)
  endif(BUILD_HOSEGMENTATION_LIB)
endif(BUILD_MSLEVELSET_LIB)



# VolumeRover2 dependencies
set(LIBS
    ${LIBS}
    ColorTable2
    ${ISOCONTOURING_LIB}
 # VolumeViewer is deprecated as of 2/13/2011.  New VolumeViewer development should be done on the version
 # in the VolumeRover2 directory for now.  The reason for this is that it is taking advantage of the CVC::App
 # functionality now, and that is only available in VolumeRover2.  In the future we should create a new
 # library for core functionality such as CVC::App.  We should keep linking to VolumeViewer library for as long
 # as we need to keep building the old VolumeViewerPage for reference.
#
# arand, 4-25-2011: it looks like we don't need this anymore
#
#    VolumeViewer
    VolMagick
    cvcraw_geometry
    VolumeRenderer
    GeometryRenderer
    QGLViewer
    cvcalgo
    cvc::xmlrpc
#    Socket
    log4cplus::log4cplus
    $<IF:$<BOOL:${WIN32}>,PThreads4W::PThreads4W,Threads::Threads>
)

# Link against VolumeGridRover if available
if(BUILD_VOLUMEGRIDROVER AND NOT DISABLE_CGAL)
  include_directories(${CMAKE_BINARY_DIR}/src/VolumeGridRover)
  set(LIBS ${LIBS} VolumeGridRover)
endif()

# Source grouping for IDEs
source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})
source_group("UI Files" FILES ${UI_FORMS})

set(MACOSX_BUNDLE_ICON_FILE ${VolumeRover_SOURCE_DIR}/dat/VolRover.icns)

# AUTOMOC and AUTOUIC will handle MOC and UI file generation automatically
add_executable(VolumeRover2 MACOSX_BUNDLE
  ${SOURCE_FILES}
  ${INCLUDE_FILES}
  ${UI_FORMS}
)

# Generate wrapper script to run VolumeRover2 with clean environment
# This avoids conflicts with snap-provided libraries on Linux
if(UNIX AND NOT APPLE)
  set(WRAPPER_SCRIPT "${CMAKE_BINARY_DIR}/bin/run-volumerover2.sh")
  file(WRITE ${WRAPPER_SCRIPT}
"#!/bin/bash
# Wrapper script to run VolumeRover2 with clean environment
# Avoids conflicts with snap-provided libraries (GTK_PATH, LD_LIBRARY_PATH pollution)

# Get the directory where this script is located
SCRIPT_DIR=\"\$(cd \"\$(dirname \"\${BASH_SOURCE[0]}\")\" && pwd)\"

# Clean environment - remove snap-related variables that cause library conflicts
exec env -i \\
  HOME=\"\$HOME\" \\
  USER=\"\$USER\" \\
  LOGNAME=\"\$LOGNAME\" \\
  PATH=\"/usr/local/bin:/usr/bin:/bin:/usr/local/sbin:/usr/sbin:/sbin\" \\
  DISPLAY=\"\$DISPLAY\" \\
  XAUTHORITY=\"\$XAUTHORITY\" \\
  XDG_RUNTIME_DIR=\"\$XDG_RUNTIME_DIR\" \\
  TERM=\"\$TERM\" \\
  LANG=\"\$LANG\" \\
  \"\${SCRIPT_DIR}/VolumeRover2\" \"\$@\"
")
  # Make the script executable using execute_process at configure time
  execute_process(COMMAND chmod +x ${WRAPPER_SCRIPT})
endif()

# Tell AUTOUIC where to find .ui files when headers are in inc/ but .ui files are in src/
set_target_properties(VolumeRover2 PROPERTIES
  AUTOUIC_SEARCH_PATHS "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Link libraries - use modern target-based linking for Qt6
if(USE_QT6)
  target_link_libraries(VolumeRover2
    PRIVATE
      Qt6::Core
      Qt6::Gui
      Qt6::Widgets
      Qt6::Xml
      Qt6::OpenGL
      Qt6::OpenGLWidgets
      Boost::thread
      Boost::date_time
      Boost::regex
      Boost::filesystem
      ${LIBS}
  )
else()
  target_link_libraries(VolumeRover2
    PRIVATE
      ${QT_LIBRARIES}
      Boost::thread
      Boost::date_time
      Boost::regex
      Boost::filesystem
      ${LIBS}
  )
endif()

# Setup additional dependencies
if(BUILD_VOLUMEGRIDROVER AND NOT DISABLE_CGAL)
  SetupGSL(VolumeRover2)
endif()

if(NOT DISABLE_CGAL)
  SetupCGAL(VolumeRover2)
endif()

if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.4)
  if(APPLE)
    set_target_properties(VolumeRover2 PROPERTIES
      OUTPUT_NAME ${CMAKE_BUNDLE_NAME})
  endif(APPLE)
  set(CMAKE_INSTALL_DESTINATION_ARGS 
    BUNDLE DESTINATION "${CMAKE_BUNDLE_LOCATION}")
endif(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.4)

install(TARGETS VolumeRover2
        RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

# ── Windows: bundle all runtime DLLs so the package is self-contained ──
if(WIN32 AND USE_QT6)
  # 1) Copy vcpkg / imported-target runtime DLLs (boost, glew, fftw, gsl, …)
  #    $<TARGET_RUNTIME_DLLS:tgt> resolves every DLL that CMake knows about
  #    through IMPORTED targets linked to VolumeRover2.  (CMake ≥ 3.21)
  install(FILES $<TARGET_RUNTIME_DLLS:VolumeRover2> TYPE BIN)

  # 2) Deploy Qt DLLs + platform plugins via the official Qt6 deploy helper.
  #    This calls windeployqt under the hood.  (Qt ≥ 6.3)
  qt_generate_deploy_app_script(
    TARGET VolumeRover2
    OUTPUT_SCRIPT _qt_deploy_script
    NO_TRANSLATIONS
    NO_COMPILER_RUNTIME        # MSVC CRT handled by InstallRequiredSystemLibraries
  )
  install(SCRIPT ${_qt_deploy_script})
endif()
