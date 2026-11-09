cmake_minimum_required(VERSION 3.16)
project(QGLViewer VERSION 3.0.0)

# Setup Qt6
if(USE_QT6)
  find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets Xml OpenGL OpenGLWidgets)
  set(QT_LIBRARIES Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Xml Qt6::OpenGL Qt6::OpenGLWidgets)
else()
  SetupQt()
endif()

# Enable Qt MOC, UIC, and RCC
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Windows-specific definitions
if(WIN32)
  add_definitions(-DWIN32 -DNOMINMAX)
  if(MSVC)
    add_compile_options(/EHsc /FS)
  endif()
endif()

# Optional stereo display support
if(ENABLE_STEREO_DISPLAY)
  add_definitions(-DENABLE_STEREO_DISPLAY)
endif()

# QGLViewer header files (from inc directory)
set(QGL_HEADERS
  ../../inc/QGLViewer/qglviewer.h
  ../../inc/QGLViewer/camera.h
  ../../inc/QGLViewer/manipulatedFrame.h
  ../../inc/QGLViewer/manipulatedCameraFrame.h
  ../../inc/QGLViewer/frame.h
  ../../inc/QGLViewer/constraint.h
  ../../inc/QGLViewer/keyFrameInterpolator.h
  ../../inc/QGLViewer/mouseGrabber.h
  ../../inc/QGLViewer/quaternion.h
  ../../inc/QGLViewer/vec.h
  ../../inc/QGLViewer/domUtils.h
  ../../inc/QGLViewer/config.h
)

# QGLViewer source files
set(SOURCES
  qglviewer.cpp
  camera.cpp
  manipulatedFrame.cpp
  manipulatedCameraFrame.cpp
  frame.cpp
  saveSnapshot.cpp
  constraint.cpp
  keyFrameInterpolator.cpp
  mouseGrabber.cpp
  quaternion.cpp
  vec.cpp
)

# UI forms (using Qt4-compatible versions for Qt6)
set(UI_FORMS ImageInterface.ui)

# VRender vectorial rendering support (EPS/PS/XFIG snapshots). Disabled by
# default: VRender is GPL-2.0-or-later and includes the Generic Polygon
# Clipper (non-commercial use only), so building it makes the binary a GPL
# combined work. See THIRD_PARTY_NOTICES.md. Pass -DNO_VECTORIAL_RENDER=OFF
# to build it anyway.
option(NO_VECTORIAL_RENDER "Disable vectorial rendering support (VRender, GPL)" ON)
if(NOT NO_VECTORIAL_RENDER)
  set(VRENDER_SOURCES
    VRender/BackFaceCullingOptimizer.cpp
    VRender/BSPSortMethod.cpp
    VRender/EPSExporter.cpp
    VRender/Exporter.cpp
    VRender/FIGExporter.cpp
    VRender/gpc.cpp
    VRender/ParserGL.cpp
    VRender/Primitive.cpp
    VRender/PrimitivePositioning.cpp
    VRender/TopologicalSortMethod.cpp
    VRender/VisibilityOptimizer.cpp
    VRender/Vector2.cpp
    VRender/Vector3.cpp
    VRender/NVector3.cpp
    VRender/VRender.cpp
  )

  set(VRENDER_HEADERS
    ../../inc/QGLViewer/VRender/AxisAlignedBox.h
    ../../inc/QGLViewer/VRender/Exporter.h
    ../../inc/QGLViewer/VRender/gpc.h
    ../../inc/QGLViewer/VRender/NVector3.h
    ../../inc/QGLViewer/VRender/Optimizer.h
    ../../inc/QGLViewer/VRender/ParserGL.h
    ../../inc/QGLViewer/VRender/Primitive.h
    ../../inc/QGLViewer/VRender/PrimitivePositioning.h
    ../../inc/QGLViewer/VRender/SortMethod.h
    ../../inc/QGLViewer/VRender/Types.h
    ../../inc/QGLViewer/VRender/Vector2.h
    ../../inc/QGLViewer/VRender/Vector3.h
    ../../inc/QGLViewer/VRender/VRender.h
  )

  list(APPEND SOURCES ${VRENDER_SOURCES})
  list(APPEND QGL_HEADERS ${VRENDER_HEADERS})
  list(APPEND UI_FORMS VRenderInterface.ui)
else()
  add_definitions(-DNO_VECTORIAL_RENDER)
endif()

# Optional GLUT support
find_package(GLUT)
if(GLUT_FOUND)
  include_directories(${GLUT_INCLUDE_DIR})
endif()

# Source groups for IDE organization
source_group("Source Files" FILES ${SOURCES})
source_group("Header Files" FILES ${QGL_HEADERS})
source_group("UI Files" FILES ${UI_FORMS})

# Create the library (shared by default, can be made static)
option(QGLVIEWER_STATIC "Build QGLViewer as a static library" OFF)
if(QGLVIEWER_STATIC)
  add_library(QGLViewer STATIC ${SOURCES} ${QGL_HEADERS} ${UI_FORMS})
  target_compile_definitions(QGLViewer PUBLIC QGLVIEWER_STATIC)
else()
  add_library(QGLViewer SHARED ${SOURCES} ${QGL_HEADERS} ${UI_FORMS})
  target_compile_definitions(QGLViewer PRIVATE CREATE_QGLVIEWER_DLL)
endif()

# Set library properties
set_target_properties(QGLViewer PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  POSITION_INDEPENDENT_CODE ON
)

# Include directories
target_include_directories(QGLViewer
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../inc>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Link Qt and other libraries
target_link_libraries(QGLViewer
  PUBLIC
    ${QT_LIBRARIES}
)

# Link OpenGL and GLU
if(WIN32)
  target_link_libraries(QGLViewer PUBLIC opengl32 glu32)
elseif(APPLE)
  # GLU is part of the OpenGL framework on macOS
  target_link_libraries(QGLViewer PUBLIC "-framework OpenGL")
else()
  # Unix/Linux
  find_package(OpenGL REQUIRED)
  target_link_libraries(QGLViewer PUBLIC OpenGL::GL OpenGL::GLU)
endif()

# Optional GLUT linking
if(GLUT_FOUND)
  target_link_libraries(QGLViewer PUBLIC ${GLUT_LIBRARIES})
  if(APPLE)
    target_link_libraries(QGLViewer PUBLIC "-framework GLUT" "-lobjc")
  endif()
endif()

# Ship the shared library with the app. Windows gets QGLViewer.dll through
# VolumeRover2's $<TARGET_RUNTIME_DLLS> install, so only the other platforms
# need it here. On Linux its dependencies sit beside it in <prefix>/lib
# (DT_RUNPATH is not inherited from the executable); macOS dyld also
# searches the app's rpaths.
if(NOT QGLVIEWER_STATIC AND NOT WIN32)
  if(NOT APPLE)
    set_target_properties(QGLViewer PROPERTIES INSTALL_RPATH "$ORIGIN")
  endif()
  install(TARGETS QGLViewer LIBRARY DESTINATION lib)
endif()

# GLEW linking (required for OpenGL extensions)
if(GLEW_FOUND)
  target_link_libraries(QGLViewer PUBLIC GLEW::GLEW)
endif()
