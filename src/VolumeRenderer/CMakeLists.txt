# 10/07/2011 -- transfix - added libCVC dependency

cmake_minimum_required(VERSION 3.16)

set(SOURCE_FILES
  ClipCube.cpp
  Extent.cpp
  FragmentProgramARBImpl.cpp
  FragmentProgramImpl.cpp
  FragmentProgramGLSLImpl.cpp
  MyExtensions.cpp
  Paletted2DImpl.cpp
  PalettedImpl.cpp
  Plane.cpp
  PolygonArray.cpp
  Polygon.cpp
  RendererBase.cpp
  Renderer.cpp
  RGBABase.cpp
  SGIColorTableImpl.cpp
  SimpleRGBA2DImpl.cpp
  SimpleRGBAImpl.cpp
  UnshadedBase.cpp
  VolumeRenderer.cpp
  VolumeRendererFactory.cpp
)
set(INCLUDE_FILES
  ../../inc/VolumeRenderer/ClipCube.h
  ../../inc/VolumeRenderer/ExtensionPointers.h
  ../../inc/VolumeRenderer/Extent.h
  ../../inc/VolumeRenderer/FragmentProgramARBImpl.h
  ../../inc/VolumeRenderer/FragmentProgramImpl.h
  ../../inc/VolumeRenderer/FragmentProgramGLSLImpl.h
  ../../inc/VolumeRenderer/LookupTables.h
  ../../inc/VolumeRenderer/MyExtensions.h
  ../../inc/VolumeRenderer/Paletted2DImpl.h
  ../../inc/VolumeRenderer/PalettedImpl.h
  ../../inc/VolumeRenderer/Plane.h
  ../../inc/VolumeRenderer/PolygonArray.h
  ../../inc/VolumeRenderer/Polygon.h
  ../../inc/VolumeRenderer/RendererBase.h
  ../../inc/VolumeRenderer/Renderer.h
  ../../inc/VolumeRenderer/RGBABase.h
  ../../inc/VolumeRenderer/SGIColorTableImpl.h
  ../../inc/VolumeRenderer/SimpleRGBA2DImpl.h
  ../../inc/VolumeRenderer/SimpleRGBAImpl.h
  ../../inc/VolumeRenderer/StaticExtensionPointers.h
  ../../inc/VolumeRenderer/UnshadedBase.h
  ../../inc/VolumeRenderer/VolumeRendererFactory.h
  ../../inc/VolumeRenderer/VolumeRenderer.h
)

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${INCLUDE_FILES})

if(APPLE)
  add_definitions(-D__APPLE__)
endif(APPLE)
if(WIN32)
  add_definitions(-DWIN32)
endif(WIN32)

find_package(OpenGL)
message("OPENGL_INCLUDE_DIR: ${OPENGL_INCLUDE_DIR}")

# For some reason, OPENGL_INCLUDE_DIR is finding the framework dir and not the actual header dir on MacOSX.
# Hence why we have the following... - Joe R. - 2010May17
set(ACTUAL_OPENGL_INCLUDE_DIR ${OPENGL_INCLUDE_DIR})
if(APPLE)
  set(ACTUAL_OPENGL_INCLUDE_DIR ${OPENGL_INCLUDE_DIR}/Headers)
endif(APPLE)

INCLUDE_DIRECTORIES(
  ../../inc
  ${ACTUAL_OPENGL_INCLUDE_DIR}
)

# the core CVC app functionality is now provided by libcvc (cvc::cvc)
set(LIBS cvc::cvc)

find_package(Cg)
if(CG_FOUND)
  message("Cg found!")
  option(VOLUMERENDERER_USING_CG "Use CG with the VolumeRenderer.  Currently needed for shaded volume rendering" ON)
  mark_as_advanced(VOLUMERENDERER_USING_CG)
  if(VOLUMERENDERER_USING_CG)
    add_definitions(-DCG)
    include_directories(
      ${CG_INCLUDE_PATH}
    )
    set(SOURCE_FILES
      ${SOURCE_FILES}
      CG_Programs.cpp
      CGImpl.cpp
      CGRGBAImpl.cpp
    )
    set(INCLUDE_FILES
      ${INCLUDE_FILES}
      ../../inc/VolumeRenderer/CG_Programs.h
      ../../inc/VolumeRenderer/CGImpl.h
      ../../inc/VolumeRenderer/CGRGBAImpl.h
    )
    set(LIBS ${LIBS} ${CG_LIBRARY} ${CG_GL_LIBRARY})
  endif(VOLUMERENDERER_USING_CG)
else(CG_FOUND)
  message("Cg NOT found!")
endif(CG_FOUND)

ADD_LIBRARY (VolumeRenderer STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

set_target_properties(VolumeRenderer PROPERTIES
  AUTOMOC OFF
  AUTOUIC OFF
  AUTORCC OFF
)

#Make sure we've got boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.34.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  message("Boost includes: ${Boost_INCLUDE_DIRS}")
else(Boost_FOUND)
  message("If you're having trouble finding boost, set environment variables "
          "BOOST_INCLUDEDIR and BOOST_LIBRARYDIR to the appropriate paths")
endif(Boost_FOUND)

set(LIBS ${LIBS} GLEW::GLEW)

target_link_libraries(VolumeRenderer ${LIBS})

