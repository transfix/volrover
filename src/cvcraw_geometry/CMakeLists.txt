cmake_minimum_required(VERSION 3.16)

file(GLOB SOURCE_FILES *.cpp)
file(GLOB INCLUDE_FILES ../../inc/cvcraw_geometry/*.[Hh])

include_directories(
  ../../inc
)

set(LIBS LBIE) 

if(NOT DISABLE_CGAL)
  find_package(CGAL)
  if(CGAL_FOUND)
    include(${CGAL_USE_FILE})
    # need the following flags in case CGAL has some special compiler needs for this compiler
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_CXX_FLAGS_INIT}")
    set(LIBS ${LIBS} ${CGAL_LIBRARIES})
    add_definitions(-DCVCRAW_GEOMETRY_ENABLE_PROJECT)
    set(INCLUDE_FILES ${INCLUDE_FILES} ../../inc/cvcraw_geometry/project_verts.h)
  endif(CGAL_FOUND)
endif(NOT DISABLE_CGAL)

option(CVCRAW_GEOMETRY_ENABLE_BUNNY "Enable built-in bunny geometry test mesh" ON)
if(CVCRAW_GEOMETRY_ENABLE_BUNNY)
  add_definitions(-DCVCRAW_GEOMETRY_ENABLE_BUNNY)
endif(CVCRAW_GEOMETRY_ENABLE_BUNNY)

source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})

add_library(cvcraw_geometry STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

set_target_properties(cvcraw_geometry PROPERTIES
  AUTOMOC OFF
  AUTOUIC OFF
  AUTORCC OFF
)

target_link_libraries(cvcraw_geometry
  PUBLIC
    Boost::thread
    Boost::date_time
    Boost::regex
    Boost::filesystem
    ${LIBS}
)
