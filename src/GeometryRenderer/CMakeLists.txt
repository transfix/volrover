cmake_minimum_required(VERSION 3.16)

file(GLOB SOURCE_FILES *.cpp)
file(GLOB INCLUDE_FILES ../../inc/GeometryRenderer/*.[Hh])

include_directories(
  ../../inc
)

SetupQt()

set(LIBS CVC cvcraw_geometry QGLViewer GLEW::GLEW) 

source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})

add_library(GeometryRenderer STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

target_link_libraries(GeometryRenderer
  PRIVATE
    ${LIBS}
    Boost::thread
    Boost::date_time
    Boost::regex
    Boost::filesystem
)
