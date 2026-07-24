cmake_minimum_required(VERSION 3.16)

SetupQt()

IF(WIN32)
  ADD_DEFINITIONS(-D__WINDOWS__)
ENDIF(WIN32)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../../inc/Segmentation/GenSeg/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${INCLUDE_FILES})

INCLUDE_DIRECTORIES(
  ../../../inc
)

ADD_LIBRARY (GenSeg STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

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

TARGET_LINK_LIBRARIES(GenSeg VolMagick cvc::xmlrpc)
