cmake_minimum_required(VERSION 3.16)

SetupQt()

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../../inc/Segmentation/SegMonomer/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${INCLUDE_FILES})

INCLUDE_DIRECTORIES(
  ../../../inc
)

ADD_LIBRARY (SegMonomer STATIC ${SOURCE_FILES} ${INCLUDE_FILES})
TARGET_LINK_LIBRARIES(SegMonomer PUBLIC cvc::xmlrpc)
