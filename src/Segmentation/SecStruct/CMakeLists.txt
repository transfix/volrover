cmake_minimum_required(VERSION 3.16)

SetupQt()

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../../inc/Segmentation/SecStruct/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${INCLUDE_FILES})

INCLUDE_DIRECTORIES(
  ../../../inc
)

OPTION(VRML_OUTPUT "Secondary structures ouput as VRML file format." OFF)
IF(VRML_OUTPUT)
	ADD_DEFINITIONS(-DVRML_OUTPUT) 
ENDIF(VRML_OUTPUT)

add_library(SecStruct STATIC ${SOURCE_FILES} ${INCLUDE_FILES})
TARGET_LINK_LIBRARIES(SecStruct PUBLIC cvc::xmlrpc)
SetupCGAL(SecStruct)
