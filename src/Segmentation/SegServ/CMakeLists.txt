CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../../inc/Segmentation/SegServ/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${INCLUDE_FILES})

INCLUDE_DIRECTORIES(
  ../../../inc
)

ADD_EXECUTABLE(SegServ ${SOURCE_FILES} ${INCLUDE_FILES})
TARGET_LINK_LIBRARIES(SegServ
  GenSeg
  SegCapsid
  SegMed
  SegMonomer
  SegSubunit
  SecStruct
  cvc::xmlrpc
  VolMagick
  m
)
