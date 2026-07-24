cmake_minimum_required(VERSION 3.16)

SetupQt()

file(GLOB SOURCE_FILES *.cpp)
file(GLOB INCLUDE_FILES ../../inc/SuperSecondaryStructures/*.[Hh])

source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})

include_directories(
  ../../inc
)

add_definitions(-DHUGE=1000000)

add_library(SuperSecondaryStructures STATIC ${SOURCE_FILES} ${INCLUDE_FILES})
target_link_libraries(SuperSecondaryStructures PUBLIC VolMagick)
target_compile_definitions(SuperSecondaryStructures PRIVATE CVCGEOM_NAMESPACE=cvcraw_geometry)
SetupCGAL(SuperSecondaryStructures)
