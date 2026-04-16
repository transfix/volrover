cmake_minimum_required(VERSION 3.16)

SetupQt()

#SET(CONTOURTILER_NAMESPACE      CONTOURTILER CACHE STRING "Namespace to use for contour tiler")

#MARK_AS_ADVANCED(FORCE CONTOURTILER_NAMESPACE)

#Make sure we've got boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.41.0 COMPONENTS filesystem system regex thread)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  message("Boost includes: ${Boost_INCLUDE_DIRS}")
  message("Boost libraries: ${Boost_LIBRARIES}")
else(Boost_FOUND)
  message("If you're having trouble finding boost, set environment variables "
          "BOOST_INCLUDEDIR and BOOST_LIBRARYDIR to the appropriate paths")
endif(Boost_FOUND)

#include_directories(/org/centers/cvc/software/share/usr.linux.x86_64/include)

#OPTION( CONTOURTILER_EXACT_ARITHMETIC "Enable exact arithmetic (slow)" OFF )
#OPTION( LOG4CPLUS_DISABLE_TRACE "Disable logging trace" ON )
#CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/../../inc/ContourTiler/config_log.h.cmake ${CMAKE_CURRENT_SOURCE_DIR}/../../inc/ContourTiler/config_log.h )

SET(SRCS 
  CGAL_hash.cpp
  Contour.cpp 
  Contour2.cpp 
  contour_graph.cpp
  Hierarchy.cpp 
  tiler_operations.cpp 
  contour_utils.cpp 
  polygon_utils.cpp 
  Correspondences.cpp 
  Tiling_region.cpp 
  Composite_tiling_region.cpp 
  All_tiling_region.cpp 
  Point_tiling_region.cpp 
  Wedge.cpp 
  print_utils.cpp
  skeleton.cpp
  Point_25.cpp
  Boundary_slice_chords.cpp
  theorems.cpp
  projection.cpp
  Tiler_workspace.cpp
  tiler.cpp
  tiler_output.cpp
  Tiles.cpp
  mtiler_operations.cpp
  segment_utils.cpp
  intersection.cpp
  cut.cpp
  augment.cpp
  Polyline.cpp
  triangle_utils.cpp
  remove_contour_intersections.cpp
  offset_polygon.cpp
  polygon_difference.cpp
  polygon_intersection.cpp
  Slice.cpp
  Slice2.cpp
  minkowski.cpp
  Z_adjustments.cpp
  interp.cpp
  arrange_polygons.cpp
  perturb.cpp
  cl_options.cpp
  ecs.cpp
  smfconvert.cpp
  dev.cpp
  xml/Ser_reader.cpp
  xml/tinyxml.cpp
  xml/tinystr.cpp
  xml/tinyxmlerror.cpp
  xml/tinyxmlparser.cpp
)

SET(ECS_SRCS 
  ecs_main.cpp 
)

SET(SMF_SRCS 
  smfconvert_main.cpp 
)

SET(TEST_SRCS 
  test_contour_utils.cpp 
  test_tile.cpp 
  test_tile2.cpp 
  test_Tiling_region.cpp 
  test_Hierarchy.cpp 
  test_reader_gnuplot.cpp
  test_mtile.cpp
  test_intersection.cpp
  test_augment.cpp
  test_sweep_line.cpp
  test_skeleton.cpp
  test_minkowski.cpp
)

INCLUDE_DIRECTORIES (
#  ../../inc/ContourTiler
#  ../../inc/CppUnitLite2
  ../../inc
)

#message ("CGAL_LIBRARIES: ${CGAL_LIBRARIES}")
#get_directory_property(test LINK_DIRECTORIES)
#message ("LINK_DIRECTORIES: ${test}")

SET(LINK_LIBS 
  log4cplus::log4cplus
#  CGAL_Core
#  CGAL 
#  CGALcore++
#  gsl 
#  gslcblas 
#  mpfr 
#  gmp 
#  Tiling
#  cppunitlite
#  ContourFilterLib
)

#find_package(log4cplus)
#if(LOG4CPLUS_FOUND)
#  include_directories(${LOG4CPLUS_INCLUDE_DIR})
#  set(LINK_LIBS ${LINK_LIBS} ${LOG4CPLUS_LIBRARIES})
#else(LOG4CPLUS_FOUND)
#  message(SEND_ERROR "ContourTiler requires log4cplus!")
#endif(LOG4CPLUS_FOUND)

#find_package(CGAL 4.0)
#if(CGAL_FOUND)
#  include(${CGAL_USE_FILE})
#  # need the following flags in case CGAL has some special compiler needs for this compiler
##  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_CXX_FLAGS_INIT}")
#  set(LIBS ${LIBS} ${CGAL_LIBRARIES})
#else(CGAL_FOUND)
#  message(SEND_ERROR "ContourTiler requires CGAL!")
#endif(CGAL_FOUND)

#ADD_DEFINITIONS(-DCONTOURTILER_NAMESPACE=${CONTOURTILER_NAMESPACE} "-DCONTOURTILER_BEGIN_NAMESPACE=namespace ${CONTOURTILER_NAMESPACE} {" -DCONTOURTILER_END_NAMESPACE=})
ADD_DEFINITIONS(-DCGAL_CFG_NO_CPP0X_RVALUE_REFERENCE -DCGAL_CFG_NO_CPP0X_ARRAY -DCGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES -DCGAL_CFG_NO_CPP0X_TUPLE)
ADD_DEFINITIONS(-DTIXML_USE_STL)

#IF(CMAKE_BUILD_TYPE MATCHES Debug)
#  ADD_DEFINITIONS(-DDEBUG)
#ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

#SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -lpthread -I/h1/edwardsj/projects/CGAL-3.6/include -std=c++0x")
#SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -lpthread -I/h1/edwardsj/projects/CGAL-3.6/include -lpthread -pthread -frounding-math")
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -lpthread -pthread -frounding-math")

ADD_LIBRARY(ContourTiler STATIC ${SRCS})
ADD_EXECUTABLE (ContourTilerBin main.cpp)
#ADD_EXECUTABLE (ContourTilerTest ${TEST_SRCS} test_main.cpp)
ADD_EXECUTABLE (ContourTilerSandbox sandbox.cpp)
ADD_EXECUTABLE (ContourTilerDecimate decimate.cpp)
ADD_EXECUTABLE (ecs ${ECS_SRCS})
ADD_EXECUTABLE (smfconvert ${SMF_SRCS})
ADD_EXECUTABLE (ContourNear contour_near_main.cpp)

# CGAL must be set up for each individual project
SetupCGAL(ContourTilerBin)
SetupCGAL(ContourTilerDecimate)
SetupCGAL(ContourTilerSandbox)
SetupCGAL(ecs)
SetupCGAL(smfconvert)
SetupCGAL(ContourNear)

SET(LINK_LIBS log4cplus::log4cplus)

TARGET_LINK_LIBRARIES(ContourTilerBin PUBLIC ContourTiler ${LINK_LIBS} Boost::filesystem Boost::regex Boost::thread)
#TARGET_LINK_LIBRARIES(ContourTilerTest ContourTiler ${LINK_LIBS} Boost::filesystem Boost::regex Boost::thread)
TARGET_LINK_LIBRARIES(ContourTilerSandbox PUBLIC ContourTiler ${LINK_LIBS} Boost::filesystem Boost::regex Boost::thread)
TARGET_LINK_LIBRARIES(ContourTilerDecimate PUBLIC ContourTiler ${LINK_LIBS} Boost::filesystem Boost::regex Boost::thread)
TARGET_LINK_LIBRARIES(ecs PUBLIC ContourTiler ${LINK_LIBS} Boost::filesystem Boost::regex Boost::thread)
TARGET_LINK_LIBRARIES(smfconvert PUBLIC ContourTiler ${LINK_LIBS} Boost::filesystem Boost::regex Boost::thread)
TARGET_LINK_LIBRARIES(ContourNear PUBLIC ContourTiler ContourFilterLib ${LINK_LIBS} Boost::filesystem Boost::regex Boost::thread)

SET_TARGET_PROPERTIES(ContourTiler PROPERTIES COMPILE_FLAGS "-rdynamic")# OUTPUT_NAME "ContourTiler")
#SET_TARGET_PROPERTIES(ContourTilerBin PROPERTIES OUTPUT_NAME "ContourTiler")

# a test target
#ADD_CUSTOM_TARGET(run_unit_tests
#  ALL
#  COMMAND ${CMAKE_BINARY_DIR}/ContourTilerTest --data ../test_data --out ./output
#  DEPENDS unit_tests
#  VERBATIM
#)

#ADD_CUSTOM_TARGET(CTRunTest
#  ALL
#  COMMAND ${CMAKE_BINARY_DIR}/bin/ContourTilerTest --data ../test_data --out output --test tile8
#  DEPENDS ContourTilerTest
#  VERBATIM
#)

# GDB commands
# run --data ./test_data --out ./debug/output --test augment1
# run --data ./test_data --out ./debug/output --test print11_12
