cmake_minimum_required(VERSION 3.16)

file(GLOB INCLUDE_FILES ../../inc/CVC/*.[Hh])
set(SOURCE_FILES
  App.cpp
  State.cpp
)

include_directories(
  ../../inc
)

# CVC_NAMESPACE is now set by libcvc's PUBLIC compile definitions (=cvc)
# No need to override here.

# The libCVC version number.
set(CVC_VERSION_MAJOR 1)
set(CVC_VERSION_MINOR 0)
set(CVC_VERSION_PATCH ${CVC_GIT_REVISION})

# The full version string
set(CVC_VERSION ${CVC_VERSION_MAJOR}.${CVC_VERSION_MINOR}.${CVC_VERSION_PATCH})
message(STATUS "libCVC version: ${CVC_VERSION}")

add_definitions(-DCVC_VERSION_STRING=\"${CVC_VERSION}\")

set(CVC_LINK_LIBS "")
option(CVC_USING_XMLRPC "If this is enabled, CVC::State will support state sharing with other processes on the network" OFF)
if(CVC_USING_XMLRPC)
  set(CVC_LINK_LIBS "XmlRPC")
  add_definitions(-DUSING_XMLRPC)
endif(CVC_USING_XMLRPC)

option(CVC_USING_HDF5 "HDF5 Support. Needed for *.cvc format but requires linking to libhdf5." ON)
mark_as_advanced(CVC_USING_HDF5)
if(CVC_USING_HDF5)
  find_package(HDF5 COMPONENTS C CXX)
  if(HDF5_FOUND)
    message("HDF5 found! Enabling CVC file support")

    add_definitions(-DCVC_USING_HDF5 ${HDF5_DEFINITIONS})
    set(INCLUDE_FILES
       ${INCLUDE_FILES}
       ../../inc/CVC/HDF5_Utilities.h
    )
    set(SOURCE_FILES
       ${SOURCE_FILES}
       HDF5_Utilities.cpp
    )

    # Collect all HDF5 include directories
    # Modern FindHDF5 uses HDF5_INCLUDE_DIRS, but also check legacy HDF5_INCLUDE_DIR
    # On some systems, manually add the serial path
    set(CVC_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS} ${HDF5_C_INCLUDE_DIRS} ${HDF5_CXX_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR})
    if(EXISTS "/usr/include/hdf5/serial")
      list(APPEND CVC_HDF5_INCLUDE_DIRS "/usr/include/hdf5/serial")
    endif()
    # Remove duplicates
    if(CVC_HDF5_INCLUDE_DIRS)
      list(REMOVE_DUPLICATES CVC_HDF5_INCLUDE_DIRS)
    endif()

    set(CVC_LINK_LIBS ${CVC_LINK_LIBS} ${HDF5_LIBRARIES})
  else(HDF5_FOUND)
    message("HDF5 not found! CVC file support disabled")
    set(CVC_HDF5_DISABLED TRUE)
  endif(HDF5_FOUND)
else(CVC_USING_HDF5)
  set(CVC_HDF5_DISABLED TRUE)
endif(CVC_USING_HDF5)

#option(CVC_LOG4CPLUS_DEFAULT "Use log4cplus output by default" OFF)
#mark_as_advanced(CVC_LOG4CPLUS_DEFAULT)
#if(CVC_LOG4CPLUS_DEFAULT)
#  find_package(Log4cplus)
#  if(LOG4CPLUS_FOUND)
#    add_definitions(-DUSING_LOG4CPLUS_DEFAULT)
#    include_directories(${LOG4CPLUS_INCLUDE_DIR})
#    set(CVC_LINK_LIBS ${CVC_LINK_LIBS} ${LOG4CPLUS_LIBRARIES})
#  else(LOG4CPLUS_FOUND)
#    message("Log4cplus not found! Log4cplus based logging disabled")
#  endif(LOG4CPLUS_FOUND)
#endif(CVC_LOG4CPLUS_DEFAULT

set(CVC_LINK_LIBS ${CVC_LINK_LIBS} log4cplus::log4cplus)

option(LOG4CPLUS_DISABLE_TRACE "Disable logging trace" OFF)
mark_as_advanced(LOG4CPLUS_DISABLE_TRACE)

option(CVC_APP_XML_PROPERTY_MAP "Read/Write XML for property map info.  If not set, uses boost INFO property tree format." OFF)
if(CVC_APP_XML_PROPERTY_MAP)
  add_definitions(-DCVC_APP_XML_PROPERTY_MAP)
endif(CVC_APP_XML_PROPERTY_MAP)
mark_as_advanced(CVC_APP_XML_PROPERTY_MAP)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/../../inc/CVC/config.h.cmake ${CMAKE_CURRENT_SOURCE_DIR}/../../inc/CVC/config.h )


source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})

add_library(CVC STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

# Add HDF5 include directories if enabled
if(CVC_USING_HDF5 AND HDF5_FOUND)
  target_include_directories(CVC PUBLIC ${CVC_HDF5_INCLUDE_DIRS})
endif()

set_target_properties(CVC PROPERTIES
  AUTOMOC OFF
  AUTOUIC OFF
  AUTORCC OFF
)

target_link_libraries(CVC
  PUBLIC
    cvc::cvc
    Boost::thread
    Boost::date_time
    Boost::regex
    Boost::filesystem
    ${CVC_LINK_LIBS}
)

# Enable the following later once we figure out how to do component installs
# to separate runtime and development distributions. - Joe R. 5/20/2011

#install(TARGETS CVC EXPORT CVC
#        RUNTIME DESTINATION bin
#        LIBRARY DESTINATION lib
#        ARCHIVE DESTINATION lib)
#install(EXPORT CVC NAMESPACE cvc_ DESTINATION lib/CVC)
