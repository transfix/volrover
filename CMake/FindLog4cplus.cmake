# FindLog4cplus.cmake
# Finds the log4cplus logging library.
#
# Defines:
#   log4cplus_FOUND        - True if log4cplus was found
#   log4cplus::log4cplus   - Imported target
#
# First tries config mode (log4cplusConfig.cmake), then falls back
# to manual detection via headers and library files.

# Try config mode first (works on Linux with liblog4cplus-dev, vcpkg, etc.)
find_package(log4cplus CONFIG QUIET)
if(log4cplus_FOUND)
  return()
endif()

# Fallback: find headers and library manually (e.g., Homebrew on macOS)
find_path(LOG4CPLUS_INCLUDE_DIR
  NAMES log4cplus/logger.h
  HINTS
    ${LOG4CPLUS_ROOT}
    ENV LOG4CPLUS_ROOT
  PATH_SUFFIXES include
)

find_library(LOG4CPLUS_LIBRARY
  NAMES log4cplus
  HINTS
    ${LOG4CPLUS_ROOT}
    ENV LOG4CPLUS_ROOT
  PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Log4cplus
  REQUIRED_VARS LOG4CPLUS_LIBRARY LOG4CPLUS_INCLUDE_DIR
)

if(Log4cplus_FOUND AND NOT TARGET log4cplus::log4cplus)
  add_library(log4cplus::log4cplus UNKNOWN IMPORTED)
  set_target_properties(log4cplus::log4cplus PROPERTIES
    IMPORTED_LOCATION "${LOG4CPLUS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${LOG4CPLUS_INCLUDE_DIR}"
  )
endif()

# Propagate to parent scope variable expected by find_package
set(log4cplus_FOUND ${Log4cplus_FOUND})
