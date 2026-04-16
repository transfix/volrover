// Compatibility macros for log4cplus 2.x migration
// The vendored log4cplus 1.0.4 had a custom FUNCTION_LOGGER macro in logger.h.
// This header provides it for the system log4cplus 2.x package.

#ifndef __CVC_LOG4CPLUS_COMPAT_H__
#define __CVC_LOG4CPLUS_COMPAT_H__

#include <log4cplus/logger.h>
#include <boost/current_function.hpp>

#define FUNCTION_LOGGER log4cplus::Logger::getInstance(BOOST_CURRENT_FUNCTION)

#endif
