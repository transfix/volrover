// Compatibility macros for log4cplus 2.x migration
// The vendored log4cplus 1.0.4 had a custom FUNCTION_LOGGER macro in logger.h.
// This header provides it for the system log4cplus 2.x package.
//
// On Windows, Qt6 defines UNICODE which causes log4cplus headers to use
// wchar_t (tstring = wstring). But vcpkg builds log4cplus without UNICODE,
// so we must ensure headers are parsed in char mode to match the library ABI.

#ifndef __CVC_LOG4CPLUS_COMPAT_H__
#define __CVC_LOG4CPLUS_COMPAT_H__

#ifdef _WIN32
#  ifdef UNICODE
#    define _CVC_SAVED_UNICODE
#    undef UNICODE
#  endif
#  ifdef _UNICODE
#    define _CVC_SAVED__UNICODE
#    undef _UNICODE
#  endif
#endif

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/configurator.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/loglevel.h>

#ifdef _WIN32
#  ifdef _CVC_SAVED_UNICODE
#    define UNICODE
#    undef _CVC_SAVED_UNICODE
#  endif
#  ifdef _CVC_SAVED__UNICODE
#    define _UNICODE
#    undef _CVC_SAVED__UNICODE
#  endif
#endif

#include <boost/current_function.hpp>

#define FUNCTION_LOGGER log4cplus::Logger::getInstance(BOOST_CURRENT_FUNCTION)

#endif
