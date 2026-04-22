/*
  Compatibility header: forwards to libcvc's cvc::hdf5_utils namespace
  and exposes it as CVC::HDF5_Utilities for legacy code.
*/

#ifndef __CVC_COMPAT_HDF5_UTILITIES_H__
#define __CVC_COMPAT_HDF5_UTILITIES_H__

#include <cvc/hdf5_utils.h>
#include <cvc/app.h>
#include <CVC/Namespace.h>

namespace CVC_NAMESPACE
{
  namespace HDF5_Utilities
  {
    // Pull in everything from cvc::hdf5_utils (includes ctx-less overloads).
    using namespace cvc::hdf5_utils;

    // Attribute accessors keyed by filename + object path
    // (not provided by libcvc; inject ctx and group/dataset resolution).
    template<class T>
    inline void getAttribute(const std::string& hdf5_filename,
                             const std::string& hdf5_objname,
                             const std::string& attrname,
                             T& value)
    {
      auto f = cvc::hdf5_utils::getH5File(hdf5_filename);
      if (cvc::hdf5_utils::isGroup(cvc::app::instance(),
                                   hdf5_filename, hdf5_objname))
        {
          auto g = cvc::hdf5_utils::getGroup(*f, hdf5_objname, false);
          cvc::hdf5_utils::getAttribute(*g, attrname, value);
        }
      else
        {
          auto d = cvc::hdf5_utils::getDataSet(*f, hdf5_objname, false);
          cvc::hdf5_utils::getAttribute(*d, attrname, value);
        }
    }

    template<class T>
    inline void setAttribute(const std::string& hdf5_filename,
                             const std::string& hdf5_objname,
                             const std::string& attrname,
                             const T& value)
    {
      auto f = cvc::hdf5_utils::getH5File(hdf5_filename);
      if (cvc::hdf5_utils::isGroup(cvc::app::instance(),
                                   hdf5_filename, hdf5_objname))
        {
          auto g = cvc::hdf5_utils::getGroup(*f, hdf5_objname, false);
          cvc::hdf5_utils::setAttribute(*g, attrname, value);
        }
      else
        {
          auto d = cvc::hdf5_utils::getDataSet(*f, hdf5_objname, false);
          cvc::hdf5_utils::setAttribute(*d, attrname, value);
        }
    }

    inline void setAttribute(const std::string& hdf5_filename,
                             const std::string& hdf5_objname,
                             const std::string& attrname,
                             const char* value)
    {
      setAttribute<std::string>(hdf5_filename, hdf5_objname,
                                attrname, std::string(value));
    }
  } // namespace HDF5_Utilities
} // namespace CVC_NAMESPACE

#endif
