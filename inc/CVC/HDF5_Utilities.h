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
    // Pull in everything from cvc::hdf5_utils
    using namespace cvc::hdf5_utils;

    // Non-ctx overloads that inject cvc::app::instance()
    inline bool isGroup(const std::string& hdf5_filename,
                        const std::string& hdf5_objname)
    {
      return cvc::hdf5_utils::isGroup(cvc::app::instance(),
                                      hdf5_filename, hdf5_objname);
    }

    inline bool isDataSet(const std::string& hdf5_filename,
                          const std::string& hdf5_objname)
    {
      return cvc::hdf5_utils::isDataSet(cvc::app::instance(),
                                        hdf5_filename, hdf5_objname);
    }

    inline bool objectExists(const std::string& hdf5_filename,
                             const std::string& hdf5_objname)
    {
      return cvc::hdf5_utils::objectExists(cvc::app::instance(),
                                           hdf5_filename, hdf5_objname);
    }

    inline void removeObject(const std::string& hdf5_filename,
                             const std::string& hdf5_objname)
    {
      cvc::hdf5_utils::removeObject(cvc::app::instance(),
                                    hdf5_filename, hdf5_objname);
    }

    inline void createHDF5File(const std::string& hdf5_filename)
    {
      cvc::hdf5_utils::createHDF5File(cvc::app::instance(), hdf5_filename);
    }

    inline void createGroup(const std::string& hdf5_filename,
                            const std::string& hdf5_objname,
                            bool replace = false)
    {
      cvc::hdf5_utils::createGroup(cvc::app::instance(),
                                   hdf5_filename, hdf5_objname, replace);
    }

    inline void createDataSet(const std::string& hdf5_filename,
                              const std::string& hdf5_objname,
                              const cvc::bounding_box& boundingBox,
                              const cvc::dimension& dimension,
                              cvc::data_type dataType,
                              const bool replace = false,
                              const bool createGroups = true)
    {
      cvc::hdf5_utils::createDataSet(cvc::app::instance(),
                                     hdf5_filename, hdf5_objname,
                                     boundingBox, dimension, dataType,
                                     replace, createGroups);
    }

    inline void createDataSet(const std::string& hdf5_filename,
                              const std::string& hdf5_objname,
                              const cvc::dimension& dimension,
                              cvc::data_type dataType,
                              const bool createGroups = true)
    {
      cvc::hdf5_utils::createDataSet(cvc::app::instance(),
                                     hdf5_filename, hdf5_objname,
                                     dimension, dataType, createGroups);
    }

    inline void createDataSet(const std::string& hdf5_filename,
                              const std::string& hdf5_objname,
                              const std::string& value,
                              bool createGroups = true)
    {
      cvc::hdf5_utils::createDataSet(cvc::app::instance(),
                                     hdf5_filename, hdf5_objname,
                                     value, createGroups);
    }

    inline cvc::dimension getObjectDimension(const std::string& hdf5_filename,
                                             const std::string& hdf5_objname)
    {
      return cvc::hdf5_utils::getObjectDimension(cvc::app::instance(),
                                                 hdf5_filename, hdf5_objname);
    }

    inline void setObjectDimension(const std::string& hdf5_filename,
                                   const std::string& hdf5_objname,
                                   const cvc::dimension& dim)
    {
      cvc::hdf5_utils::setObjectDimension(cvc::app::instance(),
                                          hdf5_filename, hdf5_objname, dim);
    }

    inline cvc::dimension getDataSetDimensionForBoundingBox(
        const std::string& hdf5_filename,
        const std::string& hdf5_objname,
        const cvc::bounding_box& subvolbox)
    {
      return cvc::hdf5_utils::getDataSetDimensionForBoundingBox(
          cvc::app::instance(), hdf5_filename, hdf5_objname, subvolbox);
    }

    inline cvc::dimension getDataSetDimension(
        const std::string& hdf5_filename,
        const std::string& hdf5_objname,
        const cvc::bounding_box& subvolbox,
        const cvc::dimension& maxdim = cvc::dimension(256, 256, 256))
    {
      return cvc::hdf5_utils::getDataSetDimension(
          cvc::app::instance(), hdf5_filename, hdf5_objname,
          subvolbox, maxdim);
    }

    inline cvc::bounding_box getObjectBoundingBox(
        const std::string& hdf5_filename,
        const std::string& hdf5_objname)
    {
      return cvc::hdf5_utils::getObjectBoundingBox(
          cvc::app::instance(), hdf5_filename, hdf5_objname);
    }

    inline void setObjectBoundingBox(const std::string& hdf5_filename,
                                     const std::string& hdf5_objname,
                                     const cvc::bounding_box& boundingBox)
    {
      cvc::hdf5_utils::setObjectBoundingBox(cvc::app::instance(),
                                            hdf5_filename, hdf5_objname,
                                            boundingBox);
    }

    inline double getDataSetMinimum(const std::string& hdf5_filename,
                                    const std::string& hdf5_objname)
    {
      return cvc::hdf5_utils::getDataSetMinimum(cvc::app::instance(),
                                                hdf5_filename, hdf5_objname);
    }

    inline double getDataSetMaximum(const std::string& hdf5_filename,
                                    const std::string& hdf5_objname)
    {
      return cvc::hdf5_utils::getDataSetMaximum(cvc::app::instance(),
                                                hdf5_filename, hdf5_objname);
    }

    inline std::string getDataSetInfo(const std::string& hdf5_filename,
                                      const std::string& hdf5_objname)
    {
      return cvc::hdf5_utils::getDataSetInfo(cvc::app::instance(),
                                             hdf5_filename, hdf5_objname);
    }

    inline cvc::data_type getDataSetType(const std::string& hdf5_filename,
                                         const std::string& hdf5_objname)
    {
      return cvc::hdf5_utils::getDataSetType(cvc::app::instance(),
                                             hdf5_filename, hdf5_objname);
    }

    inline std::vector<std::string> getChildObjects(
        const std::string& hdf5_filename,
        const std::string& hdf5_objname = "/",
        const std::string& filter = std::string())
    {
      return cvc::hdf5_utils::getChildObjects(cvc::app::instance(),
                                              hdf5_filename, hdf5_objname,
                                              filter);
    }

    // Attribute accessors (by filename + object path) -- inject ctx and
    // handle group vs dataset resolution.
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
