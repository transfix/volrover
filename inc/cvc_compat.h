/*
  cvc_compat.h — single volrover-side compatibility shim for libcvc 3.1.0+.

  Volrover is a legacy application with thousands of references to
  CVC::App, CVC::ThreadInfo, the cvcapp / cvcstate macros, etc.
  libcvc 3.1.0 removed the process-wide app::instance() singleton in
  favour of constructible cvc::app contexts whose APIs take an
  app& argument. Rather than thread an app& through the entire Qt-based
  UI, volrover owns a single cvc::app via volrover_app_instance() and
  this header provides convenience wrappers that auto-supply it.

  This is the only volrover-side shim against libcvc; previous per-type
  headers under inc/CVC/ have been folded in here. Long term, library
  modules should migrate to the libcvc API directly and drop this header.
*/

#ifndef __VOLROVER_CVC_COMPAT_H__
#define __VOLROVER_CVC_COMPAT_H__

// Suppress libcvc's PascalCase typedefs for ThreadInfo/ThreadFeedback/
// ScopedLock — we define our own subclasses below that auto-supply the
// volrover-owned cvc::app to the new (app& ctx, ...) constructors.
#ifndef CVC_COMPAT_APP_DEFINED
#define CVC_COMPAT_APP_DEFINED
#endif

#include <cvc/namespace.h>
#include <cvc/types.h>
#include <cvc/app.h>
#include <cvc/state.h>
#include <cvc/state_object.h>
#include <cvc/bounding_box.h>
#include <cvc/dimension.h>
#include <cvc/exception.h>
#include <cvc/hdf5_utils.h>

namespace CVC_NAMESPACE
{
  typedef app App;

  // Volrover legacy code lives inside `namespace cvc { ... }` and refers
  // to PascalCase aliases unqualified. libcvc provides `BoundingBox`,
  // `Dimension`, etc. inside `namespace CVC` but not inside `cvc`, so
  // mirror them here.
  typedef bounding_box       BoundingBox;
  typedef index_bounding_box IndexBoundingBox;
  typedef dimension          Dimension;
  typedef exception          Exception;

  // The volrover-owned cvc::app context. Defined in
  // src/volrover_app_instance.cpp and linked into the VolumeRover2
  // executable; legacy static libraries that reference cvcapp/cvcstate
  // resolve this symbol at final link time.
  app &volrover_app_instance();

  // libcvc 3.1.0 made thread_info / thread_feedback / scoped_lock take a
  // leading app& argument. Wrap them so legacy volrover code can keep
  // constructing them with just the (optional) info string.
  struct ThreadInfo : public app::thread_info {
    ThreadInfo() : app::thread_info(volrover_app_instance()) {}
    explicit ThreadInfo(const std::string &info)
      : app::thread_info(volrover_app_instance(), info) {}
    ThreadInfo(app &ctx, const std::string &info = "running")
      : app::thread_info(ctx, info) {}
  };

  struct ThreadFeedback : public app::thread_feedback {
    ThreadFeedback() : app::thread_feedback(volrover_app_instance()) {}
    explicit ThreadFeedback(const std::string &key)
      : app::thread_feedback(volrover_app_instance(), key) {}
    ThreadFeedback(app &ctx, const std::string &key = "")
      : app::thread_feedback(ctx, key) {}
  };

  struct ScopedLock : public app::scoped_lock {
    explicit ScopedLock(const std::string &name,
                        const std::string &info = std::string())
      : app::scoped_lock(volrover_app_instance(), name, info) {}
    ScopedLock(app &ctx, const std::string &name,
               const std::string &info = std::string())
      : app::scoped_lock(ctx, name, info) {}
  };

  // Re-export cvc::hdf5_utils as CVC::HDF5_Utilities for legacy code,
  // and inject a few filename+objectpath overloads that previously
  // existed in volrover but were dropped in libcvc.
  namespace HDF5_Utilities
  {
    using namespace cvc::hdf5_utils;

    template<class T>
    inline void getAttribute(const std::string &hdf5_filename,
                             const std::string &hdf5_objname,
                             const std::string &attrname,
                             T &value)
    {
      auto f = cvc::hdf5_utils::getH5File(hdf5_filename);
      if (cvc::hdf5_utils::isGroup(volrover_app_instance(),
                                   hdf5_filename, hdf5_objname)) {
        auto g = cvc::hdf5_utils::getGroup(*f, hdf5_objname, false);
        cvc::hdf5_utils::getAttribute(*g, attrname, value);
      } else {
        auto d = cvc::hdf5_utils::getDataSet(*f, hdf5_objname, false);
        cvc::hdf5_utils::getAttribute(*d, attrname, value);
      }
    }

    template<class T>
    inline void setAttribute(const std::string &hdf5_filename,
                             const std::string &hdf5_objname,
                             const std::string &attrname,
                             const T &value)
    {
      auto f = cvc::hdf5_utils::getH5File(hdf5_filename);
      if (cvc::hdf5_utils::isGroup(volrover_app_instance(),
                                   hdf5_filename, hdf5_objname)) {
        auto g = cvc::hdf5_utils::getGroup(*f, hdf5_objname, false);
        cvc::hdf5_utils::setAttribute(*g, attrname, value);
      } else {
        auto d = cvc::hdf5_utils::getDataSet(*f, hdf5_objname, false);
        cvc::hdf5_utils::setAttribute(*d, attrname, value);
      }
    }

    inline void setAttribute(const std::string &hdf5_filename,
                             const std::string &hdf5_objname,
                             const std::string &attrname,
                             const char *value)
    {
      setAttribute<std::string>(hdf5_filename, hdf5_objname,
                                attrname, std::string(value));
    }

    inline std::vector<std::string>
    getChildObjects(const std::string &hdf5_filename,
                    const std::string &hdf5_objname)
    {
      return cvc::hdf5_utils::getChildObjects(volrover_app_instance(),
                                              hdf5_filename, hdf5_objname);
    }

    inline bool objectExists(const std::string &hdf5_filename,
                             const std::string &hdf5_objname)
    {
      return cvc::hdf5_utils::objectExists(volrover_app_instance(),
                                           hdf5_filename, hdf5_objname);
    }
  } // namespace HDF5_Utilities
} // namespace CVC_NAMESPACE

// Shorthand to access the volrover-owned cvc::app and its state from
// anywhere. These mirror the historic cvcapp / cvcstate macros.
#ifndef cvcapp
#define cvcapp (CVC_NAMESPACE::volrover_app_instance())
#endif
#ifndef cvcstate
#define cvcstate CVC_NAMESPACE::state::instance(CVC_NAMESPACE::volrover_app_instance())
#endif

#endif // __VOLROVER_CVC_COMPAT_H__
