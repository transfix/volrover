#ifndef __VOLMAGICK_EXCEPTIONS_H__
#define __VOLMAGICK_EXCEPTIONS_H__

#include <cvc/exception.h>
#include <boost/format.hpp>
#include <string>

namespace VolMagick
{
  typedef cvc::exception Exception;

#define VOLMAGICK_DEF_EXCEPTION(name) \
  class name : public VolMagick::Exception \
  { \
  public: \
    name () : _msg("VolMagick::"#name) {} \
    name (const std::string& msg) : \
      _msg(boost::str(boost::format("VolMagick::" #name " exception: %1%") % msg)) {} \
    virtual ~name() throw() {} \
    virtual const std::string& what_str() const throw() { return _msg; } \
  private: \
    std::string _msg; \
  }

  VOLMAGICK_DEF_EXCEPTION(ReadError);
  VOLMAGICK_DEF_EXCEPTION(WriteError);
  VOLMAGICK_DEF_EXCEPTION(MemoryAllocationError);
  VOLMAGICK_DEF_EXCEPTION(SubVolumeOutOfBounds);
  VOLMAGICK_DEF_EXCEPTION(UnsupportedVolumeFileType);
  VOLMAGICK_DEF_EXCEPTION(IndexOutOfBounds);
  VOLMAGICK_DEF_EXCEPTION(NullDimension);
  VOLMAGICK_DEF_EXCEPTION(VolumePropertiesMismatch);
  VOLMAGICK_DEF_EXCEPTION(VolumeCacheDirectoryFileError);
  VOLMAGICK_DEF_EXCEPTION(InvalidBoundingBox);
}

#endif
