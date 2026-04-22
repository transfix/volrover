/*
  Compatibility header: forwards to libcvc exception types
  and adds CVC-namespace aliases for volrover code.
*/

#ifndef __CVC_COMPAT_EXCEPTION_H__
#define __CVC_COMPAT_EXCEPTION_H__

#include <cvc/exception.h>
#include <CVC/Namespace.h>

// Re-export the CVC_DEF_EXCEPTION macro from libcvc (already defined)

#ifndef CVC_COMPAT_EXCEPTION_DEFINED
#define CVC_COMPAT_EXCEPTION_DEFINED
namespace CVC_NAMESPACE
{
  // Base exception class
  typedef exception Exception;

  // Exception types with CamelCase aliases
  typedef read_error                         ReadError;
  typedef write_error                        WriteError;
  typedef memory_allocation_error             MemoryAllocationError;
  typedef index_out_of_bounds                 IndexOutOfBounds;
  typedef volume_properties_mismatch          VolumePropertiesMismatch;
  typedef volume_cache_directory_file_error   VolumeCacheDirectoryFileError;
  typedef unsupported_exception               UnsupportedVolumeFileType;
  typedef unsupported_exception               UnsupportedGeometryFileType;

  // Exceptions that don't have libcvc equivalents — define locally
  CVC_DEF_EXCEPTION(SubVolumeOutOfBounds);
  CVC_DEF_EXCEPTION(NetworkError);
  CVC_DEF_EXCEPTION(XmlRpcServerTerminate);
}
#endif // CVC_COMPAT_EXCEPTION_DEFINED

#endif

