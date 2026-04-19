/*
  Compatibility header: forwards to libcvc bounding_box type.
*/

#ifndef __CVC_COMPAT_BOUNDINGBOX_H__
#define __CVC_COMPAT_BOUNDINGBOX_H__

#include <cvc/bounding_box.h>
#include <CVC/Namespace.h>

namespace CVC_NAMESPACE
{
  template <typename T>
  using GenericBoundingBox = generic_bounding_box<T>;

  typedef bounding_box       BoundingBox;       // object space
  typedef index_bounding_box IndexBoundingBox;   // image space

  typedef invalid_bounding_box        InvalidBoundingBox;
  typedef invalid_bounding_box_string InvalidBoundingBoxString;
}

#endif
