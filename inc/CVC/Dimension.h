/*
  Compatibility header: forwards to libcvc dimension type.
*/

#ifndef __CVC_COMPAT_DIMENSION_H__
#define __CVC_COMPAT_DIMENSION_H__

#include <cvc/dimension.h>
#include <CVC/Namespace.h>

namespace CVC_NAMESPACE
{
  typedef dimension Dimension;
  typedef null_dimension NullDimension;
  typedef invalid_dimension_string InvalidDimensionString;
}

#endif
