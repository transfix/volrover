/*
  Compatibility header: forwards to libcvc state_object class.
*/

#ifndef __CVC_COMPAT_STATEOBJECT_H__
#define __CVC_COMPAT_STATEOBJECT_H__

#include <cvc/state_object.h>
#include <CVC/Namespace.h>

namespace CVC_NAMESPACE
{
  template <class This>
  using StateObject = state_object<This>;
}

#endif
