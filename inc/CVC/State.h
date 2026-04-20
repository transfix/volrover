/*
  Compatibility header: forwards to libcvc state class.
*/

#ifndef __CVC_COMPAT_STATE_H__
#define __CVC_COMPAT_STATE_H__

#include <cvc/state.h>
#include <CVC/Namespace.h>

namespace CVC_NAMESPACE
{
  typedef state State;
}

//Shorthand to access the State object from anywhere
#ifndef cvcstate
#define cvcstate CVC_NAMESPACE::state::instance()
#endif

#endif
