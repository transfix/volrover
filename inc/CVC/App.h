/*
  Compatibility header: forwards to libcvc app class.
  CamelCase aliases (App, ThreadInfo, etc.) for legacy volrover code.
*/

#ifndef __CVC_COMPAT_APP_H__
#define __CVC_COMPAT_APP_H__

#include <cvc/app.h>
#include <CVC/Namespace.h>
#include <CVC/Types.h>

namespace CVC_NAMESPACE
{
  typedef app            App;
  typedef app::thread_info     ThreadInfo;
  typedef app::thread_feedback ThreadFeedback;
  typedef app::scoped_lock     ScopedLock;
}

#endif
