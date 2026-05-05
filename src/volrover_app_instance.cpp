/*
  volrover-owned cvc::app singleton.

  libcvc 3.1.0 removed the process-wide app::instance() / state::instance()
  singletons. Volrover is a legacy application with ~1300 references to
  the cvcapp / cvcstate convenience macros; threading a cvc::app& through
  the entire (Qt-based) UI is impractical.

  Instead, volrover owns a single function-local static cvc::app and the
  cvcapp / cvcstate macros (defined in inc/cvc_compat.h) route through
  this accessor.
*/

#include <cvc_compat.h>

#include <cvc/app.h>

namespace CVC_NAMESPACE {

app &volrover_app_instance()
{
  static app instance_;
  return instance_;
}

}  // namespace CVC_NAMESPACE
