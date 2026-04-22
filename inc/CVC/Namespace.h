#ifndef __CVC_COMPAT_NAMESPACE_H__
#define __CVC_COMPAT_NAMESPACE_H__

#ifndef CVC_NAMESPACE
#define CVC_NAMESPACE cvc
#endif

// Forward-declare the canonical namespace and make CVC a real nested
// namespace that re-exports cvc. (Using `namespace CVC = cvc;` as a
// namespace-alias would conflict with libcvc's own declaration of CVC
// as a real nested namespace, and would also prevent consumer code
// from reopening `namespace CVC { ... }` to add compat typedefs.)
namespace CVC_NAMESPACE {}
namespace CVC { using namespace CVC_NAMESPACE; }

#endif
