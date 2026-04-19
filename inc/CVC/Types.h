/*
  Compatibility header: forwards to libcvc types and adds CVC-namespace aliases
  for volrover application framework types not in libcvc.
*/

#ifndef __CVC_COMPAT_TYPES_H__
#define __CVC_COMPAT_TYPES_H__

// Pull in all libcvc types (defines cvc::int64, cvc::data_type, etc.)
#include <cvc/types.h>
#include <CVC/Namespace.h>

// Add CamelCase aliases within the CVC_NAMESPACE for backward compat
namespace CVC_NAMESPACE
{
  // DataType enum
  typedef data_type DataType;

  static const unsigned int *DataTypeSizes = data_type_sizes;
  static const char **DataTypeStrings = data_type_strings;

  // LocaleBool
  typedef locale_bool LocaleBool;

  // App-framework types
  typedef signal                Signal;
  typedef map_change_signal     MapChangeSignal;
  typedef data_map              DataMap;
  typedef data_type_name_map    DataTypeNameMap;
  typedef data_type_enum_map    DataTypeEnumMap;
  typedef property_map          PropertyMap;
  typedef thread_ptr            ThreadPtr;
  typedef thread_map            ThreadMap;
  typedef thread_progress_map   ThreadProgressMap;
  typedef thread_key_map        ThreadKeyMap;
  typedef thread_info_map       ThreadInfoMap;
  typedef data_reader           DataReader;
  typedef data_reader_collection DataReaderCollection;
  typedef mutex_ptr             MutexPtr;
  typedef mutex_map_element     MutexMapElement;
  typedef mutex_map             MutexMap;
}

#endif
