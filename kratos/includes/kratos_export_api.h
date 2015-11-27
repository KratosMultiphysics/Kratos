#ifndef KRATOS_EXPORT_API_H
#define KRATOS_EXPORT_API_H

#undef KRATOS_API_EXPORT
#undef KRATOS_API_IMPORT
#ifdef _WIN32
  #define KRATOS_API_EXPORT __declspec(dllexport)
  #define KRATOS_API_IMPORT __declspec(dllimport)
#else
  #define KRATOS_API_EXPORT __attribute__((visibility("default")))
  #define KRATOS_API_IMPORT __attribute__((visibility("default")))
#endif

// This fixes MSVC not expanding __VA_ARGS__ as defined in the C99 standard
#define KRATOS_EXPAND(A) A

// This expands the API call to either import or export based on the
// number of the arguments in API() call.
#define KRATOS_API_CALL(x,T1,T2,T3,...) T3
#define KRATOS_API(...) \
  KRATOS_EXPAND(KRATOS_API_CALL(,##__VA_ARGS__,KRATOS_API_EXPORT,KRATOS_API_IMPORT))
#define KRATOS_NO_EXPORT(...)

#endif
