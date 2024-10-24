#include "undefineMacros.hpp"

#define MCGS_INTERNAL

#if _WIN32
    #define MCGS_EXPORT_SYMBOL __declspec(dllexport)
#else
    #define MCGS_EXPORT_SYMBOL __attribute__((visibility ("default")))
#endif
