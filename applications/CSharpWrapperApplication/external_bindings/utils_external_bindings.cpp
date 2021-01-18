#define EXPORT __declspec(dllexport)

#include "includes/define.h"

extern "C" {
#if defined(KRATOS_COMPILED_IN_WINDOWS)

EXPORT void __stdcall Utils_DeleteDoubleArray(double *array) {
    delete[] array;
}
#endif
}
