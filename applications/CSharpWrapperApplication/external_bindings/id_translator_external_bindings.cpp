#define EXPORT __declspec(dllexport)

#include "custom_includes/id_translator.h"
#include "includes/define.h"

using namespace CSharpKratosWrapper;

extern "C" {
#if defined(KRATOS_COMPILED_IN_WINDOWS)

EXPORT bool __stdcall Translator_HasSurfaceId(IdTranslator *instance, int kratosId) {
    return instance->hasSurfaceId(kratosId);
}

EXPORT bool __stdcall Translator_HasKratosId(IdTranslator *instance, int surfaceId) {
    return instance->hasKratosId(surfaceId);
}

EXPORT int __stdcall Translator_GetSurfaceId(IdTranslator instance, int kratosId) {
    return instance.safeGetSurfaceId(kratosId);
}

EXPORT int __stdcall Translator_GetKratosId(IdTranslator *instance, int surfaceId) {
    return instance->safeGetKratosId(surfaceId);
}
#endif
}