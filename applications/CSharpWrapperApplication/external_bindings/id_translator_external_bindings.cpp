#define EXPORT __declspec(dllexport)

#include "id_translator.h"

using namespace CSharpKratosWrapper;

extern "C" {
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
}