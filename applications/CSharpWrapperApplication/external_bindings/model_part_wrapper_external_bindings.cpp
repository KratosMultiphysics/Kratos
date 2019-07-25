#define EXPORT __declspec(dllexport)

#include "model_part_wrapper.h"

using namespace CSharpKratosWrapper;

extern "C" {
EXPORT float *__stdcall GetXCoordinates(ModelPartWrapper *instance) {
    return instance->getXCoordinates();
}

EXPORT float *__stdcall GetYCoordinates(ModelPartWrapper *instance) {
    return instance->getYCoordinates();
}

EXPORT float *__stdcall GetZCoordinates(ModelPartWrapper *instance) {
    return instance->getZCoordinates();
}

EXPORT int __stdcall GetNodesCount(ModelPartWrapper *instance) {
    return instance->getNodesCount();
}

EXPORT int *__stdcall GetTriangles(ModelPartWrapper *instance) {
    return instance->getTriangles();
}

EXPORT int __stdcall GetTrianglesCount(ModelPartWrapper *instance) {
    return instance->getTrianglesCount();
}

EXPORT void __stdcall EnableSurfaceStressResults(ModelPartWrapper *instance) {
    instance->enableSurfaceStressResults();
}

EXPORT float *__stdcall GetSurfaceStress(ModelPartWrapper *instance) {
    return instance->getSurfaceStress();
}

EXPORT void __stdcall UpdateNodePos(ModelPartWrapper *instance, int nodeId, float x, float y, float z) {
    instance->updateNodePos(nodeId, x, y, z);
}

EXPORT bool __stdcall HasSubmodelPart(ModelPartWrapper *instance, char *name) {
    return instance->hasSubmodelPart(name);
}

EXPORT ModelPartWrapper *__stdcall GetSubmodelPart(ModelPartWrapper *instance, char *name) {
    return instance->getSubmodelPart(name);
}

EXPORT void __stdcall RetrieveResults(ModelPartWrapper *instance) {
    instance->retrieveResults();
}

EXPORT void __stdcall DisposeModelPartWrapper(ModelPartWrapper *instance) {
    delete instance;
}

}
