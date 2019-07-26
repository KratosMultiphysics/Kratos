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

EXPORT void __stdcall RecreateProcessedMesh(ModelPartWrapper *instance) {
    instance->recreateProcessedMesh();
}

EXPORT ModelPartWrapper *__stdcall CreateSubmodelPart(ModelPartWrapper *instance, char *name) {
    return instance->createSubmodelPart(name);
}

EXPORT void __stdcall CreateNewNode(ModelPartWrapper *instance, int id, double x, double y, double z) {
    instance->createNewNode(id, x, y, z);
}

EXPORT void __stdcall CreateNewElement(ModelPartWrapper *instance, char *name, int id, int *nodeIds) {
    instance->createNewElement(name, id, nodeIds);
}

EXPORT void __stdcall CreateNew2dCondition(ModelPartWrapper *instance, char *name, int id, int *nodeIds) {
    instance->createNew2dCondition(name, id, nodeIds);
}

EXPORT void __stdcall RemoveNode(ModelPartWrapper *instance, int id) {
    instance->removeNode(id);
}

EXPORT void __stdcall RemoveElement(ModelPartWrapper *instance, int id) {
    instance->removeElement(id);
}

EXPORT void __stdcall RemoveCondition(ModelPartWrapper *instance, int id) {
    instance->removeCondition(id);
}

EXPORT void __stdcall AddNodes(ModelPartWrapper *instance, int *nodeIds, int nodeCount) {
    instance->addNodes(nodeIds, nodeCount);
}

EXPORT void __stdcall AddElements(ModelPartWrapper *instance, int *elementIds, int elementCount) {
    instance->addElements(elementIds, elementCount);
}

EXPORT void __stdcall AddConditions(ModelPartWrapper *instance, int *conditionIds, int conditionCount) {
    instance->addConditions(conditionIds, conditionCount);
}

EXPORT int __stdcall GetMaxElementId(ModelPartWrapper *instance) {
    return instance->getMaxElementId();
}

EXPORT int __stdcall GetMaxNodeId(ModelPartWrapper *instance) {
    return instance->getMaxNodeId();
}

EXPORT NodeType *__stdcall GetNode(ModelPartWrapper *instance, int id) {
    return instance->getNode(id);
}

EXPORT NodeType **__stdcall GetNodes(ModelPartWrapper *instance) {
    return instance->getNodes();
}

EXPORT int __stdcall GetNuberOfNodes(ModelPartWrapper *instance) {
    return instance->getNumberOfNodes();
}

EXPORT ElementType *__stdcall GetElement(ModelPartWrapper *instance, int id) {
    return instance->getElement(id);
}

EXPORT ElementType **__stdcall GetElements(ModelPartWrapper *instance) {
    return instance->getElements();
}

EXPORT int __stdcall GetNumberOfElements(ModelPartWrapper *instance) {
    return instance->getNumberOfElements();
}

EXPORT ConditionType *__stdcall GetCondition(ModelPartWrapper *instance, int id) {
    return instance->getCondition(id);
}

EXPORT ConditionType **__stdcall GetConditions(ModelPartWrapper *instance) {
    return instance->getConditions();
}

EXPORT int __stdcall GetNumberOfConditions(ModelPartWrapper *instance) {
    return instance->getNumberOfConditions();
}

EXPORT void __stdcall DisposeModelPartWrapper(ModelPartWrapper *instance) {
    delete instance;
}

}
