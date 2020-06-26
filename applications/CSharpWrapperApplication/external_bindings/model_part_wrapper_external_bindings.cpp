#define EXPORT __declspec(dllexport)

#include "custom_includes/model_part_wrapper.h"

using namespace CSharpKratosWrapper;

extern "C" {
#if defined(KRATOS_COMPILED_IN_WINDOWS)

EXPORT float *__stdcall ModelPartWrapper_GetXCoordinates(ModelPartWrapper *instance) {
    return instance->getXCoordinates();
}

EXPORT float *__stdcall ModelPartWrapper_GetYCoordinates(ModelPartWrapper *instance) {
    return instance->getYCoordinates();
}

EXPORT float *__stdcall ModelPartWrapper_GetZCoordinates(ModelPartWrapper *instance) {
    return instance->getZCoordinates();
}

EXPORT int __stdcall ModelPartWrapper_GetNodesCount(ModelPartWrapper *instance) {
    return instance->getNodesCount();
}

EXPORT int *__stdcall ModelPartWrapper_GetTriangles(ModelPartWrapper *instance) {
    return instance->getTriangles();
}

EXPORT int __stdcall ModelPartWrapper_GetTrianglesCount(ModelPartWrapper *instance) {
    return instance->getTrianglesCount();
}

EXPORT void __stdcall ModelPartWrapper_EnableSurfaceStressResults(ModelPartWrapper *instance) {
    instance->enableSurfaceStressResults();
}

EXPORT float *__stdcall ModelPartWrapper_GetSurfaceStress(ModelPartWrapper *instance) {
    return instance->getSurfaceStress();
}

EXPORT void
__stdcall ModelPartWrapper_UpdateNodePos(ModelPartWrapper *instance, int nodeId, float x, float y, float z) {
    instance->updateNodePos(nodeId, x, y, z);
}

EXPORT bool __stdcall ModelPartWrapper_HasSubmodelPart(ModelPartWrapper *instance, char *name) {
    return instance->hasSubmodelPart(name);
}

EXPORT ModelPartWrapper *__stdcall ModelPartWrapper_GetSubmodelPart(ModelPartWrapper *instance, char *name) {
    return instance->getSubmodelPart(name);
}

EXPORT void __stdcall ModelPartWrapper_RetrieveResults(ModelPartWrapper *instance) {
    instance->retrieveResults();
}

EXPORT void __stdcall ModelPartWrapper_RecreateProcessedMesh(ModelPartWrapper *instance) {
    instance->recreateProcessedMesh();
}

EXPORT ModelPartWrapper *__stdcall ModelPartWrapper_CreateSubmodelPart(ModelPartWrapper *instance, char *name) {
    return instance->createSubmodelPart(name);
}

EXPORT NodeType *
__stdcall ModelPartWrapper_CreateNewNode(ModelPartWrapper *instance, int id, double x, double y, double z) {
    return instance->createNewNode(id, x, y, z);
}

EXPORT ElementType *
__stdcall ModelPartWrapper_CreateNewElement(ModelPartWrapper *instance, char *name, int id, int *nodeIds) {
    return instance->createNewElement(name, id, nodeIds);
}

EXPORT ConditionType *
__stdcall ModelPartWrapper_CreateNew2dCondition(ModelPartWrapper *instance, char *name, int id, int *nodeIds) {
    return instance->createNew2dCondition(name, id, nodeIds);
}

EXPORT void __stdcall ModelPartWrapper_RemoveNode(ModelPartWrapper *instance, int id) {
    instance->removeNode(id);
}

EXPORT void __stdcall ModelPartWrapper_RemoveElement(ModelPartWrapper *instance, int id) {
    instance->removeElement(id);
}

EXPORT void __stdcall ModelPartWrapper_RemoveCondition(ModelPartWrapper *instance, int id) {
    instance->removeCondition(id);
}

EXPORT void __stdcall ModelPartWrapper_AddNodes(ModelPartWrapper *instance, int *nodeIds, int nodeCount) {
    instance->addNodes(nodeIds, nodeCount);
}

EXPORT void __stdcall ModelPartWrapper_AddElements(ModelPartWrapper *instance, int *elementIds, int elementCount) {
    instance->addElements(elementIds, elementCount);
}

EXPORT void
__stdcall ModelPartWrapper_AddConditions(ModelPartWrapper *instance, int *conditionIds, int conditionCount) {
    instance->addConditions(conditionIds, conditionCount);
}

EXPORT int __stdcall ModelPartWrapper_GetMaxElementId(ModelPartWrapper *instance) {
    return instance->getMaxElementId();
}

EXPORT int __stdcall ModelPartWrapper_GetMaxNodeId(ModelPartWrapper *instance) {
    return instance->getMaxNodeId();
}

EXPORT NodeType *__stdcall ModelPartWrapper_GetNode(ModelPartWrapper *instance, int id) {
    return instance->getNode(id);
}

EXPORT NodeType **__stdcall ModelPartWrapper_GetNodes(ModelPartWrapper *instance) {
    return instance->getNodes();
}

EXPORT int __stdcall ModelPartWrapper_GetNumberOfNodes(ModelPartWrapper *instance) {
    return instance->getNumberOfNodes();
}

EXPORT ElementType *__stdcall ModelPartWrapper_GetElement(ModelPartWrapper *instance, int id) {
    return instance->getElement(id);
}

EXPORT ElementType **__stdcall ModelPartWrapper_GetElements(ModelPartWrapper *instance) {
    return instance->getElements();
}

EXPORT int __stdcall ModelPartWrapper_GetNumberOfElements(ModelPartWrapper *instance) {
    return instance->getNumberOfElements();
}

EXPORT ConditionType *__stdcall ModelPartWrapper_GetCondition(ModelPartWrapper *instance, int id) {
    return instance->getCondition(id);
}

EXPORT ConditionType **__stdcall ModelPartWrapper_GetConditions(ModelPartWrapper *instance) {
    return instance->getConditions();
}

EXPORT int __stdcall ModelPartWrapper_GetNumberOfConditions(ModelPartWrapper *instance) {
    return instance->getNumberOfConditions();
}

EXPORT IdTranslator *__stdcall ModelPartWrapper_GetIdTranslator(ModelPartWrapper *instance) {
    return instance->getIdTranslator();
}

EXPORT void __stdcall ModelPartWrapper_DisposeModelPartWrapper(ModelPartWrapper *instance) {
    delete instance;
}

EXPORT bool
__stdcall ModelPartWrapper_HasNodalVariable1d(ModelPartWrapper *instance, Kratos::Variable<double> *variable) {
    return instance->hasNodalVariable1d(*variable);
}

EXPORT bool __stdcall ModelPartWrapper_HasNodalVariable3d(ModelPartWrapper *instance,
                                                          Kratos::Variable<Kratos::array_1d<double, 3>> *variable) {
    return instance->hasNodalVariable3d(*variable);
}

EXPORT double *
__stdcall ModelPartWrapper_GetNodalVariables1d(ModelPartWrapper *instance, Kratos::Variable<double> *variable) {
    return instance->getNodalVariable1d(*variable);
}

EXPORT double *__stdcall ModelPartWrapper_GetNodalVariableComponents(ModelPartWrapper *instance,
                                                                     Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > *variable) {
    return instance->getNodalVariableComponent(*variable);
}

EXPORT double *__stdcall ModelPartWrapper_GetNodalVariables3d(ModelPartWrapper *instance,
                                                              Kratos::Variable<Kratos::array_1d<double, 3>> *variable) {
    return instance->getNodalVariable3d(*variable);
}
#endif
}