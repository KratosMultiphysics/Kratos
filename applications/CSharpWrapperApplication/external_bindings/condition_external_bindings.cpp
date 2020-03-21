#define EXPORT __declspec(dllexport)

#include "includes/condition.h"
#include "includes/node.h"

using namespace std;
typedef Kratos::Condition ConditionType;
typedef Kratos::Node<3> NodeType;

extern "C" {
#if defined(KRATOS_COMPILED_IN_WINDOWS)

EXPORT int __stdcall Condition_Id(ConditionType *instance) {
    return instance->Id();
}

EXPORT NodeType **__stdcall Condition_Nodes(ConditionType *instance) {
    auto nodeVector = instance->pGetGeometry()->Points();
    auto **nodes = new NodeType *[nodeVector.size()];
    for (unsigned int i = 0; i < nodeVector.size(); i++) {
        nodes[i] = &nodeVector[i];
    }
    return nodes;
}

EXPORT double __stdcall Condition_GetVariable1d(ConditionType *instance, Kratos::Variable<double> *variable) {
    return instance->GetValue(*variable);
}

EXPORT double *
__stdcall Condition_GetVariable3d(ConditionType *instance, Kratos::Variable<Kratos::array_1d<double, 3>> *variable) {
    auto value = instance->GetValue(*variable);
    auto *result = new double[3];
    for (int i = 0; i < 3; i++) result[i] = value[i];
    return result;
}

EXPORT double __stdcall Condition_GetVariableComponent(ConditionType *instance,
                                                     Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > *variable) {
    return instance->GetValue(*variable);
}

EXPORT bool __stdcall Condition_HasVariable1d(ConditionType *instance, Kratos::Variable<double> *variable) {
    return instance->Has(*variable);
}

EXPORT bool
__stdcall Condition_HasVariable3d(ConditionType *instance, Kratos::Variable<Kratos::array_1d<double, 3>> *variable) {
    return instance->Has(*variable);
}

EXPORT bool __stdcall Condition_HasVariableComponent(ConditionType *instance,
                                                   Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > *variable) {
    return instance->Has(*variable);
}
#endif
}