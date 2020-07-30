#define EXPORT __declspec(dllexport)


#include "includes/element.h"
#include "includes/node.h"

typedef Kratos::Element ElementType;
typedef Kratos::Node<3> NodeType;

extern "C" {
#if defined(KRATOS_COMPILED_IN_WINDOWS)

EXPORT int __stdcall Element_Id(ElementType *instnace) {
    return instnace->Id();
}

EXPORT NodeType **__stdcall Element_Nodes(ElementType *instance) {
    auto nodeVector = instance->pGetGeometry()->Points();
    auto **nodes = new NodeType *[nodeVector.size()];
    for (unsigned int i = 0; i < nodeVector.size(); i++) {
        nodes[i] = &nodeVector[i];
    }
    return nodes;
}

EXPORT double __stdcall Element_GetVariable1d(ElementType *instance, Kratos::Variable<double> *variable) {
    return instance->GetValue(*variable);
}

EXPORT double *
__stdcall Element_GetVariable3d(ElementType *instance, Kratos::Variable<Kratos::array_1d<double, 3>> *variable) {
    auto value = instance->GetValue(*variable);
    auto *result = new double[3];
    for (int i = 0; i < 3; i++) result[i] = value[i];
    return result;
}

EXPORT double __stdcall Element_GetVariableComponent(ElementType *instance,
                                                     Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > *variable) {
    return instance->GetValue(*variable);
}

EXPORT bool __stdcall Element_HasVariable1d(ElementType *instance, Kratos::Variable<double> *variable) {
    return instance->Has(*variable);
}

EXPORT bool
__stdcall Element_HasVariable3d(ElementType *instance, Kratos::Variable<Kratos::array_1d<double, 3>> *variable) {
    return instance->Has(*variable);
}

EXPORT bool __stdcall Element_HasVariableComponent(ElementType *instance,
                                                   Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > *variable) {
    return instance->Has(*variable);
}
#endif
}