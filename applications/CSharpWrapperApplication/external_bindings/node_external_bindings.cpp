#define EXPORT __declspec(dllexport)

#include "includes/node.h"

using namespace std;
typedef Kratos::Node<3> NodeType;

extern "C" {
#if defined(KRATOS_COMPILED_IN_WINDOWS)
EXPORT int __stdcall Node_Id(NodeType *instance) {
    return instance->Id();
}


EXPORT double __stdcall Node_X(NodeType *instance) {
    return instance->X();
}

EXPORT double __stdcall Node_Y(NodeType *instance) {
    return instance->Y();
}

EXPORT double __stdcall Node_Z(NodeType *instance) {
    return instance->Z();
}

EXPORT double __stdcall Node_X0(NodeType *instance) {
    return instance->X0();
}

EXPORT double __stdcall Node_Y0(NodeType *instance) {
    return instance->Y0();
}

EXPORT double __stdcall Node_Z0(NodeType *instance) {
    return instance->Z0();
}

EXPORT double __stdcall Node_GetVariable1d(NodeType *instance, Kratos::Variable<double> *variable) {
    return instance->FastGetSolutionStepValue(*variable);
}

EXPORT double __stdcall Node_GetVariableComponent1d(NodeType *instance,
                                                    Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > *variable) {
    return instance->FastGetSolutionStepValue(*variable);
}

EXPORT double *
__stdcall Node_GetVariable3d(NodeType *instance, Kratos::Variable<Kratos::array_1d<double, 3>> *variable) {
    auto value = instance->FastGetSolutionStepValue(*variable);
    auto *result = new double[3];
    for (int i = 0; i < 3; i++) result[i] = value[i];
    return result;
}
#endif
}