#define EXPORT __declspec(dllexport)

#include "includes/node.h"

using namespace std;
typedef Kratos::Node<3> NodeType;

extern "C" {

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

}