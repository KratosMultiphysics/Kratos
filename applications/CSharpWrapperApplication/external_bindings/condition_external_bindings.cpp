#define EXPORT __declspec(dllexport)

#include "includes/condition.h"
#include "includes/node.h"

using namespace std;
typedef Kratos::Condition ConditionType;
typedef Kratos::Node<3> NodeType;

extern "C" {

EXPORT int __stdcall Condition_Id(ConditionType *instance) {
    return instance->Id();
}

EXPORT NodeType **__stdcall Condition_Nodes(ConditionType *instance) {
    auto nodeVector = instance->pGetGeometry()->Points();
    auto **nodes = new NodeType *[nodeVector.size()];
    for (int i = 0; i < nodeVector.size(); i++) {
        nodes[i] = &nodeVector[i];
    }
    return nodes;
}
}