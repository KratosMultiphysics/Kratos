#define EXPORT __declspec(dllexport)

using namespace std;

#include "includes/element.h"
#include "includes/node.h"

typedef Kratos::Element ElementType;
typedef Kratos::Node<3> NodeType;

extern "C" {

EXPORT int __stdcall Element_Id(ElementType *instnace) {
    return instnace->Id();
}

EXPORT NodeType **__stdcall Element_Nodes(ElementType *instance) {
    auto nodeVector = instance->pGetGeometry()->Points();
    auto **nodes = new NodeType *[nodeVector.size()];
    for (int i = 0; i < nodeVector.size(); i++) {
        nodes[i] = &nodeVector[i];
    }
    return nodes;
}
}