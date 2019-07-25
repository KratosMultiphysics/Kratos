
#ifndef KRATOSMULTIPHYSICS_KRATOS_WRAPPER_H
#define KRATOSMULTIPHYSICS_KRATOS_WRAPPER_H


#include "kratos_internals.h"
#include "includes/node.h"
#include "model_part_wrapper.h"

namespace CSharpKratosWrapper {

    using NodeType = Kratos::Node<3>;
    using ModelPart = Kratos::ModelPart;

    class KratosWrapper {

    public:

        ~KratosWrapper();

        void init(const char *MDPAFilePath, const char *JSONFilePath = NULL);

        void initWithSettings(const char *JSONFilePath = NULL);

        void calculate();

        ModelPartWrapper* getRootModelPartWrapper();

    private:
        KratosInternals mKratosInternals;
        std::vector<NodeType::Pointer> mFixedNodes;
        ModelPartWrapper* pmMainModelPartWrapper;

        void freeNodes();

    };
}

#endif //KRATOSMULTIPHYSICS_KRATOS_WRAPPER_H
