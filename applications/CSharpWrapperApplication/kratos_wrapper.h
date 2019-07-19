//
// Created by huber on 2019-07-16.
//

#ifndef KRATOSMULTIPHYSICS_KRATOS_WRAPPER_H
#define KRATOSMULTIPHYSICS_KRATOS_WRAPPER_H


#include "kratos_internals.h"
#include "mesh_converter.h"
#include "id_translator.h"
#include "includes/node.h"

namespace CSharpKratosWrapper {

    using NodeType = Kratos::Node<3>;
    using ModelPart = Kratos::ModelPart;

    class KratosWrapper {

    public:
        void init(const char *MDPAFilePath, const char *JSONFilePath = NULL);

        void initWithSettings(const char *JSONFilePath = NULL);

        void updateNodePos(const int nodeId, const float x, const float y, const float z);

        void calculate();

        float *getXCoordinates();

        float *getYCoordinates();

        float *getZCoordinates();

        int getNodesCount();

        int *getTriangles();

        int getTrianglesCount();

        void enableSurfaceReactions();

        float *getSurfaceReactions();

    private:
        int *pmKratosNodeIds;
        int *pmUnityNodeIds;
        int *pmTriangles;
        int mTrianglesCount;
        int mNodesCount;
        float *pmXCoordinates;
        float *pmYCoordinates;
        float *pmZCoordinates;
        bool mReactionResultsEnabled;
        float *pmSurfaceStress;
        KratosInternals mKratosInternals;
        IdTranslator idTranslator;
        std::vector<NodeType::Pointer> mFixedNodes;

        void saveTriangles(MeshConverter &meshConverter);

        void saveNodes(MeshConverter &meshConverter);

        void retrieveNodesData();

        void freeNodes();

    };
}

#endif //KRATOSMULTIPHYSICS_KRATOS_WRAPPER_H
