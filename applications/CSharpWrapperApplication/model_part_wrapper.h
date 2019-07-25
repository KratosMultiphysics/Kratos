
#ifndef KRATOSMULTIPHYSICS_MODEL_PART_WRAPPER_H
#define KRATOSMULTIPHYSICS_MODEL_PART_WRAPPER_H

#include "includes/model_part.h"
#include "id_translator.h"
#include "mesh_converter.h"
#include "includes/variables.h"
#include <vector>


namespace CSharpKratosWrapper {

    using NodeType = Kratos::Node<3>;
    using ModelPart = Kratos::ModelPart;

    class ModelPartWrapper {

    public:

        ModelPartWrapper(Kratos::ModelPart &mModelPart, std::vector<NodeType::Pointer> &mFixedNodes)
                : mModelPart(mModelPart), mFixedNodes(mFixedNodes), pmParentModelPart(NULL) {
            initialize();
        }

        ModelPartWrapper(Kratos::ModelPart &mModelPart, std::vector<NodeType::Pointer> &mFixedNodes,
                         ModelPartWrapper *parent)
                : mModelPart(mModelPart), mFixedNodes(mFixedNodes), pmParentModelPart(parent) {
            initialize();
        }

        ~ModelPartWrapper();


        ModelPartWrapper *getSubmodelPart(char *name);

        bool hasSubmodelPart(char *name);

        float *getXCoordinates();

        float *getYCoordinates();

        float *getZCoordinates();

        int getNodesCount();

        int *getTriangles();

        int getTrianglesCount();

        void updateNodePos(const int nodeId, const float x, const float y, const float z);

        void retrieveResults();

        void enableSurfaceStressResults();

        float *getSurfaceStress();

        ModelPart &getKratosModelPart();

    protected:
        int getMaxId();
        void updateMaxId(int maxId);

    private:
        Kratos::ModelPart &mModelPart;
        ModelPartWrapper *pmParentModelPart;
        IdTranslator idTranslator;
        std::vector<NodeType::Pointer> &mFixedNodes;

        float *pmXCoordinates;
        float *pmYCoordinates;
        float *pmZCoordinates;
        int *pmTriangles;
        int mNodesCount;
        int mTrianglesCount;
        float *pmSurfaceStress;
        bool mStressResultsEnabled;
        int mMaxId;


        void initialize();

        void saveNodes(MeshConverter &meshConverter);

        void saveTriangles(MeshConverter &meshConverter);


    };
};


#endif //KRATOSMULTIPHYSICS_MODEL_PART_WRAPPER_H
