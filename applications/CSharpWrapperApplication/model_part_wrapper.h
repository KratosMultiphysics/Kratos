
#ifndef KRATOSMULTIPHYSICS_MODEL_PART_WRAPPER_H
#define KRATOSMULTIPHYSICS_MODEL_PART_WRAPPER_H

#include "includes/model_part.h"
#include "id_translator.h"
#include "mesh_converter.h"
#include "includes/variables.h"
#include <vector>


namespace CSharpKratosWrapper {

    using NodeType = Kratos::Node<3>;
    using ElementType = Kratos::Element;
    using ConditionType = Kratos::Condition;
    using ModelPart = Kratos::ModelPart;

    class ModelPartWrapper {

    public:

        ModelPartWrapper(Kratos::ModelPart &mModelPart, std::vector<NodeType::Pointer> &fixedNodes)
                : mModelPart(mModelPart), mFixedNodes(fixedNodes), pmParentModelPart(NULL) {
            initialize();
        }

        ModelPartWrapper(Kratos::ModelPart &mModelPart, std::vector<NodeType::Pointer> &fixedNodes,
                         ModelPartWrapper *parent)
                : mModelPart(mModelPart), mFixedNodes(fixedNodes), pmParentModelPart(parent) {
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

        void recreateProcessedMesh();

        ModelPartWrapper *createSubmodelPart(char *name);

        NodeType *createNewNode(int id, double x, double y, double z);

        ElementType *createNewElement(char *name, int id, int *nodeIds);

        ConditionType *createNew2dCondition(char *name, int id, int *nodeIds);

        void removeNode(int id);

        void removeElement(int id);

        void removeCondition(int id);

        void addNodes(int *nodeIds, int nodeCount);

        void addElements(int *elementIds, int elementCount);

        void addConditions(int *conditionIds, int elementCount);

        int getMaxElementId();

        int getMaxNodeId();

        NodeType *getNode(int id);

        NodeType **getNodes();

        int getNumberOfNodes();

        ElementType *getElement(int id);

        ElementType **getElements();

        int getNumberOfElements();

        ConditionType *getCondition(int id);

        ConditionType **getConditions();

        int getNumberOfConditions();

        IdTranslator *getIdTranslator();

    protected:

        void updateMaxElementId(int maxId);

        void updateMaxNodeId(int maxId);

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
        int mMaxElementId;
        int mMaxNodeId;
        bool mInitialized;

        void initialize();

        void saveNodes(MeshConverter &meshConverter);

        void saveTriangles(MeshConverter &meshConverter);

        void deleteSkin();
    };
};


#endif //KRATOSMULTIPHYSICS_MODEL_PART_WRAPPER_H
