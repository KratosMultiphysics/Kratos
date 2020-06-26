//    _____  _____ _                  __          __                                                 _ _           _   _
//   / ____|/ ____| |                 \ \        / /                               /\               | (_)         | | (_)
//  | |    | (___ | |__   __ _ _ __ _ _\ \  /\  / / __ __ _ _ __  _ __   ___ _ __ /  \   _ __  _ __ | |_  ___ __ _| |_ _  ___  _ __
//  | |     \___ \| '_ \ / _` | '__| '_ \ \/  \/ / '__/ _` | '_ \| '_ \ / _ \ '__/ /\ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//  | |____ ____) | | | | (_| | |  | |_) \  /\  /| | | (_| | |_) | |_) |  __/ | / ____ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//   \_____|_____/|_| |_|\__,_|_|  | .__/ \/  \/ |_|  \__,_| .__/| .__/ \___|_|/_/    \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                 | |                     | |   | |                    | |   | |
//                                 |_|                     |_|   |_|                    |_|   |_|
//
//
//  License: BSD License
//   license: CSharpWrapperApplication/license.txt
//
//  Main authors:    Hubert Balcerzak

#ifndef KRATOSMULTIPHYSICS_MODEL_PART_WRAPPER_H
#define KRATOSMULTIPHYSICS_MODEL_PART_WRAPPER_H

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/model_part.h"
#include "id_translator.h"
#include "mesh_converter.h"
#include "includes/variables.h"


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

        /**
         * Fixes and updates DISPLACEMENT variable of a node, so that final position is as given. X0 + DISPLACEMENT_X = x
         * @param nodeId Surface id of the node to update
         * @param x X coordinate
         * @param y Y coordinate
         * @param z Z coordinate
         */
        void updateNodePos(const int nodeId, const float x, const float y, const float z);

        /**
         * Extracts new node positions and stress values
         */
        void retrieveResults();

        /**
         * Enables calculation of surface stress values. Call `retrieveResults` after each simulation to recalculate stress.
         */
        void enableSurfaceStressResults();

        float *getSurfaceStress();

        ModelPart &getKratosModelPart();

        /**
         * Recreates surface mesh after changes to its structure, like creation of new nodes and elements
         */
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

        double *getNodalVariable1d(Kratos::Variable<double> &variable);

        double *getNodalVariableComponent(
                Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > &variable);

        double *getNodalVariable3d(Kratos::Variable<Kratos::array_1d<double, 3>> &variable);

        bool hasNodalVariable1d(Kratos::Variable<double> &variable);

        bool hasNodalVariable3d(Kratos::Variable<Kratos::array_1d<double, 3>> &variable);

    protected:

        void updateMaxElementId(int maxId);

        void updateMaxNodeId(int maxId);

    private:
        Kratos::ModelPart &mModelPart;
        std::vector<NodeType::Pointer> &mFixedNodes;
        ModelPartWrapper *pmParentModelPart;
        IdTranslator idTranslator;

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
