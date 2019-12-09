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

// System includes
#include <stdlib.h>

// External includes

// Project includes
#include <structural_mechanics_application_variables.h>
#include <processes/tetrahedral_mesh_orientation_check.h>
#include "custom_includes/model_part_wrapper.h"

#define SKIN_SUBMODEL_PART_NAME "CSharpWrapper_skin"

using namespace CSharpKratosWrapper;

bool ModelPartWrapper::hasSubmodelPart(char *name) {
    return mModelPart.HasSubModelPart(name);
}

float *ModelPartWrapper::getXCoordinates() {
    return pmXCoordinates;
}

float *ModelPartWrapper::getYCoordinates() {
    return pmYCoordinates;
}

float *ModelPartWrapper::getZCoordinates() {
    return pmZCoordinates;
}

int ModelPartWrapper::getNodesCount() {
    return mNodesCount;
}

int *ModelPartWrapper::getTriangles() {
    auto &r_conditions_array = mModelPart.GetSubModelPart(SKIN_SUBMODEL_PART_NAME).Conditions();
    const auto it_condition_begin = r_conditions_array.begin();
    int *triangles = new int[mTrianglesCount * 3];
    for (int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_condition = it_condition_begin + i;
        auto geometry = it_condition->GetGeometry();

        for (int j = 0; j < 3; j++) {
            triangles[3 * i + j] = idTranslator.getSurfaceId(geometry.GetPoint(j).Id());
        }
    }
    return triangles;
}

int ModelPartWrapper::getTrianglesCount() {
    return mTrianglesCount;
}

void ModelPartWrapper::updateNodePos(const int nodeId, const float x, const float y, const float z) {

    // Get node
    NodeType::Pointer p_node = mModelPart.pGetNode(idTranslator.getKratosId(nodeId));

    // Fix nodes
    p_node->Fix(Kratos::DISPLACEMENT_X);
    p_node->Fix(Kratos::DISPLACEMENT_Y);
    p_node->Fix(Kratos::DISPLACEMENT_Z);

    // Arrays
    Kratos::array_1d<double, 3> &r_displacement = p_node->FastGetSolutionStepValue(Kratos::DISPLACEMENT);
    const Kratos::array_1d<double, 3> &r_initial_coordinates = p_node->GetInitialPosition().Coordinates();

    // Update coordinates
    Kratos::array_1d<double, 3> &r_coordinates = p_node->Coordinates();
    r_coordinates[0] = x;
    r_coordinates[1] = y;
    r_coordinates[2] = z;

    // Update displacement
    r_displacement[0] = x - r_initial_coordinates[0];
    r_displacement[1] = y - r_initial_coordinates[1];
    r_displacement[2] = z - r_initial_coordinates[2];

    // Adding to the database of fixed nodes
    mFixedNodes.push_back(p_node);
}

void ModelPartWrapper::initialize() {
    mStressResultsEnabled = false;

    if (mModelPart.NumberOfElements() > 0) {
        mMaxElementId = mModelPart.Elements().back().Id();
        mMaxNodeId = mModelPart.Nodes().back().Id();
        MeshConverter meshConverter;
        meshConverter.ProcessMesh(mModelPart.ElementsArray());

        saveNodes(meshConverter);

        saveTriangles(meshConverter);
        retrieveResults();
        mInitialized = true;
    } else {
        mNodesCount = 0;
        mTrianglesCount = 0;
        mInitialized = false;
    }


}

void ModelPartWrapper::saveNodes(MeshConverter &meshConverter) {
    std::vector<int> nodes = meshConverter.GetNodes();
    mNodesCount = nodes.size();

    idTranslator.init(nodes);

    pmXCoordinates = new float[mNodesCount];
    pmYCoordinates = new float[mNodesCount];
    pmZCoordinates = new float[mNodesCount];
}

void ModelPartWrapper::saveTriangles(MeshConverter &meshConverter) {
    if (mModelPart.HasSubModelPart(SKIN_SUBMODEL_PART_NAME)) return;

    ModelPart &rSkinModelPart = mModelPart.CreateSubModelPart(SKIN_SUBMODEL_PART_NAME);

    int lastId = getMaxElementId();

    std::vector<face> faces = meshConverter.GetFaces();

    mTrianglesCount = faces.size();
    pmTriangles = new int[mTrianglesCount * 3];


    for (int i = 0; i < mTrianglesCount; i++) {
        std::vector<Kratos::IndexType> nodes;
        for (int j = 0; j < 3; j++) {
            pmTriangles[3 * i + j] = faces.at(i).nodes[j];
            nodes.push_back(idTranslator.getKratosId(faces.at(i).nodes[j]));
            rSkinModelPart.AddNode(mModelPart.pGetNode(idTranslator.getKratosId(faces.at(i).nodes[j])));
        }
        lastId++;
        rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", lastId, nodes,
                                          mModelPart.pGetProperties(0));
    }
    mMaxElementId = lastId;
    updateMaxElementId(mMaxElementId);
}

void ModelPartWrapper::retrieveResults() {
    auto &rSkinModelPart = mModelPart.GetSubModelPart(SKIN_SUBMODEL_PART_NAME);
    auto &rNodesArray = rSkinModelPart.Nodes();
    const auto nodeBegin = rNodesArray.begin();

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rNodesArray.size()); ++i) {
        auto it_node = nodeBegin + i;
        const Kratos::array_1d<double, 3> &r_coordinates = it_node->Coordinates();
        int current_node_id = idTranslator.getSurfaceId(it_node->Id());
        pmXCoordinates[current_node_id] = r_coordinates[0];
        pmYCoordinates[current_node_id] = r_coordinates[1];
        pmZCoordinates[current_node_id] = r_coordinates[2];
    }
    if (mStressResultsEnabled) {
        auto &rConditionsArray = rSkinModelPart.Conditions();
        const auto it_condition_begin = rConditionsArray.begin();
        for (int i = 0; i < static_cast<int>(rConditionsArray.size()); ++i) {
            std::vector<double> result;
            auto it_condition = it_condition_begin + i;

            it_condition->GetValue(Kratos::NEIGHBOUR_ELEMENTS)[0].CalculateOnIntegrationPoints(Kratos::VON_MISES_STRESS,
                                                                                               result,
                                                                                               mModelPart.GetProcessInfo());
            pmSurfaceStress[i] = result[0];
        }
    }

}

void ModelPartWrapper::enableSurfaceStressResults() {
    mStressResultsEnabled = true;
    pmSurfaceStress = new float[mTrianglesCount];

    Kratos::TetrahedralMeshOrientationCheck tetrahedralMeshOrientationCheck(mModelPart,
                                                                            false,
                                                                            Kratos::TetrahedralMeshOrientationCheck::ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS);
    tetrahedralMeshOrientationCheck.Execute();
}

float *ModelPartWrapper::getSurfaceStress() {
    return pmSurfaceStress;
}

ModelPartWrapper::~ModelPartWrapper() {
    delete[] pmXCoordinates;
    delete[] pmYCoordinates;
    delete[] pmZCoordinates;
}

ModelPartWrapper *ModelPartWrapper::getSubmodelPart(char *name) {
    return new ModelPartWrapper(mModelPart.GetSubModelPart(name), mFixedNodes, this);
}

ModelPart &ModelPartWrapper::getKratosModelPart() {
    return mModelPart;
}

int ModelPartWrapper::getMaxElementId() {
    if (pmParentModelPart == NULL) return mMaxElementId;
    else return pmParentModelPart->getMaxElementId();
}

int ModelPartWrapper::getMaxNodeId() {
    if (pmParentModelPart == NULL) return mMaxNodeId;
    else return pmParentModelPart->getMaxNodeId();
}

void ModelPartWrapper::updateMaxElementId(int maxId) {
    mMaxElementId = std::max(maxId, mMaxElementId);
    if (pmParentModelPart != NULL) pmParentModelPart->updateMaxElementId(mMaxElementId);
}

void ModelPartWrapper::updateMaxNodeId(int maxId) {
    mMaxNodeId = std::max(maxId, mMaxNodeId);
    if (pmParentModelPart != NULL) pmParentModelPart->updateMaxNodeId(mMaxNodeId);
}

void ModelPartWrapper::recreateProcessedMesh() {
    deleteSkin();
    if (mInitialized) {
        delete pmXCoordinates;
        delete pmYCoordinates;
        delete pmZCoordinates;
        delete pmTriangles;
    }
    if (mStressResultsEnabled) delete pmSurfaceStress;
    initialize();
}

ModelPartWrapper *ModelPartWrapper::createSubmodelPart(char *name) {
    Kratos::ModelPart &newModelPart = mModelPart.CreateSubModelPart(name);
    return new ModelPartWrapper(newModelPart, mFixedNodes);
}

NodeType *ModelPartWrapper::createNewNode(int id, double x, double y, double z) {
    updateMaxNodeId(id);
    return mModelPart.CreateNewNode(id, x, y, z).get();
}

ElementType *ModelPartWrapper::createNewElement(char *name, int id, int *nodeIds) {
    std::vector<Kratos::IndexType> nodes;
    for (int i = 0; i < 4; i++) nodes.push_back(nodeIds[i]);

    updateMaxElementId(id);
    return mModelPart.CreateNewElement(name, id, nodes, mModelPart.pGetProperties(0)).get();

}

ConditionType *ModelPartWrapper::createNew2dCondition(char *name, int id, int *nodeIds) {
    std::vector<Kratos::IndexType> nodes;
    for (int i = 0; i < 4; i++) nodes.push_back(nodeIds[i]);

    updateMaxElementId(id);
    return mModelPart.CreateNewCondition(name, id, nodes, mModelPart.pGetProperties(0)).get();

}

void ModelPartWrapper::removeNode(int id) {
    mModelPart.RemoveNode(id);
}

void ModelPartWrapper::removeElement(int id) {
    mModelPart.RemoveNode(id);
}

void ModelPartWrapper::removeCondition(int id) {
    mModelPart.RemoveNode(id);
}

void ModelPartWrapper::addNodes(int *nodeIds, int nodeCount) {
    std::vector<Kratos::IndexType> nodes;
    for (int i = 0; i < nodeCount; i++) nodes.push_back(nodeIds[i]);

    mModelPart.AddNodes(nodes);
}

void ModelPartWrapper::addElements(int *elementIds, int elementCount) {
    std::vector<Kratos::IndexType> elements;
    for (int i = 0; i < elementCount; i++) elements.push_back(elementIds[i]);

    mModelPart.AddElements(elements);
}

void ModelPartWrapper::addConditions(int *conditionIds, int conditionCount) {
    std::vector<Kratos::IndexType> conditions;
    for (int i = 0; i < conditionCount; i++)conditions.push_back(conditionIds[i]);

    mModelPart.AddConditions(conditions);
}

void ModelPartWrapper::deleteSkin() {
    if (mModelPart.HasSubModelPart(SKIN_SUBMODEL_PART_NAME)) {
        auto &rSkinModelPart = mModelPart.GetSubModelPart(SKIN_SUBMODEL_PART_NAME);
        auto &rConditionsArray = rSkinModelPart.Conditions();
        std::vector<Kratos::IndexType> toRemove;
        const auto it_condition_begin = rConditionsArray.begin();
        for (int i = 0; i < static_cast<int>(rConditionsArray.size()); ++i) {
            auto it_condition = it_condition_begin + i;
//            mModelPart.RemoveConditionFromAllLevels(it_condition->Id());
            toRemove.push_back(it_condition->Id());
        }
        for (unsigned int i : toRemove) mModelPart.RemoveConditionFromAllLevels(i);
        mModelPart.RemoveSubModelPart(SKIN_SUBMODEL_PART_NAME);
    }
}

NodeType *ModelPartWrapper::getNode(int id) {
    return &mModelPart.GetNode(id);
}

NodeType **ModelPartWrapper::getNodes() {
    int numberOfNodes = mModelPart.NumberOfNodes();
    auto **nodes = new NodeType *[numberOfNodes];
    auto nodeVector = mModelPart.Nodes().GetContainer();
    for (int i = 0; i < numberOfNodes; i++) {
        nodes[i] = nodeVector[i].get();
    }

    return nodes;
}

int ModelPartWrapper::getNumberOfNodes() {
    return mModelPart.NumberOfNodes();
}

ElementType *ModelPartWrapper::getElement(int id) {
    return &mModelPart.GetElement(id);
}

ElementType **ModelPartWrapper::getElements() {
    int numberOfElements = mModelPart.NumberOfElements();
    auto **elements = new ElementType *[numberOfElements];
    auto elementVector = mModelPart.Elements().GetContainer();

    for (int i = 0; i < numberOfElements; i++) elements[i] = elementVector[i].get();

    return elements;
}

int ModelPartWrapper::getNumberOfElements() {
    return mModelPart.NumberOfElements();
}

ConditionType *ModelPartWrapper::getCondition(int id) {
    return &mModelPart.GetCondition(id);
}

ConditionType **ModelPartWrapper::getConditions() {
    int numberOfConditions = mModelPart.NumberOfConditions();
    auto **conditions = new ConditionType *[numberOfConditions];
    auto conditionVector = mModelPart.Conditions().GetContainer();

    for (int i = 0; i < numberOfConditions; i++) conditions[i] = conditionVector[i].get();
    return conditions;
}

int ModelPartWrapper::getNumberOfConditions() {
    return mModelPart.NumberOfConditions();
}

IdTranslator *ModelPartWrapper::getIdTranslator() {
    return &idTranslator;
}

double *ModelPartWrapper::getNodalVariable1d(Kratos::Variable<double> &variable) {
    auto *values = new double[mNodesCount];

    auto &rSkinModelPart = mModelPart.GetSubModelPart(SKIN_SUBMODEL_PART_NAME);
    auto &rNodesArray = rSkinModelPart.Nodes();
    const auto nodeBegin = rNodesArray.begin();
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rNodesArray.size()); ++i) {
        auto it_node = nodeBegin + i;
        int current_node_id = idTranslator.getSurfaceId(it_node->Id());
        values[current_node_id] = it_node->FastGetSolutionStepValue(variable);
    }
    return values;
}

double *ModelPartWrapper::getNodalVariableComponent(
        Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > &variable) {
    auto *values = new double[mNodesCount];

    auto &rSkinModelPart = mModelPart.GetSubModelPart(SKIN_SUBMODEL_PART_NAME);
    auto &rNodesArray = rSkinModelPart.Nodes();
    const auto nodeBegin = rNodesArray.begin();
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rNodesArray.size()); ++i) {
        auto it_node = nodeBegin + i;
        int current_node_id = idTranslator.getSurfaceId(it_node->Id());
        values[current_node_id] = it_node->FastGetSolutionStepValue(variable);
    }
    return values;
}

double *ModelPartWrapper::getNodalVariable3d(Kratos::Variable<Kratos::array_1d<double, 3>> &variable) {
    auto *values = new double[mNodesCount];

    auto &rSkinModelPart = mModelPart.GetSubModelPart(SKIN_SUBMODEL_PART_NAME);
    auto &rNodesArray = rSkinModelPart.Nodes();
    const auto nodeBegin = rNodesArray.begin();
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rNodesArray.size()); ++i) {
        auto it_node = nodeBegin + i;
        int current_node_id = idTranslator.getSurfaceId(it_node->Id());
        auto value = it_node->FastGetSolutionStepValue(variable);
        for (int j = 0; j < 3; j++) values[current_node_id * 3 + j] = value[j];
    }
    return values;
}

bool ModelPartWrapper::hasNodalVariable1d(Kratos::Variable<double> &variable) {
    return mModelPart.HasNodalSolutionStepVariable(variable);
}


bool ModelPartWrapper::hasNodalVariable3d(Kratos::Variable<Kratos::array_1d<double, 3>> &variable) {
    return mModelPart.HasNodalSolutionStepVariable(variable);
}
