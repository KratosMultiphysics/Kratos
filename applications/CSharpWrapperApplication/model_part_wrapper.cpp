
#include <structural_mechanics_application_variables.h>
#include <stdlib.h>
#include "model_part_wrapper.h"

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
    return pmTriangles;
    auto &r_conditions_array = mModelPart.GetSubModelPart(SKIN_SUBMODEL_PART_NAME).Conditions();
    const auto it_condition_begin = r_conditions_array.begin();
    int *triangles = new int[mTrianglesCount];
    for (int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_condition = it_condition_begin + i;
        auto geometry = it_condition->GetGeometry();

        for (int j = 0; j < 3; j++) {
            triangles[3 * i + j] = geometry.GetPoint(j).Id();
        }
    }
    return triangles;
}

int ModelPartWrapper::getTrianglesCount() {
    return mTrianglesCount;
}

//Update DISPLACEMENT variable of a node, so that final position is as given. X0 + DISPLACEMENT_X = x
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
    mMaxId = mModelPart.Elements().back().Id();
    MeshConverter meshConverter;
    meshConverter.ProcessMesh(mModelPart.ElementsArray());

    saveNodes(meshConverter);

    saveTriangles(meshConverter);
    mStressResultsEnabled = false;

    retrieveResults();
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

    int lastId = getMaxId();

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
    mMaxId = lastId;
    updateMaxId(mMaxId);
}

void ModelPartWrapper::retrieveResults() {
    auto &rSkinModelPart = mModelPart.GetSubModelPart(SKIN_SUBMODEL_PART_NAME);
    auto &rNodesArray = rSkinModelPart.Nodes();
    const auto nodeBegin = rNodesArray.begin();

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rNodesArray.size()); ++i) {
        auto it_node = nodeBegin + i;
        const Kratos::array_1d<double, 3> &r_coordinates = it_node->Coordinates();
        int current_node_id = idTranslator.getUnityId(it_node->Id());
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

int ModelPartWrapper::getMaxId() {
    if (pmParentModelPart == NULL) return mMaxId;
    else return pmParentModelPart->getMaxId();
}

void ModelPartWrapper::updateMaxId(int maxId) {
    mMaxId = std::max(maxId, mMaxId);
    if (pmParentModelPart != NULL) pmParentModelPart->updateMaxId(mMaxId);
}
