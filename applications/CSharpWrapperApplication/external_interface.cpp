#include "external_interface.h"

#include "includes/constitutive_law.h"

using namespace CSharpKratosWrapper;

int* CSharpInterface::pmTriangles;
int CSharpInterface::mTrianglesCount;
int CSharpInterface::mNodesCount;
float* CSharpInterface::pmXCoordinates;
float* CSharpInterface::pmYCoordinates;
float* CSharpInterface::pmZCoordinates;
KratosInternals CSharpInterface::mKratosInternals;
std::vector<NodeType::Pointer> CSharpInterface::mFixedNodes;
IdTranslator CSharpInterface::idTranslator;

//Fetch faces from mesh converter and create skin model part
void CSharpInterface::saveTriangles(MeshConverter& meshConverter) {

    ModelPart& r_skin_model_part = mKratosInternals.GetSkinModelPart();
    ModelPart& r_main_model_part = mKratosInternals.GetMainModelPart();
    int lastId = r_main_model_part.Elements().back().Id();

    std::vector<face> faces = meshConverter.GetFaces();

    mTrianglesCount = faces.size();
    pmTriangles = new int[mTrianglesCount * 3];

    for (int i = 0; i < mTrianglesCount; i++) {
        std::vector<Kratos::IndexType> nodes;
        for (int j = 0; j < 3; j++) {
            pmTriangles[3 * i + j] = faces.at(i).nodes[j];
            nodes.push_back(idTranslator.getKratosId(faces.at(i).nodes[j]));
            r_skin_model_part.AddNode(r_main_model_part.pGetNode(idTranslator.getKratosId(faces.at(i).nodes[j])));
        }
        lastId++;
        r_skin_model_part.CreateNewCondition("SurfaceCondition3D3N", lastId, nodes, r_main_model_part.pGetProperties(0));
    }
}

//Fetch surface nodes from mesh converter and initialize ID translator
void CSharpInterface::saveNodes(MeshConverter& meshConverter) {
    std::vector<int> nodes = meshConverter.GetNodes();
    mNodesCount = nodes.size();

    idTranslator.init(nodes);

    pmXCoordinates = new float[mNodesCount];
    pmYCoordinates = new float[mNodesCount];
    pmZCoordinates = new float[mNodesCount];
    retrieveNodesPos();
}

//Save recalculated surface nodes positions
void CSharpInterface::retrieveNodesPos() {
    ModelPart& r_skin_part = mKratosInternals.GetSkinModelPart();
    auto& r_nodes_array = r_skin_part.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    #pragma omp parallel for
    for (int i = 0; i<static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin+ i;
        const Kratos::array_1d<double, 3>& r_coordinates = it_node->Coordinates();
        int current_node_id = idTranslator.getUnityId(it_node->Id());
        pmXCoordinates[current_node_id] = r_coordinates[0];
        pmYCoordinates[current_node_id] = r_coordinates[1];
        pmZCoordinates[current_node_id] = r_coordinates[2];
    }
}

void CSharpInterface::freeNodes() {
    for (auto& node : mFixedNodes) {
        node->Free(Kratos::DISPLACEMENT_X);
        node->Free(Kratos::DISPLACEMENT_Y);
        node->Free(Kratos::DISPLACEMENT_Z);
    }
    mFixedNodes.clear();
}

void CSharpInterface::init(const char* MDPAFilePath, const char* JSONFilePath) {
    // Init internals
    mKratosInternals.initInternals();

    // Load settings
    const std::string json_file_name = JSONFilePath != NULL ? std::string(JSONFilePath) : "";
    mKratosInternals.loadSettingsParameters(json_file_name);

    // Init model part
    mKratosInternals.initModelPart();

    // Load MDPA file
    mKratosInternals.loadMDPA(std::string(MDPAFilePath));

    // Init dofs
    mKratosInternals.initDofs();

    // Read properties and materials
    mKratosInternals.initProperties();

    // Init solver
    mKratosInternals.initSolver();

    // Calling mesh converted
    MeshConverter meshConverter;
    meshConverter.ProcessMesh(mKratosInternals.GetMainModelPart().ElementsArray());

    // Save mesh
    saveNodes(meshConverter);
    saveTriangles(meshConverter);
}

void CSharpInterface::initWithSettings(const char* JSONFilePath) {
    // Init internals
    mKratosInternals.initInternals();

    // Load settings
    const std::string json_file_name = JSONFilePath != NULL ? std::string(JSONFilePath) : "";
    mKratosInternals.loadSettingsParameters(json_file_name);

    // Init model part
    mKratosInternals.initModelPart();

    // Load MDPA file
    const auto setting = mKratosInternals.GetSettings();
    mKratosInternals.loadMDPA(setting["solver_settings"]["model_import_settings"]["input_filename"].GetString());

    // Init dofs
    mKratosInternals.initDofs();

    // Read properties and materials
    mKratosInternals.initProperties();

    // Init solver
    mKratosInternals.initSolver();

    // Calling mesh converted
    MeshConverter meshConverter;
    meshConverter.ProcessMesh(mKratosInternals.GetMainModelPart().ElementsArray());

    // Save mesh
    saveNodes(meshConverter);
    saveTriangles(meshConverter);
}

//Update DISPLACEMENT variable of a node, so that final position is as given. X0 + DISPLACEMENT_X = x
void CSharpInterface::updateNodePos(const int nodeId, const float x, const float y, const float z) {

    // Get node
    NodeType::Pointer p_node =  mKratosInternals.GetMainModelPart().pGetNode(idTranslator.getKratosId(nodeId));

    // Fix nodes
    p_node->Fix(Kratos::DISPLACEMENT_X);
    p_node->Fix(Kratos::DISPLACEMENT_Y);
    p_node->Fix(Kratos::DISPLACEMENT_Z);

    // Arrays
    Kratos::array_1d<double, 3>& r_displacement = p_node->FastGetSolutionStepValue(Kratos::DISPLACEMENT);
    const Kratos::array_1d<double, 3>& r_initial_coordinates = p_node->GetInitialPosition().Coordinates();

    // Update coordinates
    Kratos::array_1d<double, 3>& r_coordinates = p_node->Coordinates();
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

void CSharpInterface::calculate() {
    mKratosInternals.solve();
    retrieveNodesPos();

    freeNodes();

}

float* CSharpInterface::getXCoordinates() {
    return pmXCoordinates;
}

float* CSharpInterface::getYCoordinates() {
    return pmYCoordinates;
}

float* CSharpInterface::getZCoordinates() {
    return pmZCoordinates;
}

int CSharpInterface::getNodesCount() {
    return mNodesCount;
}

int* CSharpInterface::getTriangles() {
    return pmTriangles;
}

int CSharpInterface::getTrianglesCount() {
    return mTrianglesCount;
}
