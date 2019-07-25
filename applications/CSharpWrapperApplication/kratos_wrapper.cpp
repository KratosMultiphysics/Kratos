
#include "kratos_wrapper.h"
#include "processes/tetrahedral_mesh_orientation_check.h"

using namespace CSharpKratosWrapper;

void KratosWrapper::init(const char *MDPAFilePath, const char *JSONFilePath) {
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

    // Create modelPart wrapper
    pmMainModelPartWrapper = new ModelPartWrapper(mKratosInternals.GetMainModelPart(), mFixedNodes);

}

void KratosWrapper::initWithSettings(const char *JSONFilePath) {
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

    // Create modelPart wrapper
    pmMainModelPartWrapper = new ModelPartWrapper(mKratosInternals.GetMainModelPart(), mFixedNodes);
}

void KratosWrapper::freeNodes() {
    for (auto &node : mFixedNodes) {
        node->Free(Kratos::DISPLACEMENT_X);
        node->Free(Kratos::DISPLACEMENT_Y);
        node->Free(Kratos::DISPLACEMENT_Z);
    }
    mFixedNodes.clear();
}

void KratosWrapper::calculate() {
    mKratosInternals.solve();
    freeNodes();
}

KratosWrapper::~KratosWrapper() {
    delete pmMainModelPartWrapper;
}

ModelPartWrapper *KratosWrapper::getRootModelPartWrapper() {
    return pmMainModelPartWrapper;
}
