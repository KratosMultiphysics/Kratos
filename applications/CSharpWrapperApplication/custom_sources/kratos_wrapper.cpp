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

// External includes

// Project includes
#include "custom_includes/kratos_wrapper.h"
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
