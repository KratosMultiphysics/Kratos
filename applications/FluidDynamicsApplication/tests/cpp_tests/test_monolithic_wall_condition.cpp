//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Raul Bravo
//

// System includes
#include <sstream>

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "input_output/logger.h"

// Application includes
#include "fluid_dynamics_application.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(MonolithicWallCondition2D2NOutletInflowWarning, FluidDynamicsApplicationFastSuite)
{
    // Create the test model part
    Model model;
    auto& r_model_part = model.CreateModelPart("TestModelPart", 2);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    auto p_properties = r_model_part.CreateNewProperties(0);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
        r_node.AddDof(VELOCITY_Z);
        r_node.AddDof(PRESSURE);
    }

    // Create an outlet condition with the outlet inflow contribution switched on
    auto p_condition = r_model_part.CreateNewCondition("MonolithicWallCondition2D2N", 1, {{1, 2}}, p_properties);
    p_condition->Set(OUTLET, true);
    r_model_part.GetProcessInfo().SetValue(OUTLET_INFLOW_CONTRIBUTION_SWITCH, true);

    // Check the condition and capture the log
    std::stringstream buffer;
    auto p_output = Kratos::make_shared<LoggerOutput>(buffer);
    Logger::AddOutput(p_output);
    p_condition->Check(r_model_part.GetProcessInfo());
    Logger::RemoveOutput(p_output);

    // Check that the warning is printed
    if (Testing::GetDefaultDataCommunicator().Rank() == 0) {
        KRATOS_EXPECT_NE(buffer.str().find("outlet inflow contribution"), std::string::npos);
    }
}

} // namespace Kratos::Testing
