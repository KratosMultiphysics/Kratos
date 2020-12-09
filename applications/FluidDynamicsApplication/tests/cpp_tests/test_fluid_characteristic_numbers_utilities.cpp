//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_characteristic_numbers_utilities.h"

namespace Kratos {
namespace Testing  {

namespace Internals {

void TestFluidCharacteristicNumberInitializeModelPart(
    ModelPart& rModelPart,
    const double DeltaTime)
{

    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    Properties::Pointer p_properties = rModelPart.CreateNewProperties(0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 2.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> element_nodes_1{1, 2, 3};
    std::vector<ModelPart::IndexType> element_nodes_2{3, 2, 4};
    rModelPart.CreateNewElement("Element2D3N", 1, element_nodes_1, p_properties);
    rModelPart.CreateNewElement("Element2D3N", 2, element_nodes_2, p_properties);

    // Set a fake current delta time
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, DeltaTime);

    // Set nodal data
    for (auto& rNode : rModelPart.Nodes()) {
        rNode.FastGetSolutionStepValue(VELOCITY_X) = rNode.Id() * rNode.X();
        rNode.FastGetSolutionStepValue(VELOCITY_Y) = rNode.Id() * rNode.Y();
    }
}

} // namespace internals

KRATOS_TEST_CASE_IN_SUITE(FluidCharacteristicNumbersUtilitiesCalculateLocalCFL, FluidDynamicsApplicationFastSuite)
{
    // Set the current delta time to calculate the CFL number
    const double current_dt = 1.0e-1;

    // Create the test model part
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    Internals::TestFluidCharacteristicNumberInitializeModelPart(r_model_part, current_dt);

    // Calculate the CFL number for each element
    FluidCharacteristicNumbersUtilities::CalculateLocalCFL(r_model_part);

    // Check results
    const double tolerance = 2.0e-6;
    KRATOS_CHECK_NEAR(r_model_part.GetElement(1).GetValue(CFL_NUMBER), 0.186339, tolerance);
    KRATOS_CHECK_NEAR(r_model_part.GetElement(2).GetValue(CFL_NUMBER), 0.792324, tolerance);
}

}
}