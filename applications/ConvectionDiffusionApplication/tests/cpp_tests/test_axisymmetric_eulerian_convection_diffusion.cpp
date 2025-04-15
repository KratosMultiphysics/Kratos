// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "convection_diffusion_application.h"
#include "../test_utilities/convection_diffusion_testing_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(AxisymmetricEulerianConvectionDiffusion2D3N, KratosConvectionDiffusionFastSuite)
{
    // Create the test element
    Model model;
    auto &r_test_model_part = model.CreateModelPart("TestModelPart");
    ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

    // Fill the process info container
    auto &r_process_info = r_test_model_part.GetProcessInfo();
    r_process_info.SetValue(TIME_INTEGRATION_THETA, 1.0);
    r_process_info.SetValue(DELTA_TIME, 0.1);
    r_process_info.SetValue(DYNAMIC_TAU, 1.0);

    // Element creation
    r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
    r_test_model_part.CreateNewElement("AxisymmetricEulerianConvectionDiffusion2D3N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

    // Set the nodal values
    for (auto &r_node : r_test_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DENSITY) = 1.0;
        r_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
        r_node.FastGetSolutionStepValue(SPECIFIC_HEAT) = 1.0;
        array_1d<double,3> aux_vel = ZeroVector(3);
        aux_vel[0] = r_node.X();
        aux_vel[1] = r_node.Y();
        r_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
        r_node.FastGetSolutionStepValue(TEMPERATURE) = r_node.Id();
    }

    // Calculate RHS and LHS
    Vector RHS;
    Matrix LHS;
    auto p_element = r_test_model_part.pGetElement(1);
    p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

    // Check values
    Matrix expected_LHS(3, 3);
    std::vector<double> expected_RHS = {-2.29951736772, -7.10404985495, -15.7198429261};
    expected_LHS(0, 0) = 2.95807430822; expected_LHS(0, 1) = -0.413676492019; expected_LHS(0, 2) = 0.0562653478447;
    expected_LHS(1, 0) = -0.639773806745; expected_LHS(1, 1) = 2.18994627326; expected_LHS(1, 2) = 1.12131037172;
    expected_LHS(2, 0) = -0.417168056159; expected_LHS(2, 1) = 1.14274872206; expected_LHS(2, 2) = 4.61717117937;
    KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-8)
    KRATOS_EXPECT_MATRIX_NEAR(LHS, expected_LHS, 1.0e-8)
}

KRATOS_TEST_CASE_IN_SUITE(AxisymmetricEulerianConvectionDiffusion2D4N, KratosConvectionDiffusionFastSuite)
{
    // Create the test element
    Model model;
    auto &r_test_model_part = model.CreateModelPart("TestModelPart");
    ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

    // Fill the process info container
    auto &r_process_info = r_test_model_part.GetProcessInfo();
    r_process_info.SetValue(TIME_INTEGRATION_THETA, 1.0);
    r_process_info.SetValue(DELTA_TIME, 0.1);
    r_process_info.SetValue(DYNAMIC_TAU, 1.0);

    // Element creation
    r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4};
    r_test_model_part.CreateNewElement("AxisymmetricEulerianConvectionDiffusion2D4N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

    // Set the nodal values
    for (auto &r_node : r_test_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DENSITY) = 1.0;
        r_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
        r_node.FastGetSolutionStepValue(SPECIFIC_HEAT) = 1.0;
        array_1d<double,3> aux_vel = ZeroVector(3);
        aux_vel[0] = r_node.X();
        aux_vel[1] = r_node.Y();
        r_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
        r_node.FastGetSolutionStepValue(TEMPERATURE) = r_node.Id();
    }

    // Calculate RHS and LHS
    Vector RHS;
    Matrix LHS;
    auto p_element = r_test_model_part.pGetElement(1);
    p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

    // Check values
    Matrix expected_LHS(4, 4);
    std::vector<double> expected_RHS = {-9.79117129173, -10.9554458519, -34.4046862197, -39.3898359347};
    expected_LHS(0, 0) = 3.07670654707; expected_LHS(0, 1) = 0.85904076018; expected_LHS(0, 2) = -0.0302338911763; expected_LHS(0, 3) = 1.27177122446;
    expected_LHS(1, 0) = 0.640532577755; expected_LHS(1, 1) = 3.38630181886; expected_LHS(1, 2) = 1.54644080045; expected_LHS(1, 3) = -0.27425319123;
    expected_LHS(2, 0) = -0.565733544955; expected_LHS(2, 1) = 0.985473628635; expected_LHS(2, 2) = 9.14090890724; expected_LHS(2, 3) = 1.39418644642;
    expected_LHS(3, 0) = 0.677486900786; expected_LHS(3, 1) = -0.334321094873; expected_LHS(3, 2) = 2.12547106059; expected_LHS(3, 3) = 8.25114451047;
    KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-8)
    KRATOS_EXPECT_MATRIX_NEAR(LHS, expected_LHS, 1.0e-8)
}

} // namespace Kratos::Testing.
