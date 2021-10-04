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
#include <iomanip> // for std::setprecision

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "convection_diffusion_application.h"


namespace Kratos
{
namespace Testing
{

void SetThermalCouplingConditionTestModelPart(ModelPart &rModelPart)
{
    // Set buffer size
    rModelPart.SetBufferSize(2);

    // Set convection diffusion settings
    ConvectionDiffusionSettings::Pointer p_conv_dff_set = Kratos::make_shared<ConvectionDiffusionSettings>();
    p_conv_dff_set->SetDensityVariable(DENSITY);
    p_conv_dff_set->SetDiffusionVariable(CONDUCTIVITY);
    p_conv_dff_set->SetUnknownVariable(TEMPERATURE);
    p_conv_dff_set->SetVolumeSourceVariable(HEAT_FLUX);
    p_conv_dff_set->SetSurfaceSourceVariable(FACE_HEAT_FLUX);
    p_conv_dff_set->SetProjectionVariable(PROJECTED_SCALAR1);
    p_conv_dff_set->SetConvectionVariable(CONVECTION_VELOCITY);
    p_conv_dff_set->SetMeshVelocityVariable(MESH_VELOCITY);
    p_conv_dff_set->SetVelocityVariable(VELOCITY);
    p_conv_dff_set->SetSpecificHeatVariable(SPECIFIC_HEAT);
    p_conv_dff_set->SetReactionVariable(REACTION_FLUX);
    rModelPart.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_dff_set);

    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);
    rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
    rModelPart.AddNodalSolutionStepVariable(HEAT_FLUX);
    rModelPart.AddNodalSolutionStepVariable(FACE_HEAT_FLUX);
    rModelPart.AddNodalSolutionStepVariable(PROJECTED_SCALAR1);
    rModelPart.AddNodalSolutionStepVariable(CONVECTION_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    rModelPart.AddNodalSolutionStepVariable(REACTION_FLUX);
}

KRATOS_TEST_CASE_IN_SUITE(ThermalCouplingCondition2D4N, KratosConvectionDiffusionFastSuite)
{
    // Create the test condition
    Model model;
    ModelPart &test_model_part = model.CreateModelPart("TestModelPart");
    SetThermalCouplingConditionTestModelPart(test_model_part);

    // Set the condition properties
    auto p_cond_prop = test_model_part.CreateNewProperties(0);
    p_cond_prop->SetValue(TRANSFER_COEFFICIENT, 0.5);

    // Condition creation
    test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    test_model_part.CreateNewNode(3, 1.0, 0.0, 0.0);
    test_model_part.CreateNewNode(4, 0.0, 0.0, 0.0);
    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4};
    test_model_part.CreateNewCondition("ThermalCouplingCondition2D4N", 1, cond_nodes, p_cond_prop);

    // Set the face heat flux
    double aux_temp = 1.0;
    for (auto &r_node : test_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = aux_temp;
        aux_temp *= 2.0;
    }

    // Test condition
    Condition::Pointer p_condition = test_model_part.pGetCondition(1);
    Vector RHS = ZeroVector(4);
    Matrix LHS = ZeroMatrix(4,4);
    const auto& r_process_info = test_model_part.GetProcessInfo();
    p_condition->CalculateLocalSystem(LHS, RHS, r_process_info);

    std::vector<double> expected_RHS = {0.75,1.5,-0.75,-1.5};
    std::vector<double> expected_LHS_row_0 = {0.25,0,-0.25,0};
    KRATOS_CHECK_VECTOR_NEAR(RHS, expected_RHS, 1.0e-8)
    KRATOS_CHECK_VECTOR_NEAR(row(LHS,0), expected_LHS_row_0, 1.0e-8)
}

KRATOS_TEST_CASE_IN_SUITE(ThermalCouplingCondition3D6N, KratosConvectionDiffusionFastSuite)
{
    // Create the test condition
    Model model;
    ModelPart &test_model_part = model.CreateModelPart("TestModelPart");
    SetThermalCouplingConditionTestModelPart(test_model_part);

    // Set the condition properties
    auto p_cond_prop = test_model_part.CreateNewProperties(0);
    p_cond_prop->SetValue(TRANSFER_COEFFICIENT, 0.5);

    // Condition creation
    test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    test_model_part.CreateNewNode(4, 0.0, 0.0, 0.0);
    test_model_part.CreateNewNode(5, 1.0, 0.0, 0.0);
    test_model_part.CreateNewNode(6, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5, 6};
    test_model_part.CreateNewCondition("ThermalCouplingCondition3D6N", 1, cond_nodes, p_cond_prop);

    // Set the face heat flux
    double aux_temp = 1.0;
    for (auto &i_node : test_model_part.Nodes()) {
        i_node.FastGetSolutionStepValue(TEMPERATURE) = aux_temp;
        aux_temp *= 2.0;
    }

    // Test condition
    Condition::Pointer p_condition = test_model_part.pGetCondition(1);
    Vector RHS = ZeroVector(6);
    Matrix LHS = ZeroMatrix(6,6);
    const auto& r_process_info = test_model_part.GetProcessInfo();
    p_condition->CalculateLocalSystem(LHS, RHS, r_process_info);

    std::vector<double> expected_RHS = {0.5833333333,1.166666667,2.333333333,-0.5833333333,-1.166666667,-2.333333333};
    std::vector<double> expected_LHS_row_0 = {0.08333333333,0,0,-0.08333333333,0,0};
    KRATOS_CHECK_VECTOR_NEAR(RHS, expected_RHS, 1.0e-8)
    KRATOS_CHECK_VECTOR_NEAR(row(LHS,0), expected_LHS_row_0, 1.0e-8)
}

} // namespace Testing
} // namespace Kratos.
