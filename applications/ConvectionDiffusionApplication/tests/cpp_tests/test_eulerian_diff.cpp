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


namespace Kratos
{
namespace Testing
{
    void SetEulerianDiffTestModelPart(ModelPart &rModelPart)
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

        // Create a fake properties container
        auto p_elem_prop = rModelPart.CreateNewProperties(0);

        // Fill the process info container
        auto &r_process_info = rModelPart.GetProcessInfo();
        r_process_info.SetValue(DELTA_TIME, 0.1);
    }

    KRATOS_TEST_CASE_IN_SUITE(EulerianDiff2D3N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetEulerianDiffTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
        r_test_model_part.CreateNewElement("EulerianDiffusion2D3N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(DENSITY) = 1.0;
            i_node.FastGetSolutionStepValue(HEAT_FLUX) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
            i_node.FastGetSolutionStepValue(SPECIFIC_HEAT) = 1.0;
            i_node.FastGetSolutionStepValue(PROJECTED_SCALAR1) = 1.0;
            i_node.FastGetSolutionStepValue(PROJECTED_SCALAR1, 1) = 0.5;
        }

        // Set expected solution
        std::vector<double> expected_RHS = {5.0/3.0, 5.0/3.0, 5.0/3.0};
        Matrix expected_LHS(3, 3);
        expected_LHS(0, 0) = 1.333333333333;
        expected_LHS(0, 1) = 5.0/30.0;
        expected_LHS(0, 2) = 5.0/30.0;
        expected_LHS(1, 0) = 5.0/30.0;
        expected_LHS(1, 1) = 1.083333333333;
        expected_LHS(1, 2) = 0.416666666666;
        expected_LHS(2, 0) = 5.0/30.0;
        expected_LHS(2, 1) = 0.416666666666;
        expected_LHS(2, 2) = 1.083333333333;

        // Test CalculateLocalSystem
        auto p_element = r_test_model_part.pGetElement(1);
        Vector RHS = ZeroVector(3);
        Matrix LHS = ZeroMatrix(3,3);

        const auto& r_process_info = r_test_model_part.GetProcessInfo();
        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);
        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-12)
        KRATOS_EXPECT_MATRIX_NEAR(LHS, expected_LHS, 1.0e-12)

        // Test CalculateRightHandSide
        p_element->CalculateRightHandSide(RHS, r_process_info);
        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-12)
    }

    KRATOS_TEST_CASE_IN_SUITE(EulerianDiff3D4N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetEulerianDiffTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        r_test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4};
        r_test_model_part.CreateNewElement("EulerianDiffusion3D4N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(DENSITY) = 1.0;
            i_node.FastGetSolutionStepValue(HEAT_FLUX) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
            i_node.FastGetSolutionStepValue(SPECIFIC_HEAT) = 1.0;
            i_node.FastGetSolutionStepValue(PROJECTED_SCALAR1) = 1.0;
            i_node.FastGetSolutionStepValue(PROJECTED_SCALAR1, 1) = 0.5;
        }

        // Set expected solution
        std::vector<double> expected_RHS = {0.416666666666, 0.416666666666, 0.416666666666, 0.416666666666};
        Matrix expected_LHS(4, 4);
        expected_LHS(0, 0) = 0.416666666666;
        expected_LHS(0, 1) = -4.19267e-10;
        expected_LHS(0, 2) = -4.19267e-10;
        expected_LHS(0, 3) = -4.19267e-10;
        expected_LHS(1, 0) = -4.19267e-10;
        expected_LHS(1, 1) = 0.25;
        expected_LHS(1, 2) = 0.083333333333;
        expected_LHS(1, 3) = 0.083333333333;
        expected_LHS(2, 0) = -4.19267e-10;
        expected_LHS(2, 1) = 0.083333333333;
        expected_LHS(2, 2) = 0.25;
        expected_LHS(2, 3) = 0.083333333333;
        expected_LHS(3, 0) = -4.19267e-10;
        expected_LHS(3, 1) = 0.083333333333;
        expected_LHS(3, 2) = 0.083333333333;
        expected_LHS(3, 3) = 0.25;

        // Test CalculateLocalSystem
        auto p_element = r_test_model_part.pGetElement(1);
        Vector RHS = ZeroVector(4);
        Matrix LHS = ZeroMatrix(4,4);
        const auto& r_process_info = r_test_model_part.GetProcessInfo();
        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);
        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-12)
        KRATOS_EXPECT_MATRIX_NEAR(LHS, expected_LHS, 1.0e-8)

        // Test CalculateRightHandSide
        p_element->CalculateRightHandSide(RHS, r_process_info);
        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-12)
    }

} // namespace Testing
} // namespace Kratos.
