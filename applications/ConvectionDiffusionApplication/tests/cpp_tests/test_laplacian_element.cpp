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
    void SetLaplacianElementTestModelPart(ModelPart &rModelPart)
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
    }

    KRATOS_TEST_CASE_IN_SUITE(LaplacianElement2D3N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetLaplacianElementTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
        r_test_model_part.CreateNewElement("LaplacianElement2D3N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(HEAT_FLUX) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
        }

        // Test element
        auto p_element = r_test_model_part.pGetElement(1);
        Vector RHS = ZeroVector(3);
        Matrix LHS = ZeroMatrix(3,3);
        const auto& r_process_info = r_test_model_part.GetProcessInfo();

        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

        std::vector<double> expected_RHS = {0.166667, 0.166667, 0.166667};
        Matrix expected_LHS(3,3);
        expected_LHS(0, 0) = 1.0;
        expected_LHS(0, 1) = -0.5;
        expected_LHS(0, 2) = -0.5;
        expected_LHS(1, 0) = -0.5;
        expected_LHS(1, 1) = 0.5;
        expected_LHS(1, 2) = 0.0;
        expected_LHS(2, 0) = -0.5;
        expected_LHS(2, 1) = 0.0;
        expected_LHS(2, 2) = 0.5;

        KRATOS_CHECK_VECTOR_NEAR(RHS, expected_RHS, 1.0e-4)
        KRATOS_CHECK_MATRIX_NEAR(LHS, expected_LHS, 1.0e-4)
    }

    KRATOS_TEST_CASE_IN_SUITE(LaplacianElement3D4N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetLaplacianElementTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        r_test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4};
        r_test_model_part.CreateNewElement("LaplacianElement3D4N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(HEAT_FLUX) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
        }

        // Test element
        auto p_element = r_test_model_part.pGetElement(1);
        Vector RHS = ZeroVector(4);
        Matrix LHS = ZeroMatrix(4,4);
        const auto& r_process_info = r_test_model_part.GetProcessInfo();

        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

        std::vector<double> expected_RHS = {0.0416667, 0.0416667, 0.0416667, 0.0416667};
        Matrix expected_LHS = ZeroMatrix(4,4);
        expected_LHS(0, 0) = 0.5;
        expected_LHS(0, 1) = -0.166667;
        expected_LHS(0, 2) = -0.166667;
        expected_LHS(0, 3) = -0.166667;
        expected_LHS(1, 0) = -0.166667;
        expected_LHS(1, 1) = 0.166667;
        expected_LHS(2, 0) = -0.166667;
        expected_LHS(2, 2) = 0.166667;
        expected_LHS(3, 0) = -0.166667;
        expected_LHS(3, 3) = 0.166667;

        KRATOS_CHECK_VECTOR_NEAR(RHS, expected_RHS, 1.0e-4)
        KRATOS_CHECK_MATRIX_NEAR(LHS, expected_LHS, 1.0e-4)
    }

} // namespace Testing
} // namespace Kratos.
