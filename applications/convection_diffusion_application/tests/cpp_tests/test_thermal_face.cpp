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


namespace Kratos
{
namespace Testing
{
    KRATOS_TEST_CASE_IN_SUITE(ThermalFace2D, KratosConvectionDiffusionFastSuite)
    {
        // Create the test condition
        Model model;
        ModelPart &test_model_part = model.CreateModelPart("TestModelPart");
        test_model_part.SetBufferSize(2);

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
        test_model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_dff_set);

        // Variables addition
        test_model_part.AddNodalSolutionStepVariable(DENSITY);
        test_model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);
        test_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
        test_model_part.AddNodalSolutionStepVariable(HEAT_FLUX);
        test_model_part.AddNodalSolutionStepVariable(FACE_HEAT_FLUX);
        test_model_part.AddNodalSolutionStepVariable(PROJECTED_SCALAR1);
        test_model_part.AddNodalSolutionStepVariable(CONVECTION_VELOCITY);
        test_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
        test_model_part.AddNodalSolutionStepVariable(VELOCITY);
        test_model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
        test_model_part.AddNodalSolutionStepVariable(REACTION_FLUX);

        // Set the condition properties
        Properties::Pointer p_cond_prop = test_model_part.pGetProperties(0);
        p_cond_prop->SetValue(EMISSIVITY, 1.0);
        p_cond_prop->SetValue(AMBIENT_TEMPERATURE, 293.0);
        p_cond_prop->SetValue(CONVECTION_COEFFICIENT, 20.0);

        // Condition creation
        test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        std::vector<ModelPart::IndexType> cond_nodes{1, 2};
        test_model_part.CreateNewCondition("ThermalFace2D", 1, cond_nodes, p_cond_prop);

        // Set the face heat flux
        for (auto &i_node : test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(TEMPERATURE) = 400.0;
            i_node.FastGetSolutionStepValue(FACE_HEAT_FLUX) = 200.0;
            // i_node.FastGetSolutionStepValue(TEMPERATURE) = i_node.X() * 400.0;
            // i_node.FastGetSolutionStepValue(FACE_HEAT_FLUX) = 200.0 + i_node.X() * 100.0;
        }

        // Test condition
        Condition::Pointer p_condition = test_model_part.pGetCondition(1);
        Vector RHS = ZeroVector(2);
        Matrix LHS = ZeroMatrix(2,2);
        p_condition->CalculateLocalSystem(LHS, RHS, test_model_part.GetProcessInfo());

        std::vector<double> expected_RHS = {-1486.82, -1486.82};
        std::vector<double> expected_LHS = {17.2576, 0.0, 0.0, 17.2576};
        for (unsigned int i = 0; i < 2; ++i) {
            KRATOS_CHECK_NEAR(RHS(i), expected_RHS[i], 1.0e-2);
            for (unsigned int j = 0; j < 2; ++j) {
                KRATOS_CHECK_NEAR(LHS(i,j), expected_LHS[i*2+j], 1.0e-4);
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(ThermalFace2D2N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test condition
        Model model;
        ModelPart &test_model_part = model.CreateModelPart("TestModelPart");
        test_model_part.SetBufferSize(2);

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
        test_model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_dff_set);

        // Variables addition
        test_model_part.AddNodalSolutionStepVariable(DENSITY);
        test_model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);
        test_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
        test_model_part.AddNodalSolutionStepVariable(HEAT_FLUX);
        test_model_part.AddNodalSolutionStepVariable(FACE_HEAT_FLUX);
        test_model_part.AddNodalSolutionStepVariable(PROJECTED_SCALAR1);
        test_model_part.AddNodalSolutionStepVariable(CONVECTION_VELOCITY);
        test_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
        test_model_part.AddNodalSolutionStepVariable(VELOCITY);
        test_model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
        test_model_part.AddNodalSolutionStepVariable(REACTION_FLUX);

        // Set the condition properties
        Properties::Pointer p_cond_prop = test_model_part.pGetProperties(0);
        p_cond_prop->SetValue(EMISSIVITY, 1.0);
        p_cond_prop->SetValue(AMBIENT_TEMPERATURE, 293.0);
        p_cond_prop->SetValue(CONVECTION_COEFFICIENT, 20.0);

        // Condition creation
        test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        std::vector<ModelPart::IndexType> cond_nodes{1, 2};
        test_model_part.CreateNewCondition("ThermalFace2D2N", 1, cond_nodes, p_cond_prop);

        // Set the face heat flux
        for (auto &i_node : test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(TEMPERATURE) = 400.0;
            i_node.FastGetSolutionStepValue(FACE_HEAT_FLUX) = 200.0;
            // i_node.FastGetSolutionStepValue(TEMPERATURE) = i_node.X() * 400.0;
            // i_node.FastGetSolutionStepValue(FACE_HEAT_FLUX) = 200.0 + i_node.X() * 100.0;
        }

        // Test condition
        Condition::Pointer p_condition = test_model_part.pGetCondition(1);
        Vector RHS = ZeroVector(2);
        Matrix LHS = ZeroMatrix(2,2);
        p_condition->CalculateLocalSystem(LHS, RHS, test_model_part.GetProcessInfo());

        std::vector<double> expected_RHS = {-1486.82, -1486.82};
        std::vector<double> expected_LHS = {11.5051, 5.75253, 5.75253, 11.5051};
        for (unsigned int i = 0; i < 2; ++i) {
            KRATOS_CHECK_NEAR(RHS(i), expected_RHS[i], 1.0e-2);
            for (unsigned int j = 0; j < 2; ++j) {
                KRATOS_CHECK_NEAR(LHS(i,j), expected_LHS[i*2+j], 1.0e-4);
            }
        }
    }
} // namespace Testing
} // namespace Kratos.
