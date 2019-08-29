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
    void SetEulerianConvDiffTestModelPart(ModelPart &rModelPart)
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
        r_process_info.SetValue(THETA, 1.0);
        r_process_info.SetValue(DELTA_TIME, 0.1);
        r_process_info.SetValue(DYNAMIC_TAU, 1.0);
    }

    KRATOS_TEST_CASE_IN_SUITE(EulerianConvDiff2D3NWithoutConvection, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetEulerianConvDiffTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
        r_test_model_part.CreateNewElement("EulerianConvDiff2D", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(DENSITY) = 1.0;
            i_node.FastGetSolutionStepValue(HEAT_FLUX) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
            i_node.FastGetSolutionStepValue(SPECIFIC_HEAT) = 1.0;
        }

        // Test element
        auto p_element = r_test_model_part.pGetElement(1);
        Vector RHS = ZeroVector(3);
        Matrix LHS = ZeroMatrix(3,3);
        p_element->CalculateLocalSystem(LHS, RHS, r_test_model_part.GetProcessInfo());

        std::vector<double> expected_RHS = {0.166667, 0.166667, 0.166667};
        std::vector<double> expected_LHS = {1.83333, -0.0833333, -0.0833333,
                                            -0.0833333, 1.33333, 0.416667,
                                            -0.0833333, 0.416667, 1.33333};

        for (unsigned int i = 0; i < 3; ++i) {
            KRATOS_CHECK_NEAR(RHS(i), expected_RHS[i], 1.0e-4);
            for (unsigned int j = 0; j < 3; ++j) {
                KRATOS_CHECK_NEAR(LHS(i,j), expected_LHS[i*3+j], 1.0e-4);
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(EulerianConvDiff2D3N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetEulerianConvDiffTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
        r_test_model_part.CreateNewElement("EulerianConvDiff2D", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(DENSITY) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
            i_node.FastGetSolutionStepValue(SPECIFIC_HEAT) = 1.0;
            array_1d<double,3> aux_vel = ZeroVector(3);
            aux_vel[0] = i_node.X();
            aux_vel[1] = i_node.Y();
            i_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
        }

        // Test element
        auto p_element = r_test_model_part.pGetElement(1);
        Vector RHS = ZeroVector(3);
        Matrix LHS = ZeroMatrix(3,3);
        p_element->CalculateLocalSystem(LHS, RHS, r_test_model_part.GetProcessInfo());

        std::vector<double> expected_RHS = {0.0, 0.0, 0.0};
        std::vector<double> expected_LHS = {1.72752,-0.0928335,-0.0928335,
                                            -0.197093,1.45074,0.475432,
                                            -0.197093,0.475432,1.45074};

        for (unsigned int i = 0; i < 3; ++i) {
            KRATOS_CHECK_NEAR(RHS(i), expected_RHS[i], 1.0e-4);
            for (unsigned int j = 0; j < 3; ++j) {
                KRATOS_CHECK_NEAR(LHS(i,j), expected_LHS[i*3+j], 1.0e-4);
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(EulerianConvDiff2D4N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetEulerianConvDiffTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
        r_test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4};
        r_test_model_part.CreateNewElement("EulerianConvDiff2D4N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(DENSITY) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
            i_node.FastGetSolutionStepValue(SPECIFIC_HEAT) = 1.0;
            array_1d<double,3> aux_vel = ZeroVector(3);
            aux_vel[0] = i_node.X();
            aux_vel[1] = i_node.Y();
            i_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
        }

        // Test element
        auto p_element = r_test_model_part.pGetElement(1);
        Vector RHS = ZeroVector(4);
        Matrix LHS = ZeroMatrix(4,4);
        p_element->CalculateLocalSystem(LHS, RHS, r_test_model_part.GetProcessInfo());

        std::vector<double> expected_RHS = {0.0, 0.0, 0.0, 0.0};
        std::vector<double> expected_LHS = {1.49976, 0.494317, -0.232728, 0.494317,
                                            0.430556, 1.67535, 0.680556, -0.286456,
                                            -0.360871, 0.616794, 1.87162, 0.616794,
                                            0.430556,-0.286456,0.680556,1.67535};

        for (unsigned int i = 0; i < 3; ++i) {
            KRATOS_CHECK_NEAR(RHS(i), expected_RHS[i], 1.0e-4);
            for (unsigned int j = 0; j < 3; ++j) {
                KRATOS_CHECK_NEAR(LHS(i,j), expected_LHS[i*4+j], 1.0e-4);
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(EulerianConvDiff3D4N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetEulerianConvDiffTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        r_test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4};
        r_test_model_part.CreateNewElement("EulerianConvDiff3D", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(DENSITY) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
            i_node.FastGetSolutionStepValue(SPECIFIC_HEAT) = 1.0;
            array_1d<double,3> aux_vel = ZeroVector(3);
            aux_vel[0] = i_node.X();
            aux_vel[1] = i_node.Y();
            i_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
        }

        // Test element
        auto p_element = r_test_model_part.pGetElement(1);
        Vector RHS = ZeroVector(4);
        Matrix LHS = ZeroMatrix(4,4);
        p_element->CalculateLocalSystem(LHS, RHS, r_test_model_part.GetProcessInfo());

        std::vector<double> expected_RHS = {0.0, 0.0, 0.0, 0.0};
        std::vector<double> expected_LHS = {0.646192,-0.0837032,-0.0837032,-0.0887235,
                                            -0.106429,0.355791,0.0945786,0.0860284,
                                            -0.106429,0.0945786,0.355791,0.0860284,
                                            -0.1,0.0916667,0.0916667,0.333333};

        for (unsigned int i = 0; i < 3; ++i) {
            KRATOS_CHECK_NEAR(RHS(i), expected_RHS[i], 1.0e-4);
            for (unsigned int j = 0; j < 3; ++j) {
                KRATOS_CHECK_NEAR(LHS(i,j), expected_LHS[i*4+j], 1.0e-4);
            }
        }
    }

} // namespace Testing
} // namespace Kratos.
