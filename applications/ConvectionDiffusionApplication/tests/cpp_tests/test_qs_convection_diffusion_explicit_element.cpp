// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
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
    void SetQSConvectionDiffusionExplicitTestModelPart(ModelPart &rModelPart)
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
        r_process_info.SetValue(DYNAMIC_TAU, 1.0);
        r_process_info.SetValue(OSS_SWITCH, 1);
        r_process_info.SetValue(RUNGE_KUTTA_STEP, 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(QSConvectionDiffusionExplicit2D3N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetQSConvectionDiffusionExplicitTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
        r_test_model_part.CreateNewElement("QSConvectionDiffusionExplicit2D3N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(HEAT_FLUX) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
            array_1d<double,3> aux_vel = ZeroVector(3);
            aux_vel[0] = i_node.X();
            aux_vel[1] = i_node.Y();
            i_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
            i_node.FastGetSolutionStepValue(TEMPERATURE) = i_node.X();
            i_node.FastGetSolutionStepValue(TEMPERATURE,1) = i_node.Y();
        }

        // Test element
        auto p_element = r_test_model_part.pGetElement(1);
        const auto process_info = r_test_model_part.GetProcessInfo();
        p_element->Initialize(process_info);
        p_element->AddExplicitContribution(process_info);

        std::vector<double> expected_reaction({0.5, -0.5, 0.0});
        for (unsigned int i_node = 0; i_node < r_test_model_part.NumberOfNodes(); ++i_node) {
            const auto it_node = r_test_model_part.NodesBegin() + i_node;
            KRATOS_CHECK_NEAR(it_node->FastGetSolutionStepValue(REACTION_FLUX), expected_reaction[i_node], 1.0e-6);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(QSConvectionDiffusionExplicitOSS2D3N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetQSConvectionDiffusionExplicitTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
        r_test_model_part.CreateNewElement("QSConvectionDiffusionExplicit2D3N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(HEAT_FLUX) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
            array_1d<double,3> aux_vel = ZeroVector(3);
            aux_vel[0] = i_node.X();
            aux_vel[1] = i_node.Y();
            i_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
            i_node.FastGetSolutionStepValue(TEMPERATURE) = i_node.X();
            i_node.FastGetSolutionStepValue(TEMPERATURE,1) = i_node.Y();
        }

        // Test element
        auto p_element = r_test_model_part.pGetElement(1);
        const auto process_info = r_test_model_part.GetProcessInfo();
        double unknown_proj;
        p_element->Initialize(process_info);
        p_element->Calculate(PROJECTED_SCALAR1, unknown_proj,process_info);

        std::vector<double> expected_projection({-0.5, 0.5, 0.0});
        for (unsigned int i_node = 0; i_node < r_test_model_part.NumberOfNodes(); ++i_node) {
            const auto it_node = r_test_model_part.NodesBegin() + i_node;
            KRATOS_CHECK_NEAR(it_node->GetValue(PROJECTED_SCALAR1), expected_projection[i_node], 1.0e-6);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(QSConvectionDiffusionExplicit3D4N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetQSConvectionDiffusionExplicitTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        r_test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4};
        r_test_model_part.CreateNewElement("QSConvectionDiffusionExplicit3D4N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(HEAT_FLUX) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
            array_1d<double,3> aux_vel = ZeroVector(3);
            aux_vel[0] = i_node.X();
            aux_vel[1] = i_node.Y();
            aux_vel[2] = i_node.Z();
            i_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
            i_node.FastGetSolutionStepValue(TEMPERATURE) = i_node.X();
            i_node.FastGetSolutionStepValue(TEMPERATURE,1) = i_node.Y();
        }

        // Test element
        auto p_element = r_test_model_part.pGetElement(1);
        const auto process_info = r_test_model_part.GetProcessInfo();
        p_element->Initialize(process_info);
        p_element->AddExplicitContribution(process_info);

        std::vector<double> expected_reaction({0.166667, -0.166667, 0.0, 0.0});
        for (unsigned int i_node = 0; i_node < r_test_model_part.NumberOfNodes(); ++i_node) {
            const auto it_node = r_test_model_part.NodesBegin() + i_node;
            KRATOS_CHECK_NEAR(it_node->FastGetSolutionStepValue(REACTION_FLUX), expected_reaction[i_node], 1.0e-6);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(QSConvectionDiffusionExplicitOSS3D4N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetQSConvectionDiffusionExplicitTestModelPart(r_test_model_part);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        r_test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4};
        r_test_model_part.CreateNewElement("QSConvectionDiffusionExplicit3D4N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(HEAT_FLUX) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
            array_1d<double,3> aux_vel = ZeroVector(3);
            aux_vel[0] = i_node.X();
            aux_vel[1] = i_node.Y();
            aux_vel[2] = i_node.Z();
            i_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
            i_node.FastGetSolutionStepValue(TEMPERATURE) = i_node.X();
            i_node.FastGetSolutionStepValue(TEMPERATURE,1) = i_node.Y();
        }

        // Test element
        auto p_element = r_test_model_part.pGetElement(1);
        const auto process_info = r_test_model_part.GetProcessInfo();
        double unknown_proj;
        p_element->Initialize(process_info);
        p_element->Calculate(PROJECTED_SCALAR1, unknown_proj,process_info);
        std::vector<double> expected_projection({-0.166667, 0.166667, 0.0, 0.0});
        for (unsigned int i_node = 0; i_node < r_test_model_part.NumberOfNodes(); ++i_node) {
            const auto it_node = r_test_model_part.NodesBegin() + i_node;
            KRATOS_CHECK_NEAR(it_node->GetValue(PROJECTED_SCALAR1), expected_projection[i_node], 1.0e-6);
        }
    }


} // namespace Testing
} // namespace Kratos

