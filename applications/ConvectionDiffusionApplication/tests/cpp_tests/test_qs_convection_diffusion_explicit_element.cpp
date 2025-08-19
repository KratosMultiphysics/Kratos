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
#include "tests/cpp_tests/convection_diffusion_fast_suite.h"
#include "containers/model.h"
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "convection_diffusion_application.h"
#include "../test_utilities/convection_diffusion_testing_utilities.h"

namespace Kratos::Testing
{

    KRATOS_TEST_CASE_IN_SUITE(QSConvectionDiffusionExplicit2D3N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

        // Fill the process info container
        auto &r_process_info = r_test_model_part.GetProcessInfo();
        r_process_info.SetValue(DELTA_TIME, 0.1);
        r_process_info.SetValue(DYNAMIC_TAU, 1.0);
        r_process_info.SetValue(OSS_SWITCH, 1);
        r_process_info.SetValue(RUNGE_KUTTA_STEP, 4);

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
        p_element->Initialize(r_process_info);
        p_element->AddExplicitContribution(r_process_info);

        std::vector<double> expected_reaction({0.5, -0.5, 0.0});
        for (unsigned int i_node = 0; i_node < r_test_model_part.NumberOfNodes(); ++i_node) {
            const auto it_node = r_test_model_part.NodesBegin() + i_node;
            KRATOS_EXPECT_NEAR(it_node->FastGetSolutionStepValue(REACTION_FLUX), expected_reaction[i_node], 1.0e-6);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(QSConvectionDiffusionExplicitOSS2D3N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

        // Fill the process info container
        auto &r_process_info = r_test_model_part.GetProcessInfo();
        r_process_info.SetValue(DELTA_TIME, 0.1);
        r_process_info.SetValue(DYNAMIC_TAU, 1.0);
        r_process_info.SetValue(OSS_SWITCH, 1);
        r_process_info.SetValue(RUNGE_KUTTA_STEP, 4);

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
        double unknown_proj;
        p_element->Initialize(r_process_info);
        p_element->Calculate(PROJECTED_SCALAR1, unknown_proj, r_process_info);

        std::vector<double> expected_projection({-0.5, 0.5, 0.0});
        for (unsigned int i_node = 0; i_node < r_test_model_part.NumberOfNodes(); ++i_node) {
            const auto it_node = r_test_model_part.NodesBegin() + i_node;
            KRATOS_EXPECT_NEAR(it_node->GetValue(PROJECTED_SCALAR1), expected_projection[i_node], 1.0e-6);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(QSConvectionDiffusionExplicit3D4N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

        // Fill the process info container
        auto &r_process_info = r_test_model_part.GetProcessInfo();
        r_process_info.SetValue(DELTA_TIME, 0.1);
        r_process_info.SetValue(DYNAMIC_TAU, 1.0);
        r_process_info.SetValue(OSS_SWITCH, 1);
        r_process_info.SetValue(RUNGE_KUTTA_STEP, 4);

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
        p_element->Initialize(r_process_info);
        p_element->AddExplicitContribution(r_process_info);

        std::vector<double> expected_reaction({0.166667, -0.166667, 0.0, 0.0});
        for (unsigned int i_node = 0; i_node < r_test_model_part.NumberOfNodes(); ++i_node) {
            const auto it_node = r_test_model_part.NodesBegin() + i_node;
            KRATOS_EXPECT_NEAR(it_node->FastGetSolutionStepValue(REACTION_FLUX), expected_reaction[i_node], 1.0e-6);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(QSConvectionDiffusionExplicitOSS3D4N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

        // Fill the process info container
        auto &r_process_info = r_test_model_part.GetProcessInfo();
        r_process_info.SetValue(DELTA_TIME, 0.1);
        r_process_info.SetValue(DYNAMIC_TAU, 1.0);
        r_process_info.SetValue(OSS_SWITCH, 1);
        r_process_info.SetValue(RUNGE_KUTTA_STEP, 4);

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
        double unknown_proj;
        p_element->Initialize(r_process_info);
        p_element->Calculate(PROJECTED_SCALAR1, unknown_proj,r_process_info);
        std::vector<double> expected_projection({-0.166667, 0.166667, 0.0, 0.0});
        for (unsigned int i_node = 0; i_node < r_test_model_part.NumberOfNodes(); ++i_node) {
            const auto it_node = r_test_model_part.NodesBegin() + i_node;
            KRATOS_EXPECT_NEAR(it_node->GetValue(PROJECTED_SCALAR1), expected_projection[i_node], 1.0e-6);
        }
    }

} // namespace Kratos::Testing.

