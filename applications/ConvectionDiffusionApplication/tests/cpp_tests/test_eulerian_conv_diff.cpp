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

    KRATOS_TEST_CASE_IN_SUITE(EulerianConvDiff2D3NNullConvection, KratosConvectionDiffusionFastSuite)
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

        const auto& rConstProcessInfoRef = r_test_model_part.GetProcessInfo();
        p_element->CalculateLocalSystem(LHS, RHS, rConstProcessInfoRef);

        std::vector<double> expected_RHS = {0.166667, 0.166667, 0.166667};
        Matrix expected_LHS(3,3);
        expected_LHS(0,0) = 1.83333;
        expected_LHS(0,1) = -0.0833333;
        expected_LHS(0,2) = -0.0833333;
        expected_LHS(1, 0) = -0.0833333;
        expected_LHS(1, 1) = 1.33333;
        expected_LHS(1, 2) = 0.416667;
        expected_LHS(2, 0) = -0.0833333;
        expected_LHS(2, 1) = 0.416667;
        expected_LHS(2, 2) = 1.33333;

        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-4)
        KRATOS_EXPECT_MATRIX_NEAR(LHS, expected_LHS, 1.0e-4)
    }

    KRATOS_TEST_CASE_IN_SUITE(EulerianConvDiff2D3N, KratosConvectionDiffusionFastSuite)
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
        const auto& rConstProcessInfoRef = r_test_model_part.GetProcessInfo();
        p_element->CalculateLocalSystem(LHS, RHS, rConstProcessInfoRef);

        std::vector<double> expected_RHS = {0.0, 0.0, 0.0};
        Matrix expected_LHS(3, 3);
        expected_LHS(0, 0) = 1.72752;
        expected_LHS(0, 1) = -0.0928335;
        expected_LHS(0, 2) = -0.0928335;
        expected_LHS(1, 0) = -0.197093;
        expected_LHS(1, 1) = 1.45074;
        expected_LHS(1, 2) = 0.475432;
        expected_LHS(2, 0) = -0.197093;
        expected_LHS(2, 1) = 0.475432;
        expected_LHS(2, 2) = 1.45074;

        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-4)
        KRATOS_EXPECT_MATRIX_NEAR(LHS, expected_LHS, 1.0e-4)

    }

    KRATOS_TEST_CASE_IN_SUITE(EulerianConvDiff2D4N, KratosConvectionDiffusionFastSuite)
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
        const auto& rConstProcessInfoRef = r_test_model_part.GetProcessInfo();
        p_element->CalculateLocalSystem(LHS, RHS, rConstProcessInfoRef);

        std::vector<double> expected_RHS = {0.0, 0.0, 0.0, 0.0};
        Matrix expected_LHS(4, 4);
        expected_LHS(0, 0) = 1.49976;
        expected_LHS(0, 1) = 0.494317;
        expected_LHS(0, 2) = -0.232728;
        expected_LHS(0, 3) = 0.494317;
        expected_LHS(1, 0) = 0.430556;
        expected_LHS(1, 1) = 1.67535;
        expected_LHS(1, 2) = 0.680556;
        expected_LHS(1, 3) = -0.286456;
        expected_LHS(2, 0) = -0.360871;
        expected_LHS(2, 1) = 0.616794;
        expected_LHS(2, 2) = 1.87162;
        expected_LHS(2, 3) = 0.616794;
        expected_LHS(3, 0) = 0.430556;
        expected_LHS(3, 1) = -0.286456;
        expected_LHS(3, 2) = 0.680556;
        expected_LHS(3, 3) = 1.67535;

        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-4)
        KRATOS_EXPECT_MATRIX_NEAR(LHS, expected_LHS, 1.0e-4)
    }

    KRATOS_TEST_CASE_IN_SUITE(EulerianConvDiff3D4N, KratosConvectionDiffusionFastSuite)
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
        const auto& rConstProcessInfoRef = r_test_model_part.GetProcessInfo();
        p_element->CalculateLocalSystem(LHS, RHS, rConstProcessInfoRef);

        std::vector<double> expected_RHS = {0.0, 0.0, 0.0, 0.0};
        Matrix expected_LHS(4, 4);
        expected_LHS(0, 0) = 0.646192;
        expected_LHS(0, 1) = -0.0837032;
        expected_LHS(0, 2) = -0.0837032;
        expected_LHS(0, 3) = -0.0887235;
        expected_LHS(1, 0) = -0.106429;
        expected_LHS(1, 1) = 0.355791;
        expected_LHS(1, 2) = 0.0945786;
        expected_LHS(1, 3) = 0.0860284;
        expected_LHS(2, 0) = -0.106429;
        expected_LHS(2, 1) = 0.0945786;
        expected_LHS(2, 2) = 0.355791;
        expected_LHS(2, 3) = 0.0860284;
        expected_LHS(3, 0) = -0.1;
        expected_LHS(3, 1) = 0.0916667;
        expected_LHS(3, 2) = 0.0916667;
        expected_LHS(3, 3) = 0.333333;

        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-4)
        KRATOS_EXPECT_MATRIX_NEAR(LHS, expected_LHS, 1.0e-4)
    }

} // namespace Kratos::Testing.
