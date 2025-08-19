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
#include "tests/cpp_tests/convection_diffusion_fast_suite.h"
#include "containers/model.h"
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "convection_diffusion_application.h"
#include "../test_utilities/convection_diffusion_testing_utilities.h"

namespace Kratos
{
namespace Testing
{

    KRATOS_TEST_CASE_IN_SUITE(LaplacianElement2D3N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

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

        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-4)
        KRATOS_EXPECT_MATRIX_NEAR(LHS, expected_LHS, 1.0e-4)
    }

    KRATOS_TEST_CASE_IN_SUITE(LaplacianElement3D4N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

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

        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-4)
        KRATOS_EXPECT_MATRIX_NEAR(LHS, expected_LHS, 1.0e-4)
    }

} // namespace Testing
} // namespace Kratos.
