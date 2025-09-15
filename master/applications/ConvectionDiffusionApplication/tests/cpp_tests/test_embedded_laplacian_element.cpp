// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/variables.h"

// Application includes
#include "convection_diffusion_application.h"
#include "convection_diffusion_application_variables.h"
#include "../test_utilities/convection_diffusion_testing_utilities.h"

namespace Kratos::Testing
{

    KRATOS_TEST_CASE_IN_SUITE(EmbeddedLaplacianElement2D3N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test element
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

        // Add the DISTANCE variable (not added by SetEntityUnitTestModelPart)
        r_test_model_part.AddNodalSolutionStepVariable(DISTANCE);

        // Element creation
        r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
        r_test_model_part.CreateNewElement("EmbeddedLaplacianElement2D3N", 1, elem_nodes, r_test_model_part.pGetProperties(0));

        // Set the nodal values
        for (auto &i_node : r_test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(HEAT_FLUX) = 1.0;
            i_node.FastGetSolutionStepValue(CONDUCTIVITY) = 1.0;
        }

        // Get element and process info
        auto p_element = r_test_model_part.pGetElement(1);
        Vector RHS = ZeroVector(3);
        Matrix LHS = ZeroMatrix(3,3);
        const auto& r_process_info = r_test_model_part.GetProcessInfo();

        // Set Dirichlet penalty constant and boundary value
        p_element->SetValue(PENALTY_COEFFICIENT, 1e0);
        p_element->SetValue(EMBEDDED_SCALAR, 0.0);

        // Set distances for uncut element
        p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = 1.0;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;

        // Test uncut element
        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

        std::vector<double> expected_RHS = {0.166667, 0.166667, 0.166667};
        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-4);

        // Set distances for intersected element
        p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  0.5;

        // Test intersected element
        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

        expected_RHS = {0.00617284, 0.00617284, 0.0432099};
        KRATOS_EXPECT_VECTOR_NEAR(RHS, expected_RHS, 1.0e-4);
    }

} // namespace Kratos::Testing.
