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

    KRATOS_TEST_CASE_IN_SUITE(ThermalFace2D2N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test condition
        Model model;
        ModelPart &test_model_part = model.CreateModelPart("TestModelPart");
        ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(test_model_part);

        // Set the condition properties
        auto p_cond_prop = test_model_part.pGetProperties(0);
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
        Vector RHS = ZeroVector(3);
        Matrix LHS = ZeroMatrix(3,3);
        const auto& rConstProcessRef = test_model_part.GetProcessInfo();
        p_condition->CalculateLocalSystem(LHS, RHS, rConstProcessRef);

        std::vector<double> expected_RHS = {-1486.82, -1486.82};
        std::vector<double> expected_LHS = {11.5051, 5.75253,
                                            5.75253, 11.5051};
        for (unsigned int i = 0; i < 2; ++i) {
            KRATOS_EXPECT_NEAR(RHS(i), expected_RHS[i], 1.0e-2);
            for (unsigned int j = 0; j < 2; ++j) {
                KRATOS_EXPECT_NEAR(LHS(i,j), expected_LHS[i*2+j], 1.0e-4);
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(ThermalFace3D3N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test condition
        Model model;
        ModelPart &test_model_part = model.CreateModelPart("TestModelPart");
        ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(test_model_part);

        // Set the condition properties
        auto p_cond_prop = test_model_part.pGetProperties(0);
        p_cond_prop->SetValue(EMISSIVITY, 1.0);
        p_cond_prop->SetValue(AMBIENT_TEMPERATURE, 293.0);
        p_cond_prop->SetValue(CONVECTION_COEFFICIENT, 20.0);

        // Condition creation
        test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3};
        test_model_part.CreateNewCondition("ThermalFace3D3N", 1, cond_nodes, p_cond_prop);

        // Set the face heat flux
        for (auto &i_node : test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(TEMPERATURE) = 400.0;
            i_node.FastGetSolutionStepValue(FACE_HEAT_FLUX) = 200.0;
            // i_node.FastGetSolutionStepValue(TEMPERATURE) = i_node.X() * 400.0;
            // i_node.FastGetSolutionStepValue(FACE_HEAT_FLUX) = 200.0 + i_node.X() * 100.0;
        }

        // Test condition
        Condition::Pointer p_condition = test_model_part.pGetCondition(1);
        Vector RHS = ZeroVector(3);
        Matrix LHS = ZeroMatrix(3,3);
        const auto& rConstProcessInfoRef = test_model_part.GetProcessInfo();
        p_condition->CalculateLocalSystem(LHS, RHS, rConstProcessInfoRef);

        std::vector<double> expected_RHS = {-495.606, -495.606, -495.606};
        std::vector<double> expected_LHS = {2.87627, 1.43813, 1.43813,
                                            1.43813, 2.87627, 1.43813,
                                            1.43813, 1.43813, 2.87627};
        for (unsigned int i = 0; i < 3; ++i) {
            KRATOS_EXPECT_NEAR(RHS(i), expected_RHS[i], 1.0e-3);
            for (unsigned int j = 0; j < 3; ++j) {
                KRATOS_EXPECT_NEAR(LHS(i,j), expected_LHS[i*3+j], 1.0e-5);
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(ThermalFace3D4N, KratosConvectionDiffusionFastSuite)
    {
        // Create the test condition
        Model model;
        ModelPart &test_model_part = model.CreateModelPart("TestModelPart");
        ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(test_model_part);

        // Set the condition properties
        auto p_cond_prop = test_model_part.pGetProperties(0);
        p_cond_prop->SetValue(EMISSIVITY, 1.0);
        p_cond_prop->SetValue(AMBIENT_TEMPERATURE, 293.0);
        p_cond_prop->SetValue(CONVECTION_COEFFICIENT, 20.0);

        // Condition creation
        test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        test_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
        test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4};
        test_model_part.CreateNewCondition("ThermalFace3D4N", 1, cond_nodes, p_cond_prop);

        // Set the face heat flux
        for (auto &i_node : test_model_part.Nodes()) {
            i_node.FastGetSolutionStepValue(TEMPERATURE) = 400.0;
            i_node.FastGetSolutionStepValue(FACE_HEAT_FLUX) = 200.0;
            // i_node.FastGetSolutionStepValue(TEMPERATURE) = i_node.X() * 400.0;
            // i_node.FastGetSolutionStepValue(FACE_HEAT_FLUX) = 200.0 + i_node.X() * 100.0;
        }

        // Test condition
        Condition::Pointer p_condition = test_model_part.pGetCondition(1);
        Vector RHS = ZeroVector(4);
        Matrix LHS = ZeroMatrix(4,4);
        const auto& rConstProcessInfoRef = test_model_part.GetProcessInfo();

        p_condition->CalculateLocalSystem(LHS, RHS, rConstProcessInfoRef);

        std::vector<double> expected_RHS = {-743.41, -743.41, -743.41, -743.41};
        std::vector<double> expected_LHS = {3.83502, 1.91751, 0.958756, 1.91751,
                                            1.91751, 3.83502, 1.91751, 0.958756,
                                            0.958756, 1.91751, 3.83502, 1.91751,
                                            1.91751, 0.958756, 1.91751, 3.83502};
        for (unsigned int i = 0; i < 4; ++i) {
            KRATOS_EXPECT_NEAR(RHS(i), expected_RHS[i], 1.0e-3);
            for (unsigned int j = 0; j < 4; ++j) {
                KRATOS_EXPECT_NEAR(LHS(i,j), expected_LHS[i*4+j], 1.0e-5);
            }
        }
    }

} // namespace Kratos::Testing.
