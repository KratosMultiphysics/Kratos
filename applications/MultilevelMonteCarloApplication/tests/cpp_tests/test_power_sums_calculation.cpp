//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/variables.h"

// Application includes
#include "custom_statistics/power_sums_statistics.h"
#include "multilevel_monte_carlo_application_variables.h"

namespace Kratos {
	namespace Testing {

        /**
	     * Auxiliar function to generate a triangular element to be tested.
	     */
        void GenerateModelPartToTestPowerSums(
            ModelPart& rModelPart) {

            // Variables addition
            rModelPart.AddNodalSolutionStepVariable(PRESSURE);

            // Element creation
            rModelPart.CreateNewNode(1, 1.0, 1.0, 0.0);
            rModelPart.CreateNewNode(2, 2.0, 1.0, 0.0);
            rModelPart.CreateNewNode(3, 1.5, 2.0, 0.0);
            std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
            Properties::Pointer p_properties = rModelPart.CreateNewProperties(0);
            rModelPart.CreateNewElement("Element2D3N", 1, elem_nodes, p_properties);

            // Add nodal data
            auto nodes_begin = rModelPart.NodesBegin();
            for (unsigned int i_node = 1; i_node <= rModelPart.NumberOfNodes(); ++i_node){
                auto it_node = nodes_begin + (i_node-1);
                it_node->AddDof(POWER_SUM_1);
                it_node->AddDof(POWER_SUM_2);
                it_node->AddDof(POWER_SUM_3);
                it_node->AddDof(POWER_SUM_4);
                it_node->AddDof(POWER_SUM_5);
                it_node->AddDof(POWER_SUM_6);
                it_node->AddDof(POWER_SUM_7);
                it_node->AddDof(POWER_SUM_8);
                it_node->AddDof(POWER_SUM_9);
                it_node->AddDof(POWER_SUM_10);
                it_node->SetValue(POWER_SUM_1,1.0 * i_node);
                it_node->SetValue(POWER_SUM_2,1.0 * i_node);
                it_node->SetValue(POWER_SUM_3,1.0 * i_node);
                it_node->SetValue(POWER_SUM_4,1.0 * i_node);
                it_node->SetValue(POWER_SUM_5,1.0 * i_node);
                it_node->SetValue(POWER_SUM_6,1.0 * i_node);
                it_node->SetValue(POWER_SUM_7,1.0 * i_node);
                it_node->SetValue(POWER_SUM_8,1.0 * i_node);
                it_node->SetValue(POWER_SUM_9,1.0 * i_node);
                it_node->SetValue(POWER_SUM_10,1.0 * i_node);
                it_node->FastGetSolutionStepValue(PRESSURE) = 10.0;
            }
        }

	    /**
	     * Checks the time average pressure utility.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(CalculationPowerSums, KratosMultilevelMonteCarloApplicationFastSuite)
		{
            // Create a test element inside a modelpart
            Model model;
            ModelPart& model_part = model.CreateModelPart("Main", 3);
            GenerateModelPartToTestPowerSums(model_part);
            Element::Pointer p_element = model_part.pGetElement(1);

            // Initialize the element
            p_element->Initialize();

            // Generate parameters
            Parameters parameters = Parameters(R"(
            {
                "reference_variable_name" : "PRESSURE"
            })"
            );

            // Call the power sums statistics
            PowerSumsStatistics(model_part,parameters).ExecuteFinalizeSolutionStep();

            // Check computed values over the nodes
            auto nodes_begin = model_part.NodesBegin();
            for (unsigned int i_node = 0; i_node < model_part.NumberOfNodes(); ++i_node){
                auto it_node = nodes_begin + (i_node);
                double S1_value = it_node->GetValue(POWER_SUM_1);
                double S2_value = it_node->GetValue(POWER_SUM_2);
                double S3_value = it_node->GetValue(POWER_SUM_3);
                double S4_value = it_node->GetValue(POWER_SUM_4);
                double S5_value = it_node->GetValue(POWER_SUM_5);
                double S6_value = it_node->GetValue(POWER_SUM_6);
                double S7_value = it_node->GetValue(POWER_SUM_7);
                double S8_value = it_node->GetValue(POWER_SUM_8);
                double S9_value = it_node->GetValue(POWER_SUM_9);
                double S10_value = it_node->GetValue(POWER_SUM_10);

                KRATOS_CHECK_NEAR(S1_value, (11 + i_node), 1e-10);
                KRATOS_CHECK_NEAR(S2_value, (101 + i_node), 1e-10);
                KRATOS_CHECK_NEAR(S3_value, (1001 + i_node), 1e-10);
                KRATOS_CHECK_NEAR(S4_value, (10001 + i_node), 1e-10);
                KRATOS_CHECK_NEAR(S5_value, (100001 + i_node), 1e-10);
                KRATOS_CHECK_NEAR(S6_value, (1000001 + i_node), 1e-10);
                KRATOS_CHECK_NEAR(S7_value, (10000001 + i_node), 1e-10);
                KRATOS_CHECK_NEAR(S8_value, (100000001 + i_node), 1e-10);
                KRATOS_CHECK_NEAR(S9_value, (1000000001 + i_node), 1e-10);
                KRATOS_CHECK_NEAR(S10_value, (10000000001 + i_node), 1e-10);
            }
        }

    } // namespace Testing
}  // namespace Kratos.
