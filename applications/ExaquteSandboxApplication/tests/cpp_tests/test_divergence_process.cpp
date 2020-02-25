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
#include "custom_processes/calculate_divergence_process.h"
#include "exaqute_sandbox_application_variables.h"

namespace Kratos {
	namespace Testing {

        /**
	     * Auxiliar function to generate a triangular element to be tested.
	     */
        void GenerateModelPartToTestDivergenceProcess(
            ModelPart& rModelPart) {

            // Variables addition
            rModelPart.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

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
                it_node->AddDof(VELOCITY_X);
                it_node->AddDof(VELOCITY_Y);
                it_node->AddDof(VELOCITY_Z);
                it_node->FastGetSolutionStepValue(VELOCITY_X) = 1.0 * i_node;
                it_node->FastGetSolutionStepValue(VELOCITY_Y) = 2.0 * i_node;
                it_node->FastGetSolutionStepValue(VELOCITY_Z) = 0.0 * i_node;
            }
        }

	    /**
	     * Checks the time average divergence utility.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(CalculationDivergence, ExaquteSandboxApplicationFastSuite)
		{
            // Create a test element inside a modelpart
            Model model;
            ModelPart& model_part = model.CreateModelPart("Main", 3);
            GenerateModelPartToTestDivergenceProcess(model_part);
            Element::Pointer p_element = model_part.pGetElement(1);

            // Call the divergence time average process
            CalculateDivergenceProcess(model_part).ExecuteInitialize();
            CalculateDivergenceProcess(model_part).ExecuteBeforeOutputStep();

            // Check computed values over the elements
            auto elements_begin = model_part.ElementsBegin();
            for (unsigned int i_elem = 0; i_elem < model_part.NumberOfElements(); ++i_elem){
                auto it_elem = elements_begin + i_elem;
                double divergence_value = it_elem->GetValue(DIVERGENCE);
                double velocity_seminorm_value = it_elem->GetValue(VELOCITY_H1_SEMINORM);
                KRATOS_CHECK_NEAR(divergence_value, 8.0, 1e-10);
                KRATOS_CHECK_NEAR(velocity_seminorm_value, 8.125, 1e-10);
            }
        }

    } // namespace Testing
}  // namespace Kratos.
