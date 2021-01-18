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
#include "custom_processes/weighted_divergence_calculation_process.h"
#include "exaqute_sandbox_application_variables.h"

namespace Kratos {
	namespace Testing {

        /**
	     * Auxiliar function to generate a triangular element to be tested.
	     */
        void GenerateModelPartToTestDivergence(
            ModelPart& rModelPart) {

            // Variables addition
            rModelPart.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            const int buffer_size = 2;
            const double delta_time = 0.1;
            rModelPart.SetBufferSize(buffer_size);
            rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
            rModelPart.GetProcessInfo().SetValue(TIME, 0.6);
            rModelPart.GetProcessInfo().SetValue(END_TIME, 1.0);
            rModelPart.CloneTimeStep(rModelPart.GetProcessInfo().GetValue(TIME) + delta_time);
            rModelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(TIME, 0.6);
            rModelPart.GetProcessInfo().SetValue(END_TIME, 1.0);

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

            // Add element data
            auto elements_begin = rModelPart.ElementsBegin();
            for (unsigned int i_elem = 0; i_elem < rModelPart.NumberOfElements(); ++i_elem){
                auto it_elem = elements_begin + i_elem;
                it_elem->SetValue(AVERAGED_DIVERGENCE,0.5);
                it_elem->SetValue(VELOCITY_H1_SEMINORM,0.5);
            }
        }

	    /**
	     * Checks the time average divergence utility.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(CalculationWeightedDivergence, ExaquteSandboxApplicationFastSuite)
		{
            // Create a test element inside a modelpart
            Model model;
            ModelPart& model_part = model.CreateModelPart("Main", 3);
            GenerateModelPartToTestDivergence(model_part);
            Element::Pointer p_element = model_part.pGetElement(1);

            // Initialize the element
            p_element->Initialize();

            // Call the divergence time average process
            WeightedDivergenceCalculationProcess(model_part).ExecuteFinalizeSolutionStep();

            // Check computed values over the elements
            auto elements_begin = model_part.ElementsBegin();
            for (unsigned int i_elem = 0; i_elem < model_part.NumberOfElements(); ++i_elem){
                auto it_elem = elements_begin + i_elem;
                double divergence_value = it_elem->GetValue(AVERAGED_DIVERGENCE);
                double velocity_seminorm_value = it_elem->GetValue(VELOCITY_H1_SEMINORM);
                KRATOS_CHECK_NEAR(divergence_value, 1.3416407865, 1e-10);
                KRATOS_CHECK_NEAR(velocity_seminorm_value, 1.3509256086, 1e-10);
            }
        }

    } // namespace Testing
}  // namespace Kratos.
