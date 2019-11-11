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
#include "utilities/variable_utils.h"

// Application includes
#include "custom_processes/weighted_velocity_calculation_process.h"
#include "exaqute_sandbox_application_variables.h"

namespace Kratos {
	namespace Testing {

        /**
	     * Auxiliar function to generate a triangular element to be tested.
	     */
        void GenerateModelPartToTestVelocity(
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

            // Initialize average velocity
            auto& r_nodes_array = rModelPart.Nodes();
            VariableUtils().SetNonHistoricalVariableToZero(VELOCITY_WEIGHTED, r_nodes_array);

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
	     * Checks the time average velocity utility.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(CalculationWeightedVelocity, ExaquteSandboxApplicationFastSuite)
		{
            // Create a test element inside a modelpart
            Model model;
            ModelPart& model_part = model.CreateModelPart("Main", 3);
            GenerateModelPartToTestVelocity(model_part);
            Element::Pointer p_element = model_part.pGetElement(1);

            // Initialize the element
            p_element->Initialize();

            // Call the velocity time average process
            WeightedVelocityCalculationProcess(model_part).ExecuteFinalizeSolutionStep();

            // Check computed values over the nodes
            auto nodes_begin = model_part.NodesBegin();
            for (unsigned int i_node = 0; i_node < model_part.NumberOfNodes(); ++i_node){
                auto it_node = nodes_begin + (i_node);
                array_1d<double, 3> velocity_value = it_node->GetValue(VELOCITY_WEIGHTED);
                KRATOS_CHECK_NEAR(velocity_value[0], 0.2*(i_node+1), 1e-10);
                KRATOS_CHECK_NEAR(velocity_value[1], 0.4*(i_node+1), 1e-10);
                KRATOS_CHECK_NEAR(velocity_value[2], 0.0*(i_node+1), 1e-10);
            }

        }

    } // namespace Testing
}  // namespace Kratos.
