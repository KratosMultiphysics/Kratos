// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "contact_structural_mechanics_application_variables.h"

/* Processes */
#include "custom_processes/aalm_adapt_penalty_value_process.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Point                                                     PointType;
        typedef Node<3>                                                    NodeType;
        typedef Geometry<NodeType>                                 GeometryNodeType;
        typedef Geometry<PointType>                               GeometryPointType;

        /** 
        * Checks the correct work of the AALM  dynamic penalty process
        */

        KRATOS_TEST_CASE_IN_SUITE(TestAALMProcess1, ContactStructuralApplicationFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(3);
            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            this_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            this_model_part.AddNodalSolutionStepVariable(NODAL_H);
            this_model_part.AddNodalSolutionStepVariable(NORMAL_CONTACT_STRESS);
            
            double& penalty_parameter = this_model_part.GetProcessInfo()[INITIAL_PENALTY];
            penalty_parameter = 1.0e7;
            double& max_gap_factor = this_model_part.GetProcessInfo()[MAX_GAP_FACTOR];
            max_gap_factor = 1.0;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(0,0.0,0.0,0.0);
            p_node_1->FastGetSolutionStepValue(NODAL_H) = 0.1;
            p_node_1->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.05;
            p_node_1->FastGetSolutionStepValue(WEIGHTED_GAP, 1) = -0.1;
            
            AALMAdaptPenaltyValueProcess process = AALMAdaptPenaltyValueProcess(this_model_part);
            process.Execute();
            
//             // DEBUG
//             KRATOS_WATCH(p_node_1->GetValue(INITIAL_PENALTY))
            
            const double tolerance = 1.0e-6;
            KRATOS_CHECK_NEAR(p_node_1->GetValue(INITIAL_PENALTY), 0.2 * penalty_parameter, tolerance);
        }
        
    } // namespace Testing
}  // namespace Kratos.
