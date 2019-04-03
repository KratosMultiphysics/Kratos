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
#include "containers/model.h"
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

        KRATOS_TEST_CASE_IN_SUITE(AALMProcess1, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);
            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            r_model_part.AddNodalSolutionStepVariable(NODAL_H);
            r_model_part.AddNodalSolutionStepVariable(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);
            
            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            double& penalty_parameter = process_info[INITIAL_PENALTY];
            penalty_parameter = 1.0e7;
            double& max_gap_factor = process_info[MAX_GAP_FACTOR];
            max_gap_factor = 1.0;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = r_model_part.CreateNewNode(0,0.0,0.0,0.0);
            p_node_1->SetValue(NODAL_AREA, 1.0);
            p_node_1->FastGetSolutionStepValue(NODAL_H) = 0.1;
            p_node_1->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.05;
            p_node_1->FastGetSolutionStepValue(WEIGHTED_GAP, 1) = -0.1;
            
            AALMAdaptPenaltyValueProcess process = AALMAdaptPenaltyValueProcess(r_model_part);
            process.Execute();
            
//             // DEBUG
//             KRATOS_WATCH(p_node_1->GetValue(INITIAL_PENALTY))
            
            const double tolerance = 1.0e-6;
            KRATOS_CHECK_NEAR(p_node_1->GetValue(INITIAL_PENALTY), 0.2 * penalty_parameter, tolerance);
        }
        
    } // namespace Testing
}  // namespace Kratos.
