// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes

// External includes

// Project includes
#include "geometries/triangle_3d_3.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"
#include "includes/gid_io.h"
#include "containers/model.h"

/* Processes */
#include "custom_processes/alm_variables_calculation_process.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

        void GiDIODebugALMVariables(ModelPart& ThisModelPart)
        {
            GidIO<> gid_io("TEST_ALM", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditionsOnly);
            const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);

            gid_io.InitializeMesh(label);
            gid_io.WriteMesh(ThisModelPart.GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(label, ThisModelPart.GetMesh());
            gid_io.WriteNodalResults(NODAL_H, ThisModelPart.Nodes(), label, 0);
        }

        void Create3DConditionsGeometry(ModelPart& ThisModelPart, const std::string& rConditionName)
        {
            Properties::Pointer p_cond_prop = ThisModelPart.CreateNewProperties(0);
            p_cond_prop->SetValue(YOUNG_MODULUS, 100.0);

            // First we create the nodes
            NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(2, 1.0 , 0.2 , 0.0);
            NodeType::Pointer p_node_3 = ThisModelPart.CreateNewNode(3, 0.9 , 1.3 , 0.0);
            NodeType::Pointer p_node_4 = ThisModelPart.CreateNewNode(4, 0.1 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = ThisModelPart.CreateNewNode(5, 2.2 , 0.1 , 0.0);
            NodeType::Pointer p_node_6 = ThisModelPart.CreateNewNode(6, 2.3 , 1.0 , 0.0);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{element_nodes_0} );

            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{element_nodes_1} );

            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_2( PointerVector<NodeType>{element_nodes_2} );

            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_3( PointerVector<NodeType>{element_nodes_3} );

            Condition::Pointer p_cond_0 = ThisModelPart.CreateNewCondition(rConditionName, 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond_1 = ThisModelPart.CreateNewCondition(rConditionName, 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond_2 = ThisModelPart.CreateNewCondition(rConditionName, 3, triangle_2, p_cond_prop);
            Condition::Pointer p_cond_3 = ThisModelPart.CreateNewCondition(rConditionName, 4, triangle_3, p_cond_prop);
        }

        /**
        * Checks the correct work of the ALM variables process
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(ALMVariablesProcess, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& this_model_part = this_model.CreateModelPart("Main", 2);

            this_model_part.AddNodalSolutionStepVariable(NODAL_H);

            auto& r_process_info = this_model_part.GetProcessInfo();
            r_process_info.SetValue(DOMAIN_SIZE, 3);
            r_process_info.SetValue(STEP, 1);
            r_process_info.SetValue(NL_ITERATION_NUMBER, 1);

            Create3DConditionsGeometry(this_model_part, "SurfaceCondition3D3N");

            // Assign NodalH and flags
            for (auto& r_node : this_model_part.Nodes()) {
                r_node.FastGetSolutionStepValue(NODAL_H) = static_cast<double>(r_node.Id()); // Adding some randomness
//                 r_node.FastGetSolutionStepValue(NODAL_H) = 1.0;
                r_node.Set(SLAVE);
                r_node.Set(MASTER);
            }
            for (auto& r_cond : this_model_part.Conditions()) {
                r_cond.Set(SLAVE);
                r_cond.Set(MASTER);
            }

            // Compute ALM variables
            ALMVariablesCalculationProcess alm_process(this_model_part);
            alm_process.Execute();

            const double tolerance = 1.0e-4;
            const double initial_penalty = this_model_part.GetProcessInfo().GetValue(INITIAL_PENALTY);
            const double scale_factor = this_model_part.GetProcessInfo().GetValue(SCALE_FACTOR);

//             // DEBUG
//             GiDIODebugALMVariables(this_model_part);

            KRATOS_CHECK_LESS_EQUAL(std::abs((initial_penalty - 305.858)/initial_penalty), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs((scale_factor - 305.858)/scale_factor), tolerance);
        }
    } // namespace Testing
}  // namespace Kratos.
