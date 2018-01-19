//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_flags.h"
#include "includes/mapping_variables.h"
#include "contact_structural_mechanics_application.h"
#include "contact_structural_mechanics_application_variables.h"
#include "includes/gid_io.h"
#include "utilities/variable_utils.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3>                                                    NodeType;

        void GiDIOGapDebug(ModelPart& ThisModelPart)
        {
            GidIO<> gid_io("TEST_WEIGHTED_GAP", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditionsOnly);
            const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);

            gid_io.InitializeMesh(label);
            gid_io.WriteMesh(ThisModelPart.GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(label, ThisModelPart.GetMesh());
            gid_io.WriteNodalFlags(ACTIVE, "ACTIVE", ThisModelPart.Nodes(), label);
            gid_io.WriteNodalFlags(SLAVE, "SLAVE", ThisModelPart.Nodes(), label);
            gid_io.WriteNodalFlags(ISOLATED, "ISOLATED", ThisModelPart.Nodes(), label);
            gid_io.WriteNodalResults(WEIGHTED_GAP, ThisModelPart.Nodes(), label, 0);
            gid_io.WriteNodalResults(NORMAL, ThisModelPart.Nodes(), label, 0);
        }
        
        /**
         * This method can be used to create a plane/cylinder condition set
         */
        void CreatePlaneCilynderProblem(
            ModelPart& ThisModelPart,
            const int NumberOfDivisions,
            const double Lenght,
            const double Radius,
            const double Angle,
            const double Slope = 0.0
            )
        {
            Properties::Pointer p_cond_prop = ThisModelPart.pGetProperties(0);
            
            // Creating the base geometry
            std::size_t id_node = 0;
            const double dx = Lenght/static_cast<double>(NumberOfDivisions);
            for (auto i = 0; i < NumberOfDivisions + 1; ++i){
                const double x = dx * i;
                const double y = Slope * dx * i;
                id_node++;
                NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(id_node, x , y , 0.0);
                p_node_1->Set(SLAVE, true);
                p_node_1->Set(MASTER, false);
                p_node_1->Set(ACTIVE, true);
                id_node++;
                NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(id_node, x , y , 1.0);
                p_node_2->Set(SLAVE, true);
                p_node_2->Set(MASTER, false);
                p_node_2->Set(ACTIVE, true);
            }
            
            std::size_t id_cond = 0;
            for (auto i = 0; i < NumberOfDivisions; i++)
            {
                id_cond++;
                std::vector<NodeType::Pointer> condition_nodes (4);
                condition_nodes[0] = ThisModelPart.pGetNode((2 * i)+1);
                condition_nodes[1] = ThisModelPart.pGetNode((2 * i)+2);
                condition_nodes[2] = ThisModelPart.pGetNode((2 * i)+4);
                condition_nodes[3] = ThisModelPart.pGetNode((2 * i)+3);
                Quadrilateral3D4 <NodeType> quad( condition_nodes);
                
                Condition::Pointer pcond = ThisModelPart.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", id_cond, quad, p_cond_prop);
                pcond->Set(SLAVE, true);
                pcond->Set(ISOLATED, false);
                pcond->Set(MASTER, false);
            }
            
            // Creating the base circle
            const double dtheta = Angle/static_cast<double>(NumberOfDivisions);
            for (auto i = 0; i < NumberOfDivisions + 1; ++i){
                const double x = Radius * std::sin(i * dtheta);
                const double y = Radius * (1.0 - std::cos(i * dtheta));
                id_node++;
                NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(id_node, x, y , 0.0);
                p_node_1->Set(SLAVE, false);
                p_node_1->Set(MASTER, true);
                p_node_1->Set(ACTIVE, false);
                id_node++;
                NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(id_node, x, y , 1.0);
                p_node_2->Set(SLAVE, false);
                p_node_2->Set(MASTER, true);
                p_node_2->Set(ACTIVE, false);
            }
            
            // Adding map
            IndexSet this_set;
            for (auto i = 0; i < NumberOfDivisions; i++)
            {
                id_cond++;
                this_set.AddId(id_cond);
                std::vector<NodeType::Pointer> condition_nodes (4);
                condition_nodes[3] = ThisModelPart.pGetNode((2 * (i + NumberOfDivisions + 1)+1));
                condition_nodes[2] = ThisModelPart.pGetNode((2 * (i + NumberOfDivisions + 1)+2));
                condition_nodes[1] = ThisModelPart.pGetNode((2 * (i + NumberOfDivisions + 1)+4));
                condition_nodes[0] = ThisModelPart.pGetNode((2 * (i + NumberOfDivisions + 1)+3));
                Quadrilateral3D4 <NodeType> quad( condition_nodes);
                
                Condition::Pointer pcond = ThisModelPart.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", id_cond, quad, p_cond_prop);
                pcond->Set(SLAVE, false);
                pcond->Set(MASTER, true);
            }
            
            // We set the database
            for (auto& i_cond : ThisModelPart.Conditions()) {
                if (i_cond.Is(SLAVE) == true)
                    i_cond.SetValue(INDEX_SET, Kratos::make_shared<IndexSet>(this_set));
            }
        }
        
        /** 
        * Checks the correct work of the weighted gap computation
        * Test 1
        */

        KRATOS_TEST_CASE_IN_SUITE(TestWeightedGap1, ContactStructuralApplicationFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            
            this_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            this_model_part.AddNodalSolutionStepVariable(NORMAL);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes
            const int number_of_divisions = 4;
            const double lenght = 4.0;
            const double radius = 6.0;
            const double angle = Globals::Pi/6;
            
            // We create our problem
            CreatePlaneCilynderProblem(this_model_part, number_of_divisions, lenght, radius, angle);
            
            // We compute the normals
            MortarUtilities::ComputeNodesMeanNormalModelPart(this_model_part);
            
            // We compute the explicit contribution
            VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, this_model_part.Nodes());
            for (auto& id_cond : this_model_part.Conditions())
                if (id_cond.Is(ACTIVE) == true)
                    id_cond.AddExplicitContribution(process_info);
                
            // DEBUG         
            GiDIOGapDebug(this_model_part);
            
//             const double tolerance = 1.0e-4;
//             KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_1->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_1->X(), 2) + std::pow(p_node_1->Y(), 2))), tolerance);
//             KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_2->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_2->X(), 2) + std::pow(p_node_2->Y(), 2))), tolerance);
//             KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_3->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_3->X(), 2) + std::pow(p_node_3->Y(), 2))), tolerance);
//             KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_4->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_4->X(), 2) + std::pow(p_node_4->Y(), 2))), tolerance);
        }
        
    } // namespace Testing
}  // namespace Kratos.
