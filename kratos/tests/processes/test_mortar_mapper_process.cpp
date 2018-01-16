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
#include "includes/gid_io.h"

/* Processes */
#include "processes/simple_mortar_mapper_process.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3>                                                    NodeType;

        void GiDIODebugMapper(ModelPart& ThisModelPart)
        {
            GidIO<> gid_io("TEST_MAPPER", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditionsOnly);
            const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);

            gid_io.InitializeMesh(label);
            gid_io.WriteMesh(ThisModelPart.GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(label, ThisModelPart.GetMesh());
            gid_io.WriteNodalResults(TEMPERATURE, ThisModelPart.Nodes(), label, 0);
            gid_io.WriteNodalResults(NORMAL, ThisModelPart.Nodes(), label, 0);
        }
        
        /** 
        * Checks the correct work of the simple mortar mapper
        * Test 1
        */

        KRATOS_TEST_CASE_IN_SUITE(TestSimpleMortarMapper1, KratosCoreMortarMapperFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            
            this_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
            this_model_part.AddNodalSolutionStepVariable(NORMAL);
            
            Properties::Pointer p_cond_prop = this_model_part.pGetProperties(0);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);
            
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_3;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_1;
            Triangle3D3 <NodeType> triangle_0( condition_nodes_0 );
            
            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( condition_nodes_1 );
            
            Condition::Pointer p_cond_0 = this_model_part.CreateNewCondition("Condition3D", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond_1 = this_model_part.CreateNewCondition("Condition3D", 2, triangle_1, p_cond_prop);
            
            // Adding map
            IndexSet this_set;
            this_set.AddId(2);
            p_cond_0->SetValue(INDEX_SET, boost::make_shared<IndexSet>(this_set));
            
            // Setting flags
            // SLAVE
            p_node_1->Set(SLAVE, true);
            p_node_1->Set(MASTER, false);
            p_node_2->Set(SLAVE, true);
            p_node_2->Set(MASTER, false);
            p_node_3->Set(SLAVE, true);
            p_node_3->Set(MASTER, false);
            p_cond_0->Set(SLAVE, true);
            p_cond_0->Set(MASTER, false);
            // MASTER
            p_node_4->Set(SLAVE, false);
            p_node_4->Set(MASTER, true);
            p_node_5->Set(SLAVE, false);
            p_node_5->Set(MASTER, true);
            p_node_6->Set(SLAVE, false);
            p_node_6->Set(MASTER, true);
            p_cond_1->Set(SLAVE, false);
            p_cond_1->Set(MASTER, true);
            
            // We compute the normals
            MortarUtilities::ComputeNodesMeanNormalModelPart(this_model_part);
            
            p_node_4->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_4->X(), 2) + std::pow(p_node_4->Y(), 2);
            p_node_5->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_5->X(), 2) + std::pow(p_node_5->Y(), 2);
            p_node_6->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_6->X(), 2) + std::pow(p_node_6->Y(), 2);
                         
            typedef SimpleMortarMapperProcess<3, 3, Variable<double>, Historical> MapperType;
            MapperType process = MapperType(this_model_part, TEMPERATURE);
            process.Execute();
            
            // DEBUG         
//             GiDIODebugMapper(this_model_part);
            
            const double tolerance = 1.0e-4;
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_1->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_1->X(), 2) + std::pow(p_node_1->Y(), 2))), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_2->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_2->X(), 2) + std::pow(p_node_2->Y(), 2))), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_3->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_3->X(), 2) + std::pow(p_node_3->Y(), 2))), tolerance);
        }
        
        /** 
        * Checks the correct work of the simple mortar mapper
        * Test 2
        */

        KRATOS_TEST_CASE_IN_SUITE(TestSimpleMortarMapper2, KratosCoreMortarMapperFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            
            this_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
            this_model_part.AddNodalSolutionStepVariable(NORMAL);
            
            Properties::Pointer p_cond_prop = this_model_part.pGetProperties(0);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.01);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.01);
            
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 0.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 1.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_7 = this_model_part.CreateNewNode(7, 1.0 , 1.0 , 0.02);
            NodeType::Pointer p_node_8 = this_model_part.CreateNewNode(8, 0.0 , 1.0 , 0.02);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (4);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            condition_nodes_0[3] = p_node_4;
            Quadrilateral3D4 <NodeType> quad_0( condition_nodes_0 );
            
            std::vector<NodeType::Pointer> condition_nodes_1 (4);
            condition_nodes_1[0] = p_node_8;
            condition_nodes_1[1] = p_node_7;
            condition_nodes_1[2] = p_node_6;
            condition_nodes_1[3] = p_node_5;
            Quadrilateral3D4 <NodeType> quad_1( condition_nodes_1 );
            
            Condition::Pointer p_cond_0 = this_model_part.CreateNewCondition("Condition3D4N", 1, quad_0, p_cond_prop);
            Condition::Pointer p_cond_1 = this_model_part.CreateNewCondition("Condition3D4N", 2, quad_1, p_cond_prop);
            
            // Adding map
            IndexSet this_set;
            this_set.AddId(2);
            p_cond_0->SetValue(INDEX_SET, boost::make_shared<IndexSet>(this_set));
            
            // Setting flags
            // SLAVE
            p_node_1->Set(SLAVE, true);
            p_node_1->Set(MASTER, false);
            p_node_2->Set(SLAVE, true);
            p_node_2->Set(MASTER, false);
            p_node_3->Set(SLAVE, true);
            p_node_3->Set(MASTER, false);
            p_node_4->Set(SLAVE, true);
            p_node_4->Set(MASTER, false);
            p_cond_0->Set(SLAVE, true);
            p_cond_0->Set(MASTER, false);
            // MASTER
            p_node_5->Set(SLAVE, false);
            p_node_5->Set(MASTER, true);
            p_node_6->Set(SLAVE, false);
            p_node_6->Set(MASTER, true);
            p_node_7->Set(SLAVE, false);
            p_node_7->Set(MASTER, true);
            p_node_8->Set(SLAVE, false);
            p_node_8->Set(MASTER, true);
            p_cond_1->Set(SLAVE, false);
            p_cond_1->Set(MASTER, true);
            
            // We compute the normals
            MortarUtilities::ComputeNodesMeanNormalModelPart(this_model_part);
            
            p_node_5->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_5->X(), 2) + std::pow(p_node_5->Y(), 2);
            p_node_6->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_6->X(), 2) + std::pow(p_node_6->Y(), 2);
            p_node_7->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_7->X(), 2) + std::pow(p_node_7->Y(), 2);
            p_node_8->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_8->X(), 2) + std::pow(p_node_8->Y(), 2);
                         
            typedef SimpleMortarMapperProcess<3, 4, Variable<double>, Historical> MapperType;
            MapperType process = MapperType(this_model_part, TEMPERATURE);
            process.Execute();
            
            // DEBUG         
//             GiDIODebugMapper(this_model_part);
            
            const double tolerance = 1.0e-4;
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_1->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_1->X(), 2) + std::pow(p_node_1->Y(), 2))), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_2->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_2->X(), 2) + std::pow(p_node_2->Y(), 2))), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_3->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_3->X(), 2) + std::pow(p_node_3->Y(), 2))), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_4->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_4->X(), 2) + std::pow(p_node_4->Y(), 2))), tolerance);
        }
        
        /** 
        * Checks the correct work of the simple mortar mapper
        * Test 3
        */

        KRATOS_TEST_CASE_IN_SUITE(TestSimpleMortarMapper3, KratosCoreMortarMapperFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            
            this_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
            this_model_part.AddNodalSolutionStepVariable(NORMAL);
            
            Properties::Pointer p_cond_prop = this_model_part.pGetProperties(0);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.000, 0.000, 0.000);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.000, 1.000, 1.000);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.000, 0.000, 1.000);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 1.000, 1.000, 0.000);
            
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 1.001, 0.000, 1.000);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 1.001, 1.000, 0.000);
            NodeType::Pointer p_node_7 = this_model_part.CreateNewNode(7, 0.001, 0.000, 0.000);
            NodeType::Pointer p_node_8 = this_model_part.CreateNewNode(8, 1.001, 1.000, 1.000);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_0[0] = p_node_4;
            condition_nodes_0[1] = p_node_3;
            condition_nodes_0[2] = p_node_1;
            Triangle3D3 <NodeType> triangle_0( condition_nodes_0 );
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_2;
            condition_nodes_1[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_1( condition_nodes_1 );
            
            std::vector<NodeType::Pointer> condition_nodes_2 (3);
            std::vector<NodeType::Pointer> condition_nodes_3 (3);
            condition_nodes_2[0] = p_node_7;
            condition_nodes_2[1] = p_node_5;
            condition_nodes_2[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_3( condition_nodes_2 );
            condition_nodes_3[0] = p_node_5;
            condition_nodes_3[1] = p_node_8;
            condition_nodes_3[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_4( condition_nodes_3 );
            
            Condition::Pointer p_cond_0 = this_model_part.CreateNewCondition("Condition3D", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond_1 = this_model_part.CreateNewCondition("Condition3D", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond_2 = this_model_part.CreateNewCondition("Condition3D", 3, triangle_3, p_cond_prop);
            Condition::Pointer p_cond_3 = this_model_part.CreateNewCondition("Condition3D", 4, triangle_4, p_cond_prop);
            
            // Adding map
            IndexSet this_set0, this_set1;
            this_set0.AddId(3);
            this_set0.AddId(4);
            this_set1.AddId(3);
            this_set1.AddId(4);
            p_cond_0->SetValue(INDEX_SET, boost::make_shared<IndexSet>(this_set0));
            p_cond_1->SetValue(INDEX_SET, boost::make_shared<IndexSet>(this_set1));
            
            // Setting flags
            // SLAVE
            p_node_1->Set(SLAVE, true);
            p_node_1->Set(MASTER, false);
            p_node_2->Set(SLAVE, true);
            p_node_2->Set(MASTER, false);
            p_node_3->Set(SLAVE, true);
            p_node_3->Set(MASTER, false);
            p_node_4->Set(SLAVE, true);
            p_node_4->Set(MASTER, false);
            p_cond_0->Set(SLAVE, true);
            p_cond_0->Set(MASTER, false);
            p_cond_1->Set(SLAVE, true);
            p_cond_1->Set(MASTER, false);
            // MASTER
            p_node_5->Set(SLAVE, false);
            p_node_5->Set(MASTER, true);
            p_node_6->Set(SLAVE, false);
            p_node_6->Set(MASTER, true);
            p_node_7->Set(SLAVE, false);
            p_node_7->Set(MASTER, true);
            p_node_8->Set(SLAVE, false);
            p_node_8->Set(MASTER, true);
            p_cond_2->Set(SLAVE, false);
            p_cond_2->Set(MASTER, true);
            p_cond_3->Set(SLAVE, false);
            p_cond_3->Set(MASTER, true);
            
            // We compute the normals
            MortarUtilities::ComputeNodesMeanNormalModelPart(this_model_part);
            
            p_node_5->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_5->Z(), 2) + std::pow(p_node_5->Y(), 2);
            p_node_6->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_6->Z(), 2) + std::pow(p_node_6->Y(), 2);
            p_node_7->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_7->Z(), 2) + std::pow(p_node_7->Y(), 2);
            p_node_8->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_8->Z(), 2) + std::pow(p_node_8->Y(), 2);
                         
            typedef SimpleMortarMapperProcess<3, 3, Variable<double>, Historical> MapperType;
            MapperType process = MapperType(this_model_part, TEMPERATURE);
            process.Execute();
            
            // DEBUG         
//             GiDIODebugMapper(this_model_part);
            
            const double tolerance = 1.0e-3;
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_1->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_1->Z(), 2) + std::pow(p_node_1->Y(), 2))), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_2->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_2->Z(), 2) + std::pow(p_node_2->Y(), 2))), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_3->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_3->Z(), 2) + std::pow(p_node_3->Y(), 2))), tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(p_node_4->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_4->Z(), 2) + std::pow(p_node_4->Y(), 2))), tolerance);
        }
    } // namespace Testing
}  // namespace Kratos.
