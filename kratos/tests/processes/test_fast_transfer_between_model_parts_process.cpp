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
#include "geometries/triangle_3d_3.h"

/* Processes */
#include "processes/fast_transfer_between_model_parts_process.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3> NodeType;
        
        /**
        * Checks the correct work of the fast_transfer_between_model_parts_process
        * Test 1
        */

        KRATOS_TEST_CASE_IN_SUITE(FastTransferBetweenModelPartsProcess1, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& origin_model_part = current_model.CreateModelPart("Origin");
            ModelPart& destination_model_part = current_model.CreateModelPart("Destination");

            Properties::Pointer p_cond_prop = origin_model_part.pGetProperties(0);

            auto& process_info = origin_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            // First we create the nodes
            NodeType::Pointer p_node_1 = origin_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_2 = origin_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_3 = origin_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);

            NodeType::Pointer p_node_4 = origin_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_5 = origin_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_6 = origin_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_3;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_1;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            Condition::Pointer p_cond_0 = origin_model_part.CreateNewCondition("Condition3D", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond_1 = origin_model_part.CreateNewCondition("Condition3D", 2, triangle_1, p_cond_prop);

            // This will copy all
            FastTransferBetweenModelPartsProcess process = FastTransferBetweenModelPartsProcess(destination_model_part, origin_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL);
            process.Execute();

            std::size_t count = 0;
            for (auto& node : origin_model_part.Nodes()) {
                ++count;
                KRATOS_CHECK_EQUAL(node.Id(), destination_model_part.GetNode(count).Id());
            }

            count = 0;
            for (auto& cond : origin_model_part.Conditions()) {
                ++count;
                KRATOS_CHECK_EQUAL(cond.Id(), destination_model_part.GetCondition(count).Id());
            }
        }

        /**
        * Checks the correct work of the fast_transfer_between_model_parts_process
        * Test 2 (with flags)
        */

        KRATOS_TEST_CASE_IN_SUITE(FastTransferBetweenModelPartsProcess2, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& origin_model_part = current_model.CreateModelPart("Origin");
            ModelPart& destination_model_part = current_model.CreateModelPart("Destination");

            Properties::Pointer p_cond_prop = origin_model_part.pGetProperties(0);

            auto& process_info = origin_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            // First we create the nodes
            NodeType::Pointer p_node_1 = origin_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_2 = origin_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_3 = origin_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);

            NodeType::Pointer p_node_4 = origin_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_5 = origin_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_6 = origin_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_3;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_1;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            Condition::Pointer p_cond_0 = origin_model_part.CreateNewCondition("Condition3D", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond_1 = origin_model_part.CreateNewCondition("Condition3D", 2, triangle_1, p_cond_prop);

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

            // This will copy all
            FastTransferBetweenModelPartsProcess process = FastTransferBetweenModelPartsProcess(destination_model_part, origin_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, MASTER);
            process.Execute();

            std::size_t count = 0;
            for (auto& node : origin_model_part.Nodes()) {
                ++count;
                if (node.Is(MASTER))
                    KRATOS_CHECK_EQUAL(node.Id(), destination_model_part.GetNode(count).Id());
            }

            count = 0;
            for (auto& cond : origin_model_part.Conditions()) {
                ++count;
                if (cond.Is(MASTER))
                    KRATOS_CHECK_EQUAL(cond.Id(), destination_model_part.GetCondition(count).Id());
            }
        }

        /**
        * Checks the correct work of the fast_transfer_between_model_parts_process
        * Test 3 (clone/replicate)
        */

        KRATOS_TEST_CASE_IN_SUITE(FastTransferBetweenModelPartsProcess3, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& origin_model_part = current_model.CreateModelPart("Origin");
            ModelPart& destination_model_part = current_model.CreateModelPart("Destination");

            Properties::Pointer p_cond_prop = origin_model_part.pGetProperties(0);

            auto& process_info = origin_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            // First we create the nodes
            NodeType::Pointer p_node_1 = origin_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_2 = origin_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_3 = origin_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);

            NodeType::Pointer p_node_4 = origin_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_5 = origin_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_6 = origin_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_3;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_1;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            Condition::Pointer p_cond_0 = origin_model_part.CreateNewCondition("Condition3D", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond_1 = origin_model_part.CreateNewCondition("Condition3D", 2, triangle_1, p_cond_prop);

            // This will copy all
            FastTransferBetweenModelPartsProcess process = FastTransferBetweenModelPartsProcess(destination_model_part, origin_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, Flags(), true);
            process.Execute();

            std::size_t count = 0;
            for (auto& node : destination_model_part.Nodes()) {
                ++count;
                KRATOS_CHECK_EQUAL(node.Id(), origin_model_part.GetNode(count).Id() + 6);
            }

            count = 0;
            for (auto& cond : destination_model_part.Conditions()) {
                ++count;
                KRATOS_CHECK_EQUAL(cond.Id(), origin_model_part.GetCondition(count).Id() + 2);
            }
        }

        /**
        * Checks the correct work of the fast_transfer_between_model_parts_process
        * Test 4 (clone/replicate with flags)
        */

        KRATOS_TEST_CASE_IN_SUITE(FastTransferBetweenModelPartsProcess4, KratosCoreFastSuite)
        {
            Model current_model;
            
            ModelPart& origin_model_part = current_model.CreateModelPart("Origin");
            ModelPart& destination_model_part = current_model.CreateModelPart("Destination");
            
            Properties::Pointer p_cond_prop = origin_model_part.pGetProperties(0);
            
            auto& process_info = origin_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = origin_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_2 = origin_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
            NodeType::Pointer p_node_3 = origin_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);
            
            NodeType::Pointer p_node_4 = origin_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_5 = origin_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
            NodeType::Pointer p_node_6 = origin_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_3;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_1;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );
            
            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );
            
            Condition::Pointer p_cond_0 = origin_model_part.CreateNewCondition("Condition3D", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond_1 = origin_model_part.CreateNewCondition("Condition3D", 2, triangle_1, p_cond_prop);
            
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
                         
            // This will copy all
            FastTransferBetweenModelPartsProcess process = FastTransferBetweenModelPartsProcess(destination_model_part, origin_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, MASTER, true);
            process.Execute();

            std::size_t count = 0;
            for (auto& node : origin_model_part.Nodes()) {
                ++count;
                if (node.Is(MASTER))
                    KRATOS_CHECK_EQUAL(node.Id() + 6, destination_model_part.GetNode(count + 6).Id());
            }

            count = 0;
            for (auto& cond : origin_model_part.Conditions()) {
                ++count;
                if (cond.Is(MASTER))
                    KRATOS_CHECK_EQUAL(cond.Id() + 2, destination_model_part.GetCondition(count  + 2).Id());
            }
        }

    } // namespace Testing
}  // namespace Kratos.
