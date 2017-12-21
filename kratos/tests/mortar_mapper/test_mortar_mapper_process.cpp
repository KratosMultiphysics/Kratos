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
        typedef Point                                                     PointType;
        typedef Node<3>                                                    NodeType;
        typedef Geometry<NodeType>                                 GeometryNodeType;
        typedef Geometry<PointType>                               GeometryPointType;

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
            
//             // DEBUG            
//             GidIO<> gid_io("TEST_MAPPER", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditionsOnly);
//             const int nl_iter = process_info[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
// 
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(this_model_part.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, this_model_part.GetMesh());
//             gid_io.WriteNodalResults(TEMPERATURE, this_model_part.Nodes(), label, 0);
            
            const double tolerance = 1.0e-6;
//             KRATOS_CHECK_NEAR(p_node_1->FastGetSolutionStepValue(TEMPERATURE), std::pow(p_node_1->X(), 2) + std::pow(p_node_1->Y(), 2), tolerance);
//             KRATOS_CHECK_NEAR(p_node_2->FastGetSolutionStepValue(TEMPERATURE), std::pow(p_node_2->X(), 2) + std::pow(p_node_2->Y(), 2), tolerance);
//             KRATOS_CHECK_NEAR(p_node_3->FastGetSolutionStepValue(TEMPERATURE), std::pow(p_node_3->X(), 2) + std::pow(p_node_3->Y(), 2), tolerance);
        }
        
    } // namespace Testing
}  // namespace Kratos.
