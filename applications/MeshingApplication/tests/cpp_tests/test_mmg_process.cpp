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
#ifdef INCLUDE_MMG
// System includes

// External includes

// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"
#include "includes/gid_io.h"
#include "meshing_application.h"

/* Processes */
#include "custom_processes/mmg_process.h"
#include "processes/find_nodal_h_process.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3>                                                    NodeType;
        
        void GiDIODebugMMG(ModelPart& ThisModelPart)
        {
            GidIO<> gid_io("TEST_MMG", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
            const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);

            gid_io.InitializeMesh(label);
            gid_io.WriteMesh(ThisModelPart.GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(label, ThisModelPart.GetMesh());
            gid_io.WriteNodalResults(NODAL_H, ThisModelPart.Nodes(), label, 0);
        }
        
        /** 
        * Checks the correct work of the level set MMG process
        * Test triangle 
        */

        KRATOS_TEST_CASE_IN_SUITE(TestMMGProcess1, KratosMMGProcessFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            
            this_model_part.AddNodalSolutionStepVariable(NODAL_H);
            
            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_0( element_nodes_0 );
            
            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;
            Triangle2D3 <NodeType> triangle_1( element_nodes_1 );
            
            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_2( element_nodes_2 );
            
            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_3( element_nodes_3 );
            
            Element::Pointer p_elem_0 = this_model_part.CreateNewElement("Element2D3N", 1, triangle_0, p_elem_prop);
            Element::Pointer p_elem_1 = this_model_part.CreateNewElement("Element2D3N", 2, triangle_1, p_elem_prop);
            Element::Pointer p_elem_2 = this_model_part.CreateNewElement("Element2D3N", 3, triangle_2, p_elem_prop);
            Element::Pointer p_elem_3 = this_model_part.CreateNewElement("Element2D3N", 4, triangle_3, p_elem_prop);
            
            // Set DISTANCE and other variables
            Vector ref_metric(3);
            ref_metric[0] = 1.0;
            ref_metric[1] = 0;
            ref_metric[2] = 1.0;
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->SetValue(MMG_METRIC, ref_metric);
            }
                         
            // Compute remesh
            Parameters params = Parameters(R"({ "echo_level" : 0 })" );
            MmgProcess<2> mmg_process = MmgProcess<2>(this_model_part, params);
            mmg_process.Execute();
            
            // Compute NodalH
            FindNodalHProcess process = FindNodalHProcess(this_model_part);
            process.Execute();
            
//             // DEBUG         
//             GiDIODebugMMG(this_model_part);
            
            const double tolerance = 1.0e-4;
            for (auto& i_node : this_model_part.Nodes())
                KRATOS_CHECK_LESS_EQUAL(std::abs(i_node.FastGetSolutionStepValue(NODAL_H) - 1.0/std::sqrt(2.0))*std::sqrt(2.0), tolerance);
        }
        
        /** 
        * Checks the correct work of the level set MMG process
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(TestMMGProcess2, KratosMMGProcessFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            
            this_model_part.AddNodalSolutionStepVariable(NODAL_H);
            this_model_part.AddNodalSolutionStepVariable(DISTANCE);
            this_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);
            this_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
            
            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6 , 1.0 , 1.0 , 0.0);
            
            NodeType::Pointer p_node_7 = this_model_part.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_8 = this_model_part.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_9 = this_model_part.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_10 = this_model_part.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_11 = this_model_part.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_12 = this_model_part.CreateNewNode(12 , 2.0 , 0.0 , 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (4);
            element_nodes_0[0] = p_node_12;
            element_nodes_0[1] = p_node_10;
            element_nodes_0[2] = p_node_8;
            element_nodes_0[3] = p_node_9;
            Tetrahedra3D4 <NodeType> tetrahedra_0( element_nodes_0 );
            
            std::vector<NodeType::Pointer> element_nodes_1 (4);
            element_nodes_1[0] = p_node_4;
            element_nodes_1[1] = p_node_6;
            element_nodes_1[2] = p_node_9;
            element_nodes_1[3] = p_node_7;
            Tetrahedra3D4 <NodeType> tetrahedra_1( element_nodes_1 );
            
            std::vector<NodeType::Pointer> element_nodes_2 (4);
            element_nodes_2[0] = p_node_11;
            element_nodes_2[1] = p_node_7;
            element_nodes_2[2] = p_node_9;
            element_nodes_2[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_2( element_nodes_2 );
            
            std::vector<NodeType::Pointer> element_nodes_3 (4);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_3;
            element_nodes_3[2] = p_node_8;
            element_nodes_3[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_3( element_nodes_3 );
            
            std::vector<NodeType::Pointer> element_nodes_4 (4);
            element_nodes_4[0] = p_node_4;
            element_nodes_4[1] = p_node_6;
            element_nodes_4[2] = p_node_7;
            element_nodes_4[3] = p_node_3;
            Tetrahedra3D4 <NodeType> tetrahedra_4( element_nodes_4 );
            
            std::vector<NodeType::Pointer> element_nodes_5 (4);
            element_nodes_5[0] = p_node_2;
            element_nodes_5[1] = p_node_3;
            element_nodes_5[2] = p_node_5;
            element_nodes_5[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_5( element_nodes_5 );
            
            std::vector<NodeType::Pointer> element_nodes_6 (4);
            element_nodes_6[0] = p_node_10;
            element_nodes_6[1] = p_node_9;
            element_nodes_6[2] = p_node_6;
            element_nodes_6[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_6( element_nodes_6 );
            
            std::vector<NodeType::Pointer> element_nodes_7 (4);
            element_nodes_7[0] = p_node_7;
            element_nodes_7[1] = p_node_8;
            element_nodes_7[2] = p_node_3;
            element_nodes_7[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_7( element_nodes_7 );
            
            std::vector<NodeType::Pointer> element_nodes_8 (4);
            element_nodes_8[0] = p_node_7;
            element_nodes_8[1] = p_node_8;
            element_nodes_8[2] = p_node_6;
            element_nodes_8[3] = p_node_9;
            Tetrahedra3D4 <NodeType> tetrahedra_8( element_nodes_8 );
            
            std::vector<NodeType::Pointer> element_nodes_9 (4);
            element_nodes_9[0] = p_node_4;
            element_nodes_9[1] = p_node_1;
            element_nodes_9[2] = p_node_6;
            element_nodes_9[3] = p_node_3;
            Tetrahedra3D4 <NodeType> tetrahedra_9( element_nodes_9 );
            
            std::vector<NodeType::Pointer> element_nodes_10 (4);
            element_nodes_10[0] = p_node_9;
            element_nodes_10[1] = p_node_12;
            element_nodes_10[2] = p_node_11;
            element_nodes_10[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_10( element_nodes_10 );
            
            std::vector<NodeType::Pointer> element_nodes_11 (4);
            element_nodes_11[0] = p_node_3;
            element_nodes_11[1] = p_node_2;
            element_nodes_11[2] = p_node_1;
            element_nodes_11[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_11( element_nodes_11 );
            
            Element::Pointer p_elem_0 = this_model_part.CreateNewElement("Element3D4N", 1, tetrahedra_0, p_elem_prop);
            Element::Pointer p_elem_1 = this_model_part.CreateNewElement("Element3D4N", 2, tetrahedra_1, p_elem_prop);
            Element::Pointer p_elem_2 = this_model_part.CreateNewElement("Element3D4N", 3, tetrahedra_2, p_elem_prop);
            Element::Pointer p_elem_3 = this_model_part.CreateNewElement("Element3D4N", 4, tetrahedra_3, p_elem_prop);
            Element::Pointer p_elem_4 = this_model_part.CreateNewElement("Element3D4N", 5, tetrahedra_4, p_elem_prop);
            Element::Pointer p_elem_5 = this_model_part.CreateNewElement("Element3D4N", 6, tetrahedra_5, p_elem_prop);
            Element::Pointer p_elem_6 = this_model_part.CreateNewElement("Element3D4N", 7, tetrahedra_6, p_elem_prop);
            Element::Pointer p_elem_7 = this_model_part.CreateNewElement("Element3D4N", 8, tetrahedra_7, p_elem_prop);
            Element::Pointer p_elem_8 = this_model_part.CreateNewElement("Element3D4N", 9, tetrahedra_8, p_elem_prop);
            Element::Pointer p_elem_9 = this_model_part.CreateNewElement("Element3D4N", 10, tetrahedra_9, p_elem_prop);
            Element::Pointer p_elem_10 = this_model_part.CreateNewElement("Element3D4N", 11, tetrahedra_10, p_elem_prop);
            Element::Pointer p_elem_11 = this_model_part.CreateNewElement("Element3D4N", 12, tetrahedra_11, p_elem_prop);
            
            // Set DISTANCE and other variables
            Vector ref_metric = ZeroVector(6);
            ref_metric[0] = 1.0;
            ref_metric[3] = 1.0;
            ref_metric[5] = 1.0;
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->SetValue(MMG_METRIC, ref_metric);
            }
                      
            // Compute remesh
            Parameters params = Parameters(R"({ "echo_level" : 0 })" );
            MmgProcess<3> mmg_process = MmgProcess<3>(this_model_part, params);
            mmg_process.Execute();
            
            // Compute NodalH
            FindNodalHProcess process = FindNodalHProcess(this_model_part);
            process.Execute();
            
//             // DEBUG         
//             GiDIODebugMMG(this_model_part);
            
            double max = 0.0;
            for (auto& i_node : this_model_part.Nodes())
                if (i_node.FastGetSolutionStepValue(NODAL_H) > max) 
                    max = i_node.FastGetSolutionStepValue(NODAL_H);
                
            const double tolerance = 1.0e-2;
            KRATOS_CHECK_LESS_EQUAL(std::abs(max - 1.0/std::sqrt(2.0))/max, tolerance);
        }
    } // namespace Testing
}  // namespace Kratos.
#endif
