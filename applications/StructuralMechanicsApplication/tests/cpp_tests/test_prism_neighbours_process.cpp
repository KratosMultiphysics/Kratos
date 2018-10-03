// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/prism_3d_6.h"
#include "testing/testing.h"
#include "includes/gid_io.h"
#include "includes/kernel.h"

/* Processes */
#include "custom_processes/prism_neighbours_process.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

        void PrismNeighboursProcessGiDIODebug(ModelPart& ThisModelPart)
        {
            GidIO<> gid_io("TEST_NEIGHBOUR_PRISM", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
            const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);

            gid_io.InitializeMesh(label);
            gid_io.WriteMesh(ThisModelPart.GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(label, ThisModelPart.GetMesh());
        }

        void PrismNeighboursProcessCreateModelPart(
            ModelPart& ThisModelPart,
            const std::size_t NumberNeighbours
            )
        {
            Properties::Pointer p_elem_prop = ThisModelPart.pGetProperties(0);

            // First we create the nodes
            NodeType::Pointer p_node_1  = ThisModelPart.CreateNewNode( 1,  0.0 ,  0.0 , 0.0);
            NodeType::Pointer p_node_2  = ThisModelPart.CreateNewNode( 2,  1.0 ,  0.0 , 0.0);
            NodeType::Pointer p_node_3  = ThisModelPart.CreateNewNode( 3,  0.0 ,  1.0 , 0.0);
            NodeType::Pointer p_node_4  = ThisModelPart.CreateNewNode( 4,  0.0 ,  0.0 , 1.0);
            NodeType::Pointer p_node_5  = ThisModelPart.CreateNewNode( 5,  1.0 ,  0.0 , 1.0);
            NodeType::Pointer p_node_6  = ThisModelPart.CreateNewNode( 6,  0.0 ,  1.0 , 1.0);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (6);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            element_nodes_0[3] = p_node_4;
            element_nodes_0[4] = p_node_5;
            element_nodes_0[5] = p_node_6;
            Prism3D6<NodeType> prism_0( PointerVector<NodeType>{element_nodes_0} );

            Element::Pointer p_elem_0 = ThisModelPart.CreateNewElement("Element3D6N", 1, prism_0, p_elem_prop);

            if (NumberNeighbours > 0) {
                NodeType::Pointer p_node_7  = ThisModelPart.CreateNewNode( 7,  1.0 ,  1.0 , 0.0);
                NodeType::Pointer p_node_10 = ThisModelPart.CreateNewNode(10,  1.0 ,  1.0 , 1.0);

                std::vector<NodeType::Pointer> element_nodes_1 (6);
                element_nodes_1[0] = p_node_2;
                element_nodes_1[1] = p_node_3;
                element_nodes_1[2] = p_node_7;
                element_nodes_1[3] = p_node_5;
                element_nodes_1[4] = p_node_6;
                element_nodes_1[5] = p_node_10;
                Prism3D6 <NodeType> prism_1( PointerVector<NodeType>{element_nodes_1} );

                Element::Pointer p_elem_1 = ThisModelPart.CreateNewElement("Element3D6N", 2, prism_1, p_elem_prop);
            }

            if (NumberNeighbours > 1) {
                NodeType::Pointer p_node_8  = ThisModelPart.CreateNewNode( 8,  1.0 , -1.0 , 0.0);
                NodeType::Pointer p_node_11 = ThisModelPart.CreateNewNode(11,  1.0 , -1.0 , 1.0);

                std::vector<NodeType::Pointer> element_nodes_2 (6);
                element_nodes_2[0] = p_node_1;
                element_nodes_2[1] = p_node_2;
                element_nodes_2[2] = p_node_8;
                element_nodes_2[3] = p_node_4;
                element_nodes_2[4] = p_node_5;
                element_nodes_2[5] = p_node_11;
                Prism3D6 <NodeType> prism_2( PointerVector<NodeType>{element_nodes_2} );

                Element::Pointer p_elem_2 = ThisModelPart.CreateNewElement("Element3D6N", 3, prism_2, p_elem_prop);
            }

            if (NumberNeighbours > 2) {
                NodeType::Pointer p_node_9  = ThisModelPart.CreateNewNode( 9, -1.0 ,  1.0 , 0.0);
                NodeType::Pointer p_node_12 = ThisModelPart.CreateNewNode(12, -1.0 ,  1.0 , 1.0);

                std::vector<NodeType::Pointer> element_nodes_3 (6);
                element_nodes_3[0] = p_node_1;
                element_nodes_3[1] = p_node_3;
                element_nodes_3[2] = p_node_9;
                element_nodes_3[3] = p_node_4;
                element_nodes_3[4] = p_node_6;
                element_nodes_3[5] = p_node_12;
                Prism3D6 <NodeType> prism_3( PointerVector<NodeType>{element_nodes_3} );

                Element::Pointer p_elem_3 = ThisModelPart.CreateNewElement("Element3D6N", 4, prism_3, p_elem_prop);
            }

//             // DEBUG
//             PrismNeighboursProcessGiDIODebug(ThisModelPart);
        }

        /**
        * Checks the correct work of the prism neighbour process
        * Test 3 neighbours
        */

        KRATOS_TEST_CASE_IN_SUITE(PrismNeighboursProcess1, KratosStructuralMechanicsFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            PrismNeighboursProcessCreateModelPart(this_model_part, 3);

            PrismNeighboursProcess prism_neighbours_process = PrismNeighboursProcess(this_model_part);
            prism_neighbours_process.Execute();

            auto pneigh = (this_model_part.Elements().begin())->GetValue(NEIGHBOUR_NODES);
            KRATOS_CHECK_EQUAL(pneigh[0].Id(), 7);
            KRATOS_CHECK_EQUAL(pneigh[1].Id(), 9);
            KRATOS_CHECK_EQUAL(pneigh[2].Id(), 8);
        }

        /**
        * Checks the correct work of the prism neighbour process
        * Test 2 neighbours
        */
        KRATOS_TEST_CASE_IN_SUITE(PrismNeighboursProcess2, KratosStructuralMechanicsFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            PrismNeighboursProcessCreateModelPart(this_model_part, 2);

            PrismNeighboursProcess prism_neighbours_process = PrismNeighboursProcess(this_model_part);
            prism_neighbours_process.Execute();

            auto pneigh = (this_model_part.Elements().begin())->GetValue(NEIGHBOUR_NODES);
            KRATOS_CHECK_EQUAL(pneigh[0].Id(), 7);
            KRATOS_CHECK_EQUAL(pneigh[1].Id(), 2);
            KRATOS_CHECK_EQUAL(pneigh[2].Id(), 8);
        }

        /**
        * Checks the correct work of the prism neighbour process
        * Test 1 neighbours
        */
        KRATOS_TEST_CASE_IN_SUITE(PrismNeighboursProcess3, KratosStructuralMechanicsFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            PrismNeighboursProcessCreateModelPart(this_model_part, 1);

            PrismNeighboursProcess prism_neighbours_process = PrismNeighboursProcess(this_model_part);
            prism_neighbours_process.Execute();

            auto pneigh = (this_model_part.Elements().begin())->GetValue(NEIGHBOUR_NODES);
            KRATOS_CHECK_EQUAL(pneigh[0].Id(), 7);
            KRATOS_CHECK_EQUAL(pneigh[1].Id(), 2);
            KRATOS_CHECK_EQUAL(pneigh[2].Id(), 3);
        }

        /**
        * Checks the correct work of the prism neighbour process
        * Test 0 neighbours (It will return the ID of the original element)
        */
        KRATOS_TEST_CASE_IN_SUITE(PrismNeighboursProcess4, KratosStructuralMechanicsFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            PrismNeighboursProcessCreateModelPart(this_model_part, 0);

            PrismNeighboursProcess prism_neighbours_process = PrismNeighboursProcess(this_model_part);
            prism_neighbours_process.Execute();

            auto pneigh = (this_model_part.Elements().begin())->GetValue(NEIGHBOUR_NODES);
            KRATOS_CHECK_EQUAL(pneigh[0].Id(), 1);
            KRATOS_CHECK_EQUAL(pneigh[1].Id(), 2);
            KRATOS_CHECK_EQUAL(pneigh[2].Id(), 3);
        }

    } // namespace Testing
}  // namespace Kratos.
