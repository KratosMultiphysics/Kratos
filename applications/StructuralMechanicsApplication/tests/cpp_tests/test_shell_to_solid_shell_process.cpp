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
#include "testing/testing.h"
#include "includes/gid_io.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/mortar_utilities.h"

/* Processes */
#include "custom_processes/shell_to_solid_shell_process.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

        void ShellToSolidShellProcessGiDIODebug(ModelPart& ThisModelPart)
        {
            GidIO<> gid_io("TEST_SHELL_TO_SOLID", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditions);
            const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);

            gid_io.InitializeMesh(label);
            gid_io.WriteMesh(ThisModelPart.GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(label, ThisModelPart.GetMesh());
            gid_io.WriteNodalResults(NORMAL, ThisModelPart.Nodes(), label, 0);
            gid_io.WriteNodalResultsNonHistorical(NORMAL, ThisModelPart.Nodes(), label);
            gid_io.WriteNodalResultsNonHistorical(THICKNESS, ThisModelPart.Nodes(), label);
            gid_io.WriteNodalResultsNonHistorical(NODAL_AREA, ThisModelPart.Nodes(), label);
        }

        void ShellToSolidShellProcessCreateModelPart(ModelPart& ThisModelPart)
        {
            Properties::Pointer p_elem_prop = ThisModelPart.pGetProperties(0);

            p_elem_prop->SetValue(THICKNESS, 0.2);

            // First we create the nodes
            NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = ThisModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = ThisModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = ThisModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = ThisModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{element_nodes_0} );
            ThisModelPart.CreateNewElement("Element3D3N", 1, triangle_0, p_elem_prop);

            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{element_nodes_1} );
            ThisModelPart.CreateNewElement("Element3D3N", 2, triangle_1, p_elem_prop);

            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_2( PointerVector<NodeType>{element_nodes_2} );
            ThisModelPart.CreateNewElement("Element3D3N", 3, triangle_2, p_elem_prop);

            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_2;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_3( PointerVector<NodeType>{element_nodes_3} );
            ThisModelPart.CreateNewElement("Element3D3N", 4, triangle_3, p_elem_prop);
        }

        /**
        * Checks the correct work of the shell to solid process
        * Test 1 layer
        */

        KRATOS_TEST_CASE_IN_SUITE(TestShellToSolidShellProcess1, KratosStructuralMechanicsFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part =  current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(NORMAL);
            ShellToSolidShellProcessCreateModelPart(this_model_part);

            Parameters parameters = Parameters(R"(
            {
                "element_name"    : "SolidShellElementSprism3D6N",
                "model_part_name" : "",
                "computing_model_part_name" : "",
                "number_of_layers": 1
            })" );

            ShellToSolidShellProcess<3> prism_neighbours_process = ShellToSolidShellProcess<3>(this_model_part, parameters);
            prism_neighbours_process.Execute();

//             // DEBUG
//             ShellToSolidShellProcessGiDIODebug(this_model_part);

            for (auto& elem : this_model_part.Elements())
                KRATOS_CHECK_EQUAL(elem.GetGeometry().size(), 6);
        }

        /**
        * Checks the correct work of the shell to solid process
        * Test 2 layer
        */

        KRATOS_TEST_CASE_IN_SUITE(TestShellToSolidShellProcess2, KratosStructuralMechanicsFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part =  current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(NORMAL);
            ShellToSolidShellProcessCreateModelPart(this_model_part);

            Parameters parameters = Parameters(R"(
            {
                "element_name"    : "SolidShellElementSprism3D6N",
                "model_part_name" : "",
                "computing_model_part_name" : "",
                "number_of_layers": 2
            })" );

            ShellToSolidShellProcess<3> prism_neighbours_process = ShellToSolidShellProcess<3>(this_model_part, parameters);
            prism_neighbours_process.Execute();

//             // DEBUG
//             ShellToSolidShellProcessGiDIODebug(this_model_part);

            for (auto& elem : this_model_part.Elements())
                KRATOS_CHECK_EQUAL(elem.GetGeometry().size(), 6);
        }

        /**
        * Checks the correct work of the shell to solid process
        * Test 2 layer with external conditions
        */

        KRATOS_TEST_CASE_IN_SUITE(TestShellToSolidShellProcess3, KratosStructuralMechanicsFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(NORMAL);
            ShellToSolidShellProcessCreateModelPart(this_model_part);

            Parameters parameters = Parameters(R"(
            {
                "element_name"    : "SolidShellElementSprism3D6N",
                "model_part_name" : "",
                "computing_model_part_name" : "",
                "create_submodelparts_external_layers": true,
                "number_of_layers": 2
            })" );

            ShellToSolidShellProcess<3> prism_neighbours_process = ShellToSolidShellProcess<3>(this_model_part, parameters);
            prism_neighbours_process.Execute();

            // We compute the normal
            MortarUtilities::ComputeNodesMeanNormalModelPart(this_model_part.GetSubModelPart("Upper_"));
            MortarUtilities::ComputeNodesMeanNormalModelPart(this_model_part.GetSubModelPart("Lower_"));

//             // DEBUG
//             ShellToSolidShellProcessGiDIODebug(this_model_part);

            for (auto& elem : this_model_part.GetSubModelPart("Upper_").Conditions())
                KRATOS_CHECK_EQUAL(elem.GetGeometry().size(), 3);
            for (auto& elem : this_model_part.GetSubModelPart("Lower_").Conditions())
                KRATOS_CHECK_EQUAL(elem.GetGeometry().size(), 3);

            for (auto& node : this_model_part.GetSubModelPart("Upper_").Nodes())
                KRATOS_CHECK_NEAR(node.FastGetSolutionStepValue(NORMAL)[2],  1.0, 1.0e-12);
            for (auto& node : this_model_part.GetSubModelPart("Lower_").Nodes())
                KRATOS_CHECK_NEAR(node.FastGetSolutionStepValue(NORMAL)[2], -1.0, 1.0e-12);
        }

    } // namespace Testing
}  // namespace Kratos.
