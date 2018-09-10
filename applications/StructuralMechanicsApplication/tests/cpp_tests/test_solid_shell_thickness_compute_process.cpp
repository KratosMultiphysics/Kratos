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
#include "includes/checks.h"
#include "geometries/prism_3d_6.h"

/* Processes */
#include "custom_processes/solid_shell_thickness_compute_process.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

        void SolidShellProcessCreateModelPart(ModelPart& ThisModelPart)
        {
            Properties::Pointer p_elem_prop = ThisModelPart.pGetProperties(0);

            // First we create the nodes
            NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = ThisModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = ThisModelPart.CreateNewNode(4, 0.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_5 = ThisModelPart.CreateNewNode(5, 1.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_6 = ThisModelPart.CreateNewNode(6, 1.0 , 1.0 , 1.0);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (6);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            element_nodes_0[3] = p_node_4;
            element_nodes_0[4] = p_node_5;
            element_nodes_0[5] = p_node_6;
            Prism3D6 <NodeType> prism_0( PointerVector<NodeType>{element_nodes_0} );
            ThisModelPart.CreateNewElement("Element3D6N", 1, prism_0, p_elem_prop);
        }

        /**
        * Checks the correct work of the thickness compute for solid shells
        */

        KRATOS_TEST_CASE_IN_SUITE(TestSolidShellThicknessCompute, KratosStructuralMechanicsFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);

            SolidShellProcessCreateModelPart(this_model_part);

            SolidShellThickComputeProcess thickness_process = SolidShellThickComputeProcess(this_model_part);
            thickness_process.Execute();

            for (auto& node : this_model_part.Nodes())
                KRATOS_CHECK_NEAR(node.GetValue(THICKNESS), 1.0, std::numeric_limits<double>::epsilon());
        }

    } // namespace Testing
}  // namespace Kratos.
