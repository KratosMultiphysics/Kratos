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

/* Processes */
#include "custom_processes/assign_nodal_elements_to_nodes.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3> NodeType;
        
        void AssignNodalElementsToNodesCreateModelPart(ModelPart& ThisModelPart)
        {
            Properties::Pointer p_elem_prop = ThisModelPart.pGetProperties(0);

            // First we create the nodes
            NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = ThisModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = ThisModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = ThisModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = ThisModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);
        }

        /**
        * Checks the correct work of the AssignNodalElementsToNodes process
        * Test 1. No rotation
        */

        KRATOS_TEST_CASE_IN_SUITE(TestAssignNodalElementsToNodes1, KratosStructuralMechanicsFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            this_model_part.AddNodalSolutionStepVariable(VELOCITY);
            this_model_part.AddNodalSolutionStepVariable(ACCELERATION);
            AssignNodalElementsToNodesCreateModelPart(this_model_part);

            Parameters parameters = Parameters(R"(
            {
                "nodal_mass"                     : 1.0,
                "nodal_stiffness"                : [1.0, null, 0.0]
            })" );

            AssignNodalElementsToNodes assign_nodal_elements_to_nodes = AssignNodalElementsToNodes(this_model_part, parameters);
            assign_nodal_elements_to_nodes.Execute();

            KRATOS_CHECK_EQUAL(this_model_part.Elements().size(), this_model_part.Nodes().size());

            for (auto& elem : this_model_part.Elements()) {
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_MASS], 1.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_STIFFNESS][0], 1.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_STIFFNESS][1], 0.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_STIFFNESS][2], 0.0);
            }
        }

        /**
        * Checks the correct work of the AssignNodalElementsToNodes process
        * Test 2. WIth rotation
        */

        KRATOS_TEST_CASE_IN_SUITE(TestAssignNodalElementsToNodes2, KratosStructuralMechanicsFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            this_model_part.AddNodalSolutionStepVariable(VELOCITY);
            this_model_part.AddNodalSolutionStepVariable(ACCELERATION);
            this_model_part.AddNodalSolutionStepVariable(ROTATION);
            this_model_part.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION);
            this_model_part.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION);
            AssignNodalElementsToNodesCreateModelPart(this_model_part);

            Parameters parameters = Parameters(R"(
            {
                "nodal_mass"                     : 1.0,
                "nodal_inertia"                  : [1.0, 1.0, 1.0],
                "nodal_stiffness"                : [1.0, null, 0.0],
                "nodal_rotational_stiffness"     : [1.0, null, 1.0]
            })" );

            AssignNodalElementsToNodes assign_nodal_elements_to_nodes = AssignNodalElementsToNodes(this_model_part, parameters);
            assign_nodal_elements_to_nodes.Execute();

            KRATOS_CHECK_EQUAL(this_model_part.Elements().size(), this_model_part.Nodes().size());

            for (auto& elem : this_model_part.Elements()) {
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_MASS], 1.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_STIFFNESS][0], 1.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_STIFFNESS][1], 0.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_STIFFNESS][2], 0.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_INERTIA][0], 1.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_INERTIA][1], 1.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_INERTIA][2], 1.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_ROTATIONAL_STIFFNESS][0], 1.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_ROTATIONAL_STIFFNESS][1], 0.0);
                KRATOS_CHECK_EQUAL(elem.GetProperties()[NODAL_ROTATIONAL_STIFFNESS][2], 1.0);
            }
        }

    } // namespace Testing
}  // namespace Kratos.
