// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

/* Processes */
#include "custom_processes/assign_nodal_elements_to_nodes_process.h"
#include "structural_mechanics_application_variables.h"
#include "structural_mechanics_fast_suite.h"

namespace Kratos::Testing
{
void AssignNodalElementsToNodesProcessCreateModelPart(ModelPart& ThisModelPart)
{
    Properties::Pointer p_elem_prop = ThisModelPart.pGetProperties(0);

    // First we create the nodes
    auto p_node_1 = ThisModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node_2 = ThisModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node_3 = ThisModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    auto p_node_4 = ThisModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    auto p_node_5 = ThisModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    auto p_node_6 = ThisModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);
}

/**
* Checks the correct work of the AssignNodalElementsToNodesProcess process
* Test 1. No rotation
*/
KRATOS_TEST_CASE_IN_SUITE(AssignNodalElementsToNodesProcess1, KratosStructuralMechanicsFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main");

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    AssignNodalElementsToNodesProcessCreateModelPart(r_model_part);

    Parameters parameters = Parameters(R"(
    {
        "nodal_mass"                     : 1.0,
        "nodal_stiffness"                : [1.0, null, 0.0]
    })" );

    AssignNodalElementsToNodesProcess assign_nodal_elements_to_nodes_process(r_model_part, parameters);
    assign_nodal_elements_to_nodes_process.Execute();

    KRATOS_EXPECT_EQ(r_model_part.Elements().size(), r_model_part.Nodes().size());

    for (auto& r_elem : r_model_part.Elements()) {
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_MASS], 1.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_DISPLACEMENT_STIFFNESS][0], 1.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_DISPLACEMENT_STIFFNESS][1], 0.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_DISPLACEMENT_STIFFNESS][2], 0.0);
    }
}

/**
* Checks the correct work of the AssignNodalElementsToNodesProcess process
* Test 2. WIth rotation
*/
KRATOS_TEST_CASE_IN_SUITE(AssignNodalElementsToNodesProcess2, KratosStructuralMechanicsFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 2);

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(ROTATION);
    r_model_part.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION);
    AssignNodalElementsToNodesProcessCreateModelPart(r_model_part);

    Parameters parameters = Parameters(R"(
    {
        "nodal_mass"                     : 1.0,
        "nodal_inertia"                  : [1.0, 1.0, 1.0],
        "nodal_stiffness"                : [1.0, null, 0.0],
        "nodal_rotational_stiffness"     : [1.0, null, 1.0]
    })" );

    AssignNodalElementsToNodesProcess assign_nodal_elements_to_nodes_process(r_model_part, parameters);
    assign_nodal_elements_to_nodes_process.Execute();

    KRATOS_EXPECT_EQ(r_model_part.Elements().size(), r_model_part.Nodes().size());

    for (auto& r_elem : r_model_part.Elements()) {
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_MASS], 1.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_DISPLACEMENT_STIFFNESS][0], 1.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_DISPLACEMENT_STIFFNESS][1], 0.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_DISPLACEMENT_STIFFNESS][2], 0.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_INERTIA][0], 1.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_INERTIA][1], 1.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_INERTIA][2], 1.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_ROTATIONAL_STIFFNESS][0], 1.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_ROTATIONAL_STIFFNESS][1], 0.0);
        KRATOS_EXPECT_EQ(r_elem.GetProperties()[NODAL_ROTATIONAL_STIFFNESS][2], 1.0);
    }
}

}  // namespace Kratos::Testing.
