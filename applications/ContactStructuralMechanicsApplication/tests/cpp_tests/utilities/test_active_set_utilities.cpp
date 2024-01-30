// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "contact_structural_mechanics_application_variables.h"
#include "custom_utilities/active_set_utilities.h"

namespace Kratos::Testing
{
/**
* Checks the correct work of the acttive set utilities
* Test ComputePenaltyFrictionlessActiveSet
*/
KRATOS_TEST_CASE_IN_SUITE(ComputePenaltyFrictionlessActiveSet, KratosContactStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);
    ModelPart& r_contact_model_part = r_model_part.CreateSubModelPart("Contact");

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);

    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;
    r_process_info[DELTA_TIME] = 1.0;
    r_process_info[INITIAL_PENALTY] = 1.0e3;

    // First we create the nodes
    auto p_node1 = r_contact_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node2 = r_contact_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node3 = r_contact_model_part.CreateNewNode(3, 2.0 , 0.0 , 0.0);

    // Set flags
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Set(SLAVE);
        r_node.Set(ACTIVE);
    }

    // Set values
    p_node1->FastGetSolutionStepValue(WEIGHTED_GAP) = -1.0e-3;
    p_node2->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
    p_node3->FastGetSolutionStepValue(WEIGHTED_GAP) = 1.0e-3;

    // Call utilities
    const auto is_converged = ActiveSetUtilities::ComputePenaltyFrictionlessActiveSet(r_model_part);

    // Check flags
    KRATOS_EXPECT_EQ(is_converged, 2);
    KRATOS_EXPECT_TRUE(p_node1->Is(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node2->IsNot(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node3->IsNot(ACTIVE));
}

/**
* Checks the correct work of the acttive set utilities
* Test ComputePenaltyFrictionalActiveSet
*/
KRATOS_TEST_CASE_IN_SUITE(ComputePenaltyFrictionalActiveSet, KratosContactStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);
    ModelPart& r_contact_model_part = r_model_part.CreateSubModelPart("Contact");

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
    r_model_part.AddNodalSolutionStepVariable(WEIGHTED_SLIP);

    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;
    r_process_info[DELTA_TIME] = 1.0;
    r_process_info[INITIAL_PENALTY] = 1.0e3;
    r_process_info[TANGENT_FACTOR] = 1.0;

    // First we create the nodes
    auto p_node1 = r_contact_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node2 = r_contact_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node3 = r_contact_model_part.CreateNewNode(3, 2.0 , 0.0 , 0.0);

    // Set flags
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Set(SLAVE);
        r_node.Set(ACTIVE);
        r_node.Set(SLIP);
        r_node.SetValue(FRICTION_COEFFICIENT, 1.0);
    }

    // Set values
    p_node1->FastGetSolutionStepValue(WEIGHTED_GAP) = -1.0e-3;
    p_node1->FastGetSolutionStepValue(WEIGHTED_SLIP_X) = -1.0e-3;
    p_node2->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
    p_node2->FastGetSolutionStepValue(WEIGHTED_SLIP_X) = 0.0;
    p_node3->FastGetSolutionStepValue(WEIGHTED_GAP) = 1.0e-3;
    p_node3->FastGetSolutionStepValue(WEIGHTED_SLIP_X) = 1.0e-3;

    // Call utilities
    auto is_converged = ActiveSetUtilities::ComputePenaltyFrictionalActiveSet(r_model_part);

    // Check flags
    KRATOS_EXPECT_EQ(is_converged[0], 2);
    KRATOS_EXPECT_EQ(is_converged[1], 1);
    KRATOS_EXPECT_TRUE(p_node1->Is(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node1->IsNot(SLIP));
    KRATOS_EXPECT_TRUE(p_node2->IsNot(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node2->IsNot(SLIP));
    KRATOS_EXPECT_TRUE(p_node3->IsNot(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node3->IsNot(SLIP));

    // Set flags
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Set(ACTIVE);
        r_node.Set(SLIP);
    }

    // Call utilities (pure slip)
    is_converged = ActiveSetUtilities::ComputePenaltyFrictionalActiveSet(r_model_part, true);

    // Check flags
    KRATOS_EXPECT_EQ(is_converged[0], 2);
    KRATOS_EXPECT_EQ(is_converged[1], 0);
    KRATOS_EXPECT_TRUE(p_node1->Is(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node1->Is(SLIP));
    KRATOS_EXPECT_TRUE(p_node2->IsNot(ACTIVE));
    KRATOS_EXPECT_TRUE(!p_node2->IsDefined(SLIP));
    KRATOS_EXPECT_TRUE(p_node3->IsNot(ACTIVE));
    KRATOS_EXPECT_TRUE(!p_node3->IsDefined(SLIP));
}

/**
* Checks the correct work of the acttive set utilities
* Test ComputeALMFrictionlessActiveSet
*/
KRATOS_TEST_CASE_IN_SUITE(ComputeALMFrictionlessActiveSet, KratosContactStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);
    ModelPart& r_contact_model_part = r_model_part.CreateSubModelPart("Contact");

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
    r_model_part.AddNodalSolutionStepVariable(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);

    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;
    r_process_info[DELTA_TIME] = 1.0;
    r_process_info[INITIAL_PENALTY] = 1.0e3;
    r_process_info[SCALE_FACTOR] = 1.0;

    // First we create the nodes
    auto p_node1 = r_contact_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node2 = r_contact_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node3 = r_contact_model_part.CreateNewNode(3, 2.0 , 0.0 , 0.0);

    // Set flags
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Set(SLAVE);
        r_node.Set(ACTIVE);
    }

    // Set values
    p_node1->FastGetSolutionStepValue(WEIGHTED_GAP) = -1.0e-3;
    p_node1->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) = -1.0e3;
    p_node2->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
    p_node2->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) = 0.0;
    p_node3->FastGetSolutionStepValue(WEIGHTED_GAP) = 1.0e-3;
    p_node3->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) = -1.0e3;

    // Call utilities
    const auto is_converged = ActiveSetUtilities::ComputeALMFrictionlessActiveSet(r_model_part);

    // Check flags
    KRATOS_EXPECT_EQ(is_converged, 1);
    KRATOS_EXPECT_TRUE(p_node1->Is(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node2->IsNot(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node3->Is(ACTIVE));
}

/**
* Checks the correct work of the acttive set utilities
* Test ComputeALMFrictionlessComponentsActiveSet
*/
KRATOS_TEST_CASE_IN_SUITE(ComputeALMFrictionlessComponentsActiveSet, KratosContactStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);
    ModelPart& r_contact_model_part = r_model_part.CreateSubModelPart("Contact");

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
    r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
    r_model_part.AddNodalSolutionStepVariable(NORMAL);

    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;
    r_process_info[DELTA_TIME] = 1.0;
    r_process_info[INITIAL_PENALTY] = 1.0e3;
    r_process_info[SCALE_FACTOR] = 1.0;

    // First we create the nodes
    auto p_node1 = r_contact_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node2 = r_contact_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node3 = r_contact_model_part.CreateNewNode(3, 2.0 , 0.0 , 0.0);

    // Set flags
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Set(SLAVE);
        r_node.Set(ACTIVE);
    }

    // Set values
    p_node1->FastGetSolutionStepValue(WEIGHTED_GAP) = -1.0e-3;
    p_node1->FastGetSolutionStepValue(NORMAL_Y) = 1.0;
    p_node1->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_Y) = -1.0e3;
    p_node2->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
    p_node2->FastGetSolutionStepValue(NORMAL_Y) = 1.0;
    p_node2->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_Y) = 0.0;
    p_node3->FastGetSolutionStepValue(WEIGHTED_GAP) = 1.0e-3;
    p_node3->FastGetSolutionStepValue(NORMAL_Y) = 1.0;
    p_node3->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_Y) = -1.0e3;

    // Call utilities
    const auto is_converged = ActiveSetUtilities::ComputeALMFrictionlessComponentsActiveSet(r_model_part);

    // Check flags
    KRATOS_EXPECT_EQ(is_converged, 1);
    KRATOS_EXPECT_TRUE(p_node1->Is(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node2->IsNot(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node3->Is(ACTIVE));
}

/**
* Checks the correct work of the acttive set utilities
* Test ComputeALMFrictionalActiveSet
*/
KRATOS_TEST_CASE_IN_SUITE(ComputeALMFrictionalActiveSet, KratosContactStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);
    ModelPart& r_contact_model_part = r_model_part.CreateSubModelPart("Contact");

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
    r_model_part.AddNodalSolutionStepVariable(WEIGHTED_SLIP);
    r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
    r_model_part.AddNodalSolutionStepVariable(NORMAL);

    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;
    r_process_info[DELTA_TIME] = 1.0;
    r_process_info[INITIAL_PENALTY] = 1.0e3;
    r_process_info[SCALE_FACTOR] = 1.0;
    r_process_info[TANGENT_FACTOR] = 1.0;

    // First we create the nodes
    auto p_node1 = r_contact_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node2 = r_contact_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node3 = r_contact_model_part.CreateNewNode(3, 2.0 , 0.0 , 0.0);

    // Set flags
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Set(SLAVE);
        r_node.Set(ACTIVE);
        r_node.Set(SLIP);
        r_node.SetValue(FRICTION_COEFFICIENT, 1.0);
    }

    // Set values
    p_node1->FastGetSolutionStepValue(WEIGHTED_GAP) = -1.0e-3;

    p_node1->FastGetSolutionStepValue(WEIGHTED_SLIP_X) = -1.0e-3;
    p_node1->FastGetSolutionStepValue(NORMAL_Y) = 1.0;
    p_node1->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_Y) = -1.0e3;
    p_node2->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
    p_node2->FastGetSolutionStepValue(WEIGHTED_SLIP_X) = 0.0;
    p_node2->FastGetSolutionStepValue(NORMAL_Y) = 1.0;
    p_node2->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_Y) = 0.0;
    p_node3->FastGetSolutionStepValue(WEIGHTED_GAP) = 1.0e-3;
    p_node3->FastGetSolutionStepValue(WEIGHTED_SLIP_X) = 1.0e-3;
    p_node3->FastGetSolutionStepValue(NORMAL_Y) = 1.0;
    p_node3->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_Y) = -1.0e3;

    // Call utilities
    auto is_converged = ActiveSetUtilities::ComputeALMFrictionalActiveSet(r_model_part);

    // Check flags
    KRATOS_EXPECT_EQ(is_converged[0], 1);
    KRATOS_EXPECT_EQ(is_converged[1], 2);
    KRATOS_EXPECT_TRUE(p_node1->Is(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node1->IsNot(SLIP));
    KRATOS_EXPECT_TRUE(p_node2->IsNot(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node2->IsNot(SLIP));
    KRATOS_EXPECT_TRUE(p_node3->Is(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node3->IsNot(SLIP));

    // Set flags
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Set(ACTIVE);
        r_node.Set(SLIP);
    }

    // Call utilities (pure slip)
    is_converged = ActiveSetUtilities::ComputeALMFrictionalActiveSet(r_model_part, true);

    // Check flags
    KRATOS_EXPECT_EQ(is_converged[0], 1);
    KRATOS_EXPECT_EQ(is_converged[1], 0);
    KRATOS_EXPECT_TRUE(p_node1->Is(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node1->Is(SLIP));
    KRATOS_EXPECT_TRUE(p_node2->IsNot(ACTIVE));
    KRATOS_EXPECT_TRUE(!p_node2->IsDefined(SLIP));
    KRATOS_EXPECT_TRUE(p_node3->Is(ACTIVE));
    KRATOS_EXPECT_TRUE(p_node3->Is(SLIP));
}
}  // namespace Kratos::Testing.
