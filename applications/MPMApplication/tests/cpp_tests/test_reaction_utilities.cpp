//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Crescenzio
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"

// Application includes
#include "testing/testing.h"
#include "custom_utilities/reaction_utilities.h"
#include "mpm_application_variables.h"

namespace Kratos::Testing {

    /**
     * Compute sum of REACTION variable over a set of nodes
     */
    KRATOS_TEST_CASE_IN_SUITE(MPMComputeGridConformingReaction, KratosMPMFastSuite)
    {
        Model model;

        ModelPart& model_part = model.CreateModelPart("GridModelPart");

        model_part.AddNodalSolutionStepVariable(REACTION);

        Properties::Pointer p_prop = model_part.CreateNewProperties(0);

        auto p_node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        auto p_node_3 = model_part.CreateNewNode(3, 2.0, 0.0, 0.0);
        auto p_node_4 = model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
        auto p_node_5 = model_part.CreateNewNode(5, 1.0, 1.0, 0.0);
        auto p_node_6 = model_part.CreateNewNode(6, 2.0, 1.0, 0.0);

        model_part.CreateNewElement("Element2D4N", 1, {1, 2, 5, 4}, p_prop);
        model_part.CreateNewElement("Element2D4N", 2, {2, 3, 6, 5}, p_prop);

        model_part.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_prop);
        model_part.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_prop);
        model_part.CreateNewCondition("LineCondition2D2N", 3, {{3, 6}}, p_prop);

        ModelPart& sub_model_part = model_part.CreateSubModelPart("ReactionModelPart");
        sub_model_part.AddNodes({1, 2, 3, 6});
        sub_model_part.AddConditions({1, 2, 3});

        model_part.GetNode(1).FastGetSolutionStepValue(REACTION_X) = 5.0;
        model_part.GetNode(1).FastGetSolutionStepValue(REACTION_Y) = 10.0;
        model_part.GetNode(1).FastGetSolutionStepValue(REACTION_Z) = 15.0;

        model_part.GetNode(2).FastGetSolutionStepValue(REACTION_X) = 2.0;
        model_part.GetNode(2).FastGetSolutionStepValue(REACTION_Y) = 4.0;
        model_part.GetNode(2).FastGetSolutionStepValue(REACTION_Z) = 6.0;

        model_part.GetNode(3).FastGetSolutionStepValue(REACTION_X) = 1.0;
        model_part.GetNode(3).FastGetSolutionStepValue(REACTION_Y) = 1.0;
        model_part.GetNode(3).FastGetSolutionStepValue(REACTION_Z) = 1.0;

        model_part.GetNode(4).FastGetSolutionStepValue(REACTION_X) = 10.0;
        model_part.GetNode(4).FastGetSolutionStepValue(REACTION_Y) = 4.0;
        model_part.GetNode(4).FastGetSolutionStepValue(REACTION_Z) = 8.0;

        model_part.GetNode(6).FastGetSolutionStepValue(REACTION_X) = 0.5;
        model_part.GetNode(6).FastGetSolutionStepValue(REACTION_Y) = 8.0;
        model_part.GetNode(6).FastGetSolutionStepValue(REACTION_Z) = 2.5;

        array_1d<double, 3> reaction_force = ReactionUtilities::CalculateGridConformingReaction(model.GetModelPart("GridModelPart.ReactionModelPart"));

        KRATOS_EXPECT_NEAR(reaction_force[0],  8.5, 1e-6);
        KRATOS_EXPECT_NEAR(reaction_force[1], 23.0, 1e-6);
        KRATOS_EXPECT_NEAR(reaction_force[2], 24.5, 1e-6);
    }

    /**
     * Checks the embedded drag computation utility.
     */
    KRATOS_TEST_CASE_IN_SUITE(MPMComputeNonConformingReaction, KratosMPMFastSuite)
    {
        Model model;

        ModelPart& model_part = model.CreateModelPart("GridModelPart");

        model_part.AddNodalSolutionStepVariable(REACTION);

        Properties::Pointer p_prop = model_part.CreateNewProperties(0);

        auto p_node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        auto p_node_3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        auto p_node_4 = model_part.CreateNewNode(4, 1.0, 1.0, 0.0);

        model_part.CreateNewElement("Element2D4N", 1, {1, 2, 4, 3}, p_prop);

        ModelPart& mpm_model_part = model.CreateModelPart("MPMModelPart");

        array_1d<double, 3> mp_coord_1{ 0.3, 0.3, 0.0 };
        array_1d<double, 3> mp_coord_2{ 0.3, 0.7, 0.0 };
        array_1d<double, 3> mp_coord_3{ 0.5, 0.5, 0.0 };
        array_1d<double, 3> mp_coord_4{ 0.7, 0.3, 0.0 };
        array_1d<double, 3> mp_coord_5{ 0.7, 0.7, 0.0 };

        auto p_quad_1 = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(
            model_part.GetElement(1).pGetGeometry(), mp_coord_1, 1.0);
        auto p_quad_2 = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(
            model_part.GetElement(1).pGetGeometry(), mp_coord_2, 1.0);
        auto p_quad_3 = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(
            model_part.GetElement(1).pGetGeometry(), mp_coord_3, 1.0);
        auto p_quad_4 = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(
            model_part.GetElement(1).pGetGeometry(), mp_coord_4, 1.0);
        auto p_quad_5 = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(
            model_part.GetElement(1).pGetGeometry(), mp_coord_5, 1.0);

        auto p_cond_1 = mpm_model_part.CreateNewCondition("MPMParticlePenaltyDirichletCondition", 1, p_quad_1, p_prop);
        auto p_cond_2 = mpm_model_part.CreateNewCondition("MPMParticlePenaltyDirichletCondition", 2, p_quad_2, p_prop);
        auto p_cond_3 = mpm_model_part.CreateNewCondition("MPMParticlePenaltyDirichletCondition", 3, p_quad_3, p_prop);
        auto p_cond_4 = mpm_model_part.CreateNewCondition("MPMParticlePenaltyDirichletCondition", 4, p_quad_4, p_prop);
        auto p_cond_5 = mpm_model_part.CreateNewCondition("MPMParticlePenaltyDirichletCondition", 5, p_quad_5, p_prop);

        ModelPart& sub_model_part = mpm_model_part.CreateSubModelPart("MPMReactionModelPart");
        sub_model_part.AddConditions({1, 2, 4});

        array_1d<double, 3> mpc_force_1{ 10.0, 0.5,  0.0 };
        array_1d<double, 3> mpc_force_2{  4.0, 5.0,  6.0 };
        array_1d<double, 3> mpc_force_3{  1.0, 3.0,  0.7 };
        array_1d<double, 3> mpc_force_4{  8.3, 7.0,  3.0 };
        array_1d<double, 3> mpc_force_5{  2.0, 1.0, 12.0 };

        const auto& process_info = mpm_model_part.GetProcessInfo();
        p_cond_1->SetValuesOnIntegrationPoints(MPC_CONTACT_FORCE, { mpc_force_1 }, process_info);
        p_cond_2->SetValuesOnIntegrationPoints(MPC_CONTACT_FORCE, { mpc_force_2 }, process_info);
        p_cond_3->SetValuesOnIntegrationPoints(MPC_CONTACT_FORCE, { mpc_force_3 }, process_info);
        p_cond_4->SetValuesOnIntegrationPoints(MPC_CONTACT_FORCE, { mpc_force_4 }, process_info);
        p_cond_5->SetValuesOnIntegrationPoints(MPC_CONTACT_FORCE, { mpc_force_5 }, process_info);

        array_1d<double, 3> reaction_force = ReactionUtilities::CalculateNonConformingReaction(model.GetModelPart("MPMModelPart.MPMReactionModelPart"));

        KRATOS_EXPECT_NEAR(reaction_force[0], 22.3, 1e-6);
        KRATOS_EXPECT_NEAR(reaction_force[1], 12.5, 1e-6);
        KRATOS_EXPECT_NEAR(reaction_force[2],  9.0, 1e-6);
    }

}  // namespace Kratos::Testing
