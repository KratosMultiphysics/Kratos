//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Vicente Mataix Ferrandiz
//
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "includes/stream_serializer.h"

namespace Kratos {
    namespace Testing {

    typedef Node<3> NodeType;

    /**
     *  Here the assigment operator is test
     */

    KRATOS_TEST_CASE_IN_SUITE(NodeAssignOperator, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("test");
        model_part.AddNodalSolutionStepVariable(DISTANCE);
        model_part.AddNodalSolutionStepVariable(VELOCITY);

        auto p_node = model_part.CreateNewNode(1, 1., 0, 0);

        p_node->FastGetSolutionStepValue(DISTANCE) = 12.1;
        p_node->FastGetSolutionStepValue(VELOCITY_X) = 32.4;
        p_node->Set(ACTIVE, true);

        NodeType copy_of_node(2,1,0,0);
        copy_of_node = *p_node;

        KRATOS_CHECK_EQUAL(copy_of_node.Id(), 1);
        KRATOS_CHECK_DOUBLE_EQUAL(copy_of_node.FastGetSolutionStepValue(DISTANCE), 12.1);
        KRATOS_CHECK_DOUBLE_EQUAL(copy_of_node.FastGetSolutionStepValue(VELOCITY_X), 32.4);
        KRATOS_CHECK_DOUBLE_EQUAL(copy_of_node.FastGetSolutionStepValue(VELOCITY_Y), 0.00);
        KRATOS_CHECK_DOUBLE_EQUAL(copy_of_node.FastGetSolutionStepValue(VELOCITY_Z), 0.00);
        KRATOS_CHECK(copy_of_node.Is(ACTIVE));
    }

    /**
     *  Here the clone operator is test
     */
    KRATOS_TEST_CASE_IN_SUITE(NodeCloneOperator, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("test");
        model_part.AddNodalSolutionStepVariable(DISTANCE);
        model_part.AddNodalSolutionStepVariable(VELOCITY);

        auto p_node = model_part.CreateNewNode(1, 1., 0, 0);

        p_node->FastGetSolutionStepValue(DISTANCE) = 12.1;
        p_node->FastGetSolutionStepValue(VELOCITY_X) = 32.4;
        p_node->Set(ACTIVE, true);

        NodeType::Pointer p_clone_of_node = p_node->Clone();

        KRATOS_CHECK_EQUAL(p_clone_of_node->Id(), 1);
        KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_node->FastGetSolutionStepValue(DISTANCE), 12.1);
        KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_node->FastGetSolutionStepValue(VELOCITY_X), 32.4);
        KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_node->FastGetSolutionStepValue(VELOCITY_Y), 0.00);
        KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_node->FastGetSolutionStepValue(VELOCITY_Z), 0.00);
        KRATOS_CHECK(p_clone_of_node->Is(ACTIVE));
    }

    /**
     *  Here the DoF serialization
     */
    KRATOS_TEST_CASE_IN_SUITE(NodeSerialization, KratosCoreFastSuite) 
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("test");
        r_model_part.AddNodalSolutionStepVariable(DISTANCE);
        r_model_part.AddNodalSolutionStepVariable(VELOCITY);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(REACTION);

        auto p_node_to_be_saved = r_model_part.CreateNewNode(1, 1., 0, 0);
        auto p_node_to_be_loaded = Node<3>::Pointer(nullptr);

        p_node_to_be_saved->Fix(DISTANCE);
        p_node_to_be_saved->AddDof(DISPLACEMENT_X, REACTION_X);

        StreamSerializer serializer;

        serializer.save("NodalData", p_node_to_be_saved);
        serializer.load("NodalData", p_node_to_be_loaded);

        KRATOS_CHECK_EQUAL(p_node_to_be_saved->Id(), p_node_to_be_loaded->Id());
        KRATOS_CHECK(p_node_to_be_loaded->IsFixed(DISTANCE));
        Dof<double>& dof = p_node_to_be_loaded->AddDof(DISTANCE);
        dof.GetSolutionStepValue() = 2.345;
        KRATOS_CHECK_EQUAL(p_node_to_be_loaded->GetSolutionStepValue(DISTANCE), 2.345);
        KRATOS_CHECK(p_node_to_be_loaded->GetDof(DISPLACEMENT_X).IsFree());
        KRATOS_CHECK(p_node_to_be_loaded->GetDof(DISPLACEMENT_X).HasReaction());
        KRATOS_CHECK_EQUAL(p_node_to_be_loaded->GetDof(DISPLACEMENT_X).GetReaction(), REACTION_X);
    }
    
    /**
     *  Here the pGetDof position
     */
    KRATOS_TEST_CASE_IN_SUITE(NodepGetDofPosition, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("test");
        model_part.AddNodalSolutionStepVariable(VELOCITY);

        auto p_node = model_part.CreateNewNode(1, 1., 0, 0);

        p_node->AddDof(VELOCITY_X);
        p_node->AddDof(VELOCITY_Y);
        p_node->AddDof(VELOCITY_Z);

        // GetDofPosition
        KRATOS_CHECK_EQUAL(p_node->GetDofPosition(VELOCITY_X), 0);
        KRATOS_CHECK_EQUAL(p_node->GetDofPosition(VELOCITY_Y), 1);
        KRATOS_CHECK_EQUAL(p_node->GetDofPosition(VELOCITY_Z), 2);

        // GetDof with position
        KRATOS_CHECK((p_node->pGetDof(VELOCITY_X, 0))->GetVariable() == VELOCITY_X);
        KRATOS_CHECK((p_node->pGetDof(VELOCITY_Y, 1))->GetVariable() == VELOCITY_Y);
        KRATOS_CHECK((p_node->pGetDof(VELOCITY_Z, 2))->GetVariable() == VELOCITY_Z);
    }
}  // namespace Testing.
}  // namespace Kratos.
