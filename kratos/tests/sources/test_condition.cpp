//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "containers/model.h"


namespace Kratos {
    namespace Testing {
        typedef Node<3> NodeType;
        typedef std::size_t IndexType;

        /**
        *  Here the clone operator is test
        */
        KRATOS_TEST_CASE_IN_SUITE(ConditionCloneOperator, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& model_part = current_model.CreateModelPart("test");

            // Definition of nodes
            auto p_node1 = model_part.CreateNewNode(1, 1.0, 0.0, 0.0);
            auto p_node2 = model_part.CreateNewNode(2, 1.0, 1.0, 0.0);
            auto p_node3 = model_part.CreateNewNode(3, 1.0, 1.0, 0.0);

            // Definition of properties
            auto p_prop = model_part.pGetProperties(1);

            // List onf nodes
            std::vector<NodeType::Pointer> list_nodes(3);
            list_nodes[0] = p_node1;
            list_nodes[1] = p_node2;
            list_nodes[2] = p_node3;

            auto p_cond = model_part.CreateNewCondition("Condition3D", 1, PointerVector<NodeType>{list_nodes}, p_prop);

            p_cond->SetValue(DISTANCE, 12.1);
            p_cond->SetValue(VELOCITY_X, 32.4);
            p_cond->Set(ACTIVE, true);

            Condition::Pointer p_clone_of_cond = p_cond->Clone(2, PointerVector<NodeType>{list_nodes});

            KRATOS_CHECK_EQUAL(p_clone_of_cond->Id(), 2);
            KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_cond->GetValue(DISTANCE), 12.1);
            KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_cond->GetValue(VELOCITY_X), 32.4);
            KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_cond->GetValue(VELOCITY_Y), 0.00);
            KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_cond->GetValue(VELOCITY_Z), 0.00);
            KRATOS_CHECK(p_clone_of_cond->Is(ACTIVE));
        }
    }  // namespace Testing.
}  // namespace Kratos.
