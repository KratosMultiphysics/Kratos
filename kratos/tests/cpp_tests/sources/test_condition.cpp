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
            auto p_prop = model_part.CreateNewProperties(1);

            auto p_cond = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, {{1, 2, 3}}, p_prop);

            p_cond->SetValue(DISTANCE, 12.1);
            p_cond->SetValue(VELOCITY, ZeroVector(3));
            p_cond->SetValue(VELOCITY_X, 32.4);
            p_cond->Set(ACTIVE, true);

            Condition::Pointer p_clone_of_cond = p_cond->Clone(2, p_cond->GetGeometry());

            KRATOS_EXPECT_EQ(p_clone_of_cond->Id(), 2);
            KRATOS_EXPECT_DOUBLE_EQ(p_clone_of_cond->GetValue(DISTANCE), 12.1);
            KRATOS_EXPECT_DOUBLE_EQ(p_clone_of_cond->GetValue(VELOCITY_X), 32.4);
            KRATOS_EXPECT_DOUBLE_EQ(p_clone_of_cond->GetValue(VELOCITY_Y), 0.0);
            KRATOS_EXPECT_DOUBLE_EQ(p_clone_of_cond->GetValue(VELOCITY_Z), 0.0);
            KRATOS_EXPECT_TRUE(p_clone_of_cond->Is(ACTIVE));
        }
    }  // namespace Testing.
}  // namespace Kratos.
