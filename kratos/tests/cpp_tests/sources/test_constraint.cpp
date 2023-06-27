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
#include "constraints/linear_master_slave_constraint.h"

namespace Kratos {
    namespace Testing {
        typedef Node NodeType;
        typedef std::size_t IndexType;

        /**
        *  Here the clone operator is test
        */
        KRATOS_TEST_CASE_IN_SUITE(ConstraintCloneOperator, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("test");

            r_model_part.SetBufferSize(3);
            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(REACTION);

            // Definition of nodes
            auto p_node1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            auto p_node2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);

            for (auto& r_node : r_model_part.Nodes()) {
                r_node.AddDof(DISPLACEMENT_X, REACTION_X);
                r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
                r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
            }

            // Create constraint
            auto p_const = r_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *p_node1, DISPLACEMENT_X, *p_node2, DISPLACEMENT_X, 1.0, 0);

            p_const->SetValue(DISTANCE, 12.1);
            p_const->SetValue(VELOCITY_X, 32.4);
            p_const->Set(ACTIVE, true);

            auto p_clone_of_const = p_const->Clone(2);

            KRATOS_CHECK_EQUAL(p_clone_of_const->Id(), 2);
            KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_const->GetValue(DISTANCE), 12.1);
            KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_const->GetValue(VELOCITY_X), 32.4);
            KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_const->GetValue(VELOCITY_Y), 0.00);
            KRATOS_CHECK_DOUBLE_EQUAL(p_clone_of_const->GetValue(VELOCITY_Z), 0.00);
            KRATOS_CHECK(p_clone_of_const->Is(ACTIVE));
        }

        /**
        *  Here the several operators are test
        */
        KRATOS_TEST_CASE_IN_SUITE(ConstraintSeveralOperations, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& r_model_part = current_model.CreateModelPart("test");

            r_model_part.SetBufferSize(3);
            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(REACTION);

            // Definition of nodes
            auto p_node1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            auto p_node2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);

            for (auto& r_node : r_model_part.Nodes()) {
                r_node.AddDof(DISPLACEMENT_X, REACTION_X);
                r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
                r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
            }

            // Create constraint
            auto p_const = r_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *p_node1, DISPLACEMENT_X, *p_node2, DISPLACEMENT_X, 1.0, 0.1);

            Matrix matrix;
            Vector vector;
            auto& r_process_info = r_model_part.GetProcessInfo();
            p_const->GetLocalSystem(matrix, vector, r_process_info);

            const auto& r_slave_dofs = p_const->GetSlaveDofsVector();
            const auto& r_master_dofs = p_const->GetMasterDofsVector();

            // Create empty constraint
            MasterSlaveConstraint::Pointer p_new_const = Kratos::make_shared<LinearMasterSlaveConstraint>(2);
            p_new_const->SetLocalSystem(matrix, vector, r_process_info);
            p_new_const->SetDofList(r_slave_dofs, r_master_dofs, r_process_info);

            Matrix new_matrix;
            Vector new_vector;
            p_new_const->GetLocalSystem(new_matrix, new_vector, r_process_info);

            for (IndexType i = 0; i < vector.size(); ++i) {
                KRATOS_CHECK_NEAR(vector[i], new_vector[i], 1.0e-12);
            }
            for (IndexType i = 0; i < matrix.size1(); ++i) {
                for (IndexType j = 0; j < matrix.size2(); ++j) {
                    KRATOS_CHECK_NEAR(matrix(i, j), new_matrix(i, j), 1.0e-12);
                }
            }
        }
    }  // namespace Testing.
}  // namespace Kratos.
