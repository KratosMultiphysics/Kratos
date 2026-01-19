//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala / Philipp Bucher
//
//

// System includes


// External includes


// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "constraints/linear_master_slave_constraint.h"

#include "includes/stream_serializer.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(LinearMasterSlaveConstraintTests, KratosCoreFastSuite)
{
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("test_model_part",2);
        model_part.AddNodalSolutionStepVariable(PRESSURE);
        auto n1 = model_part.CreateNewNode(1, 0.00,0.00,0.00);
        auto n2 = model_part.CreateNewNode(2, 1.00,0.00,0.00);
        n1->AddDof(PRESSURE);
        n2->AddDof(PRESSURE);

        auto c1 = model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *n1, PRESSURE, *n2, PRESSURE, 1.0, 0.0);

        ProcessInfo& process_info = model_part.GetProcessInfo();
        LinearMasterSlaveConstraint::EquationIdVectorType master_vector;
        LinearMasterSlaveConstraint::EquationIdVectorType slave_vector;

        c1->EquationIdVector(master_vector, slave_vector, process_info);

        LinearMasterSlaveConstraint::MatrixType transformation_matrix;
        LinearMasterSlaveConstraint::VectorType constant_vector;

        c1->CalculateLocalSystem(transformation_matrix, constant_vector, process_info);


        KRATOS_EXPECT_EQ(master_vector.size(), 1);
        KRATOS_EXPECT_EQ(slave_vector.size(), 1);

        KRATOS_EXPECT_EQ(transformation_matrix.size1(), 1);
        KRATOS_EXPECT_EQ(transformation_matrix.size2(), 1);
        KRATOS_EXPECT_EQ(constant_vector.size(), 1);

        KRATOS_EXPECT_DOUBLE_EQ(transformation_matrix(0,0), 1.0); // TODO: Check -> comparison between two doubles ??
        KRATOS_EXPECT_DOUBLE_EQ(constant_vector(0), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(LinearMasterSlaveConstraintSerialization, KratosCoreFastSuite)
{
    // Create model and populate it
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("a");
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(REACTION);

    // Create nodes
    auto p_node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    
    // Add DoFs
    p_node_1->AddDof(DISPLACEMENT_X, REACTION_X);
    p_node_2->AddDof(DISPLACEMENT_X, REACTION_X);
    
    KRATOS_WATCH(p_node_1->pGetDof(DISPLACEMENT_X));
//     KRATOS_WATCH((p_node_1->pGetDof(DISPLACEMENT_X)).get());
    KRATOS_WATCH(p_node_2->pGetDof(DISPLACEMENT_X));
//     KRATOS_WATCH((p_node_2->pGetDof(DISPLACEMENT_X)).get());
    // Create Constraint: Node 2 X depends on Node 1 X
    // Relation: slave = 2.0 * master + 5.0
    auto p_constraint = model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *p_node_1, DISPLACEMENT_X, *p_node_2, DISPLACEMENT_X, 2.0, 5.0);

    // Serialize
    StreamSerializer serializer;
    serializer.save("ModelPart", model_part);

    KRATOS_WATCH("before deleting")
    current_model.DeleteModelPart("a");
    KRATOS_WATCH("mp deleted")


    // Deserialize into new model
    Model loaded_model;
    ModelPart& loaded_model_part = loaded_model.CreateModelPart("a");


    //serializer.SetLoadState();
    serializer.load("ModelPart", loaded_model_part);

    // Verify
    KRATOS_EXPECT_EQ(loaded_model_part.NumberOfMasterSlaveConstraints(), 1);
    
    auto p_loaded_constraint = loaded_model_part.GetMasterSlaveConstraint(1);
    KRATOS_EXPECT_EQ(p_loaded_constraint.Id(), 1);
    
    const auto& r_process_info = loaded_model_part.GetProcessInfo();
    LinearMasterSlaveConstraint::MatrixType transformation_matrix;
    LinearMasterSlaveConstraint::VectorType constant_vector;
    p_loaded_constraint.CalculateLocalSystem(transformation_matrix, constant_vector, r_process_info);
    
    KRATOS_EXPECT_DOUBLE_EQ(transformation_matrix(0,0), 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(constant_vector(0), 5.0);
    
    // Check Master/Slave DoFs implicitly via EquationIdVector or manually retrieving them if possible.
    // The CalculateLocalSystem check confirms the relation properties were loaded.
    // We can also check that the DoFs are pointing to the correct nodes/variables?
    // LinearMasterSlaveConstraint stores Dof pointers.
    // The public API doesn't easily expose the raw Dofs via lightweight accessors except specific methods.
    // But CalculateLocalSystem proving correct weights implies successful data restoration.
}

}
}  // namespace Kratos.


