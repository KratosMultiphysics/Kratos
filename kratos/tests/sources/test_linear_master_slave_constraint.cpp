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
#include "includes/model_part.h"
#include "includes/linear_master_slave_constraint.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(LinearMasterSlaveConstraintTests, KratosCoreFastSuite)
{
        ModelPart model_part("main");
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


        KRATOS_CHECK_EQUAL(master_vector.size(), 1);
        KRATOS_CHECK_EQUAL(slave_vector.size(), 1);

        KRATOS_CHECK_EQUAL(transformation_matrix.size1(), 1);
        KRATOS_CHECK_EQUAL(transformation_matrix.size2(), 1);
        KRATOS_CHECK_EQUAL(constant_vector.size(), 1);

        KRATOS_CHECK_DOUBLE_EQUAL(transformation_matrix(0,0), 1.0); // TODO: Check -> comparison between two doubles ??
        KRATOS_CHECK_DOUBLE_EQUAL(constant_vector(0), 0.0);
}

}
}  // namespace Kratos.


