//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// Project includes
#include "testing/testing.h"
#include "future/containers/define_linear_algebra_serial.h"
#include "future/containers/linear_system.h"
#include "future/linear_operators/linear_operator.h"
#include "future/preconditioners/preconditioner.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(Preconditioner, KratosCoreFutureSuite)
{
    // Set up a fake linear system to test the preconditioner
    auto p_LHS = Kratos::make_shared<CsrMatrix<> >();
    DenseVector<double> aux_data(5);
    aux_data[0] = 3.0;
    aux_data[1] = 7.0;
    aux_data[2] = 2.0;
    aux_data[3] = 5.0;
    aux_data[4] = 1.0;
    SystemVector<> input(aux_data);
    SystemVector<> output(5);

    // Set a linear operator from the fake LHS matrix
    Future::LinearOperator<Future::SerialLinearAlgebraTraits>::UniquePointer p_lhs_lin_op = Kratos::make_unique<Future::SparseMatrixLinearOperator<Future::SerialLinearAlgebraTraits> >(p_LHS);

    // Set up the preconditioner
    Future::Preconditioner<Future::SerialLinearAlgebraTraits> preconditioner;

    // Check preconditioner features
    KRATOS_EXPECT_FALSE(preconditioner.RequiresAdditionalData());
    KRATOS_EXPECT_FALSE(preconditioner.HasAdditionalData());

    // Call the preconditioner sequence
    preconditioner.Initialize(p_lhs_lin_op);
    preconditioner.InitializeSolutionStep(p_lhs_lin_op);
    preconditioner.Apply(input, output);
    preconditioner.FinalizeSolutionStep(p_lhs_lin_op);
    preconditioner.Clear();

    // Check the obtained results
    std::vector<double> ref_output = {3.0, 7.0, 2.0, 5.0, 1.0};
    KRATOS_EXPECT_VECTOR_NEAR(output.data(), ref_output, 1e-12);
}

}