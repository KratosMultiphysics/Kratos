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
#include "future/containers/linear_operator.h"
#include "future/containers/sparse_matrix_linear_operator.h"
#include "containers/system_vector.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(LinearOperatorWithCsr, KratosCoreFutureSuite)
{
    // Set the input and output vectors
    DenseVector<double> aux_data(5);
    aux_data[0] = 3.0;
    aux_data[1] = 7.0;
    aux_data[2] = 2.0;
    aux_data[3] = 5.0;
    aux_data[4] = 1.0;
    SystemVector<> input_vector(aux_data);
    SystemVector<> output_vector(5);
    output_vector.SetValue(0.0);

    // Set up the CSR matrix
    typename CsrMatrix<double>::MatrixMapType matrix_map;
    matrix_map[{0, 0}] = 10; matrix_map[{0, 3}] = -2;
    matrix_map[{1, 1}] = 8; matrix_map[{1, 2}] = -1;
    matrix_map[{2, 2}] = 5; matrix_map[{2, 1}] = -1; matrix_map[{2, 4}] = -2;
    matrix_map[{3, 3}] = 7; matrix_map[{3, 0}] = -2; matrix_map[{3, 4}] = -1;
    matrix_map[{4, 4}] = 6; matrix_map[{4, 2}] = -2; matrix_map[{4, 3}] = -1;
    CsrMatrix<double> csr_matrix(matrix_map);

    // Set up the linear operator from the CSR matrix
    const Future::SparseMatrixLinearOperator<SystemVector<>,CsrMatrix<>> linear_operator(csr_matrix);

    // Apply the linear operator to an input vector
    linear_operator.SpMV(input_vector, output_vector);

    // Check linear operator features
    KRATOS_EXPECT_EQ(linear_operator.NumRows(), 5);
    KRATOS_EXPECT_EQ(linear_operator.NumCols(), 5);
    KRATOS_EXPECT_FALSE(linear_operator.IsMatrixFree());

    // Check the output values
    KRATOS_EXPECT_NEAR(output_vector[0], 20.0, 1e-12);
    KRATOS_EXPECT_NEAR(output_vector[1], 54.0, 1e-12);
    KRATOS_EXPECT_NEAR(output_vector[2], 1.0, 1e-12);
    KRATOS_EXPECT_NEAR(output_vector[3], 28.0, 1e-12);
    KRATOS_EXPECT_NEAR(output_vector[4], -3.0, 1e-12);
}

}