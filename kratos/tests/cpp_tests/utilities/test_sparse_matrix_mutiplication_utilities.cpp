//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "utilities/sparse_matrix_multiplication_utility.h"

namespace Kratos {
namespace Testing {

/// The sparse matrix type
typedef typename UblasSpace<double, CompressedMatrix, Vector>::MatrixType SparseMatrixType;

KRATOS_TEST_CASE_IN_SUITE(AssembleSparseMatrixByBlocks, KratosCoreFastSuite)
{
    SparseMatrixType identity2x2(2, 2);
    for (IndexType i = 0; i < 2; ++i) {
        identity2x2.push_back(i, i, 1.0);
    }

    DenseMatrix<SparseMatrixType*> matrices_p_blocks(2,2);
    matrices_p_blocks(0,0) = &identity2x2;
    matrices_p_blocks(1,0) = &identity2x2;
    matrices_p_blocks(0,1) = &identity2x2;
    matrices_p_blocks(1,1) = &identity2x2;
    DenseMatrix<double> contribution_coefficients(2,2);
    contribution_coefficients(0,0) = 1.0;
    contribution_coefficients(1,0) = -1.0;
    contribution_coefficients(0,1) = -1.0;
    contribution_coefficients(1,1) = 1.0;

    SparseMatrixType solution_matrix;
    SparseMatrixMultiplicationUtility::AssembleSparseMatrixByBlocks(solution_matrix, matrices_p_blocks, contribution_coefficients);

    const double tolerance = 1.0e-16;
    KRATOS_EXPECT_NEAR(solution_matrix(0,0), 1.0, tolerance);
    KRATOS_EXPECT_NEAR(solution_matrix(1,1), 1.0, tolerance);
    KRATOS_EXPECT_NEAR(solution_matrix(2,2), 1.0, tolerance);
    KRATOS_EXPECT_NEAR(solution_matrix(3,3), 1.0, tolerance);
    KRATOS_EXPECT_NEAR(solution_matrix(2,0), -1.0, tolerance);
    KRATOS_EXPECT_NEAR(solution_matrix(3,1), -1.0, tolerance);
    KRATOS_EXPECT_NEAR(solution_matrix(0,2), -1.0, tolerance);
    KRATOS_EXPECT_NEAR(solution_matrix(1,3), -1.0, tolerance);

    double total = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            total += solution_matrix(i,j);
        }
    }
    KRATOS_EXPECT_NEAR(total, 0.0, tolerance);
}

// Test whether products of matrices with different types compile.
KRATOS_TEST_CASE_IN_SUITE(HeterogeneousProduct, KratosCoreFastSuite)
{
    using SparseUtils = SparseMatrixMultiplicationUtility;
    constexpr int system_size = 2;
    constexpr double double_tolerance = 1e-16;
    constexpr float single_tolerance = 1e-7;

    // Make diagonal matrices as operands.
    TUblasSparseSpace<double>::MatrixType double_precision(system_size, system_size, system_size);
    TUblasSparseSpace<float>::MatrixType single_precision(system_size, system_size, system_size);

    double_precision.index1_data()[0] = 0;
    single_precision.index1_data()[0] = 0;
    for (int i_row=0; i_row<system_size; ++i_row) {
        double_precision.index1_data()[i_row + 1] = i_row + 1;
        single_precision.index1_data()[i_row + 1] = i_row + 1;

        double_precision.index2_data()[i_row] = i_row;
        single_precision.index2_data()[i_row] = i_row;

        double_precision.value_data()[i_row] = 2 * i_row;
        single_precision.value_data()[i_row] = 3 * i_row;
    } // for i_row in range(system_size)
    double_precision.set_filled(system_size + 1, system_size);
    single_precision.set_filled(system_size + 1, system_size);

    // Output matrix has identical type as the left operand.
    {
        TUblasSparseSpace<double>::MatrixType product;
        SparseUtils::MatrixMultiplication(double_precision, single_precision, product);
        for (int i_row=0; i_row<system_size; ++i_row) {
            KRATOS_EXPECT_EQ(product.index1_data()[i_row], i_row);
            KRATOS_EXPECT_EQ(product.index2_data()[i_row], i_row);
            KRATOS_EXPECT_NEAR(product.value_data()[i_row], double(2 * i_row) * float(3 * i_row), double_tolerance);
        }
        KRATOS_EXPECT_EQ(product.index1_data()[system_size], system_size);
    }

    // Output matrix has identical type as the right operand.
    {
        TUblasSparseSpace<float>::MatrixType product;
        SparseUtils::MatrixMultiplication(double_precision, single_precision, product);
        for (int i_row=0; i_row<system_size; ++i_row) {
            KRATOS_EXPECT_EQ(product.index1_data()[i_row], i_row);
            KRATOS_EXPECT_EQ(product.index2_data()[i_row], i_row);
            KRATOS_EXPECT_NEAR(product.value_data()[i_row], double(2 * i_row) * float(3 * i_row), single_tolerance);
        }
        KRATOS_EXPECT_EQ(product.index1_data()[system_size], system_size);
    }

    // Both operands as well as the result have different types.
    {
        TUblasSparseSpace<int>::MatrixType product;
        SparseUtils::MatrixMultiplication(double_precision, single_precision, product);
        for (int i_row=0; i_row<system_size; ++i_row) {
            KRATOS_EXPECT_EQ(product.index1_data()[i_row], i_row);
            KRATOS_EXPECT_EQ(product.index2_data()[i_row], i_row);
            KRATOS_EXPECT_EQ(product.value_data()[i_row], double(2 * i_row) * float(3 * i_row));
        }
        KRATOS_EXPECT_EQ(product.index1_data()[system_size], system_size);
    }
}

}   // namespace Testing
}  // namespace Kratos.
