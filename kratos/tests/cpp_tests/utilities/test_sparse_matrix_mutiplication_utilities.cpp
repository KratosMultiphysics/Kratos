//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
    KRATOS_CHECK_NEAR(solution_matrix(0,0), 1.0, tolerance);
    KRATOS_CHECK_NEAR(solution_matrix(1,1), 1.0, tolerance);
    KRATOS_CHECK_NEAR(solution_matrix(2,2), 1.0, tolerance);
    KRATOS_CHECK_NEAR(solution_matrix(3,3), 1.0, tolerance);
    KRATOS_CHECK_NEAR(solution_matrix(2,0), -1.0, tolerance);
    KRATOS_CHECK_NEAR(solution_matrix(3,1), -1.0, tolerance);
    KRATOS_CHECK_NEAR(solution_matrix(0,2), -1.0, tolerance);
    KRATOS_CHECK_NEAR(solution_matrix(1,3), -1.0, tolerance);

    double total = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            total += solution_matrix(i,j);
        }
    }
    KRATOS_CHECK_NEAR(total, 0.0, tolerance);
}

}   // namespace Testing
}  // namespace Kratos.
