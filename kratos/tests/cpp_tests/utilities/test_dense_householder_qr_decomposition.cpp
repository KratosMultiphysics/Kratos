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

// System includes
#include <limits>

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "utilities/dense_householder_qr_decomposition.h"
#include "testing/testing.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(DenseHouseholderQRDecomposition, KratosCoreFastSuite)
{
    // Allocate matrices
    Matrix A_matrix(2,2);
    Matrix Q_matrix(2,2);
    Matrix R_matrix(2,2);

    // Set the matrix to be decomposed
    A_matrix(0,0) = 0.57690;
    A_matrix(0,1) = 0.28760;
    A_matrix(1,0) = 0.72886;
    A_matrix(1,1) = 0.40541;
    Matrix A_copy = A_matrix; // Note that A will be modified. We keep a copy for the testing.

    // Calculate the QR decomposition
    using DenseSpace = UblasSpace<double, Matrix, Vector>;
    DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
    qr_decomposition.Compute(A_matrix, Q_matrix, R_matrix);

    // Check decomposition is correct
    constexpr double tolerance = 1e-10;
    const Matrix QR_matrix = prod(Q_matrix, R_matrix);
    KRATOS_EXPECT_MATRIX_NEAR(QR_matrix, A_copy, tolerance);

    // Check values
    KRATOS_EXPECT_NEAR(R_matrix(0,1), -0.49637670012, tolerance);
    KRATOS_EXPECT_NEAR(R_matrix(1,1), 0.0260998022652, tolerance);
    KRATOS_EXPECT_NEAR(Q_matrix(0,0), -0.620627440497, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DenseHouseholderQRSolveVector, KratosCoreFastSuite)
{
    // Allocate matrices
    Vector x_vector(4);
    Vector b_vector(4);
    Matrix A_matrix(4,4);

    // Set the system to be solved
    A_matrix(0, 0) = 0.0;
    A_matrix(0, 1) = 0.979749;
    A_matrix(0, 2) = 0.494393;
    A_matrix(0, 3) = 0.23073;
    A_matrix(1, 0) = 1.79224;
    A_matrix(1, 1) = 0.198842;
    A_matrix(1, 2) = 0.074485;
    A_matrix(1, 3) = 1.45717;
    A_matrix(2, 0) = 1.6039;
    A_matrix(2, 1) = 0.673926;
    A_matrix(2, 2) = 2.63817;
    A_matrix(2, 3) = 1.0287;
    A_matrix(3, 0) = 0.366503;
    A_matrix(3, 1) = 3.02634;
    A_matrix(3, 2) = 1.24104;
    A_matrix(3, 3) = 3.62022;

    b_vector[0] = 0.0;
    b_vector[1] = 1.0;
    b_vector[2] = 2.0;
    b_vector[3] = 3.0;

    // Calculate the QR decomposition
    using DenseSpace = UblasSpace<double, Matrix, Vector>;
    DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
    qr_decomposition.Compute(A_matrix);
    qr_decomposition.Solve(b_vector, x_vector);

    // Check solve result
    Vector x_ref(4);
    x_ref[0] = -0.328709086647;
    x_ref[1] = -0.605571254325;
    x_ref[2] = 0.668500192182;
    x_ref[3] = 1.13901969982;
    constexpr double tolerance = 1e-10;
    KRATOS_EXPECT_VECTOR_NEAR(x_vector, x_ref, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DenseHouseholderQRSolveMatrix, KratosCoreFastSuite)
{
    // Allocate matrices
    Matrix X_matrix(4,2);
    Matrix B_matrix(4,2);
    Matrix A_matrix(4,4);

    // Set the system to be solved
    A_matrix(0, 0) = 0.0;
    A_matrix(0, 1) = 0.979749;
    A_matrix(0, 2) = 0.494393;
    A_matrix(0, 3) = 0.23073;
    A_matrix(1, 0) = 1.79224;
    A_matrix(1, 1) = 0.198842;
    A_matrix(1, 2) = 0.074485;
    A_matrix(1, 3) = 1.45717;
    A_matrix(2, 0) = 1.6039;
    A_matrix(2, 1) = 0.673926;
    A_matrix(2, 2) = 2.63817;
    A_matrix(2, 3) = 1.0287;
    A_matrix(3, 0) = 0.366503;
    A_matrix(3, 1) = 3.02634;
    A_matrix(3, 2) = 1.24104;
    A_matrix(3, 3) = 3.62022;

    B_matrix(0,0) = 0.0;
    B_matrix(1,0) = 1.0;
    B_matrix(2,0) = 2.0;
    B_matrix(3,0) = 3.0;
    B_matrix(0,1) = 0.0;
    B_matrix(1,1) = 1.0;
    B_matrix(2,1) = 2.0;
    B_matrix(3,1) = 3.0;

    // Calculate the QR decomposition
    using DenseSpace = UblasSpace<double, Matrix, Vector>;
    DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
    qr_decomposition.Compute(A_matrix);
    qr_decomposition.Solve(B_matrix, X_matrix);

    // Check solve result
    Matrix X_ref(4,2);
    X_ref(0,0) = -0.328709086647;
    X_ref(1,0) = -0.605571254325;
    X_ref(2,0) = 0.668500192182;
    X_ref(3,0) = 1.13901969982;
    X_ref(0,1) = -0.328709086647;
    X_ref(1,1) = -0.605571254325;
    X_ref(2,1) = 0.668500192182;
    X_ref(3,1) = 1.13901969982;
    constexpr double tolerance = 1e-10;
    KRATOS_EXPECT_MATRIX_NEAR(X_matrix, X_ref, tolerance);
}

} // namespace Testing
}  // namespace Kratos.
