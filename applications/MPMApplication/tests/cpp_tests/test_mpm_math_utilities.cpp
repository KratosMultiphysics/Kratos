//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "custom_utilities/mpm_math_utilities.h"

namespace Kratos
{
namespace Testing
{

    // Number of Iteration and computation tolerance
    static constexpr unsigned int num_iteration = 100;
    static constexpr double comp_tolerance = 1.0e-9;

    // Generate 3x3 test matrix
    Matrix CreateTestMatrix3x3()
    {
        Matrix matrix = ZeroMatrix(3,3);
        matrix(0,0) = -2.0;
        matrix(0,1) = -4.0;
        matrix(0,2) =  2.0;
        matrix(1,0) = -2.0;
        matrix(1,1) =  1.0;
        matrix(1,2) =  2.0;
        matrix(2,0) =  4.0;
        matrix(2,1) =  2.0;
        matrix(2,2) =  5.0;
        return matrix;
    }

    // Generate 2x2 test matrix
    Matrix CreateTestMatrix2x2()
    {
        Matrix matrix = ZeroMatrix(2,2);
        matrix(0,0) =  2.0;
        matrix(0,1) =  2.0;
        matrix(1,0) =  5.0;
        matrix(1,1) = -1.0;
        return matrix;
    }

    // Generate 3x3 test matrix
    Matrix CreateSymmetricTestMatrix3x3()
    {
        Matrix matrix = ZeroMatrix(3,3);
        matrix(0,0) =  3.0;
        matrix(1,1) =  0.0;
        matrix(2,2) =  3.0;
        matrix(0,1) =  2.0;
        matrix(0,2) =  4.0;
        matrix(1,2) =  2.0;
        matrix(1,0) =  matrix(0,1);
        matrix(2,0) =  matrix(0,2);
        matrix(2,1) =  matrix(1,2);
        return matrix;
    }

    // Generate 3x3 test matrix
    Matrix CreateSymmetricTest2Matrix3x3()
    {
        Matrix matrix = ZeroMatrix(3,3);
        matrix(0,0) =  2.0;
        matrix(1,1) =  8.0;
        matrix(2,2) =  6.0;
        matrix(0,1) =  4.0;
        matrix(0,2) =  0.0;
        matrix(1,2) =  0.0;
        matrix(1,0) =  matrix(0,1);
        matrix(2,0) =  matrix(0,2);
        matrix(2,1) =  matrix(1,2);
        return matrix;
    }

    // Generate 3x3 test matrix
    Matrix CreateSymmetricTest3Matrix3x3()
    {
        Matrix matrix = ZeroMatrix(3,3);
        matrix(0,0) =  1.0;
        matrix(1,1) =  0.0;
        matrix(2,2) =  1.0;
        matrix(0,1) =  1.0;
        matrix(0,2) =  0.0;
        matrix(1,2) =  1.0;
        matrix(1,0) =  matrix(0,1);
        matrix(2,0) =  matrix(0,2);
        matrix(2,1) =  matrix(1,2);
        return matrix;
    }

    // Generate Vector of size 3
    Vector CreateTestVector3()
    {
        Vector vector = ZeroVector(3);
        std::fill(vector.begin(), vector.end(), 1.0);
        return vector;
    }

    // Generate Vector of size 6
    Vector CreateTestVector6()
    {
        Vector vector = ZeroVector(6);
        std::fill(vector.begin(), vector.end(), 1.0);
        return vector;
    }


    /**
    * Check whether the computation of eigenvalues and eigenvectors are performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMMathUtilsQRFactorizationCalculation, KratosMPMFastSuite)
    {
        // Initialize vectors and matrices
        Matrix Q = ZeroMatrix(3,3);
        Matrix R = ZeroMatrix(3,3);

        // QR Factorization
        Matrix A = CreateSymmetricTest3Matrix3x3();
        MPMMathUtilities<double>::QRFactorization(A, Q, R);

        Matrix Q_ref(3,3);
        Q_ref(0,0) = -7.071068e-01;
        Q_ref(0,1) =  4.082483e-01;
        Q_ref(0,2) = -5.773503e-01;
        Q_ref(1,0) = -7.071068e-01;
        Q_ref(1,1) = -4.082483e-01;
        Q_ref(1,2) =  5.773503e-01;
        Q_ref(2,0) =  0.000000e+00;
        Q_ref(2,1) =  8.164966e-01;
        Q_ref(2,2) =  5.773503e-01;
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(Q, Q_ref, 1e-6);

        Matrix R_ref(3,3);
        R_ref(0,0) = -1.414214e+00;
        R_ref(0,1) = -7.071068e-01;
        R_ref(0,2) = -7.071068e-01;
        R_ref(1,0) =  0.000000e+00;
        R_ref(1,1) =  1.224745e+00;
        R_ref(1,2) =  4.082483e-01;
        R_ref(2,0) =  0.000000e+00;
        R_ref(2,1) =  0.000000e+00;
        R_ref(2,2) =  1.154701e+00;
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(R, R_ref, 1e-6);

        // Check False
        A.resize(4,3,false);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            MPMMathUtilities<double>::QRFactorization(A, Q, R),
            " GIVEN MATRIX IS NOT A SQUARE MATRIX: QRFactorization calculation");

    }


    /**
    * Check whether the computation of eigenvalues and eigenvectors are performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMMathUtilsEigenValueVectorsCalculation, KratosMPMFastSuite)
    {
        // Initialize vectors and matrices
        Vector eigen_values_2  = ZeroVector(2);
        Vector eigen_values_3  = ZeroVector(3);
        Matrix eigen_vectors_3 = ZeroMatrix(3,3);

        // 1. Compute EigenValues
        Matrix A = CreateTestMatrix3x3();
        noalias(eigen_values_3) = MPMMathUtilities<double>::EigenValues(A);

        std::vector<double> eigen_values_3_ref = {-5.0, 3.0, 6.0};
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR( eigen_values_3, eigen_values_3_ref, 1e-6);

        Matrix B = CreateTestMatrix2x2();
        noalias(eigen_values_2) = MPMMathUtilities<double>::EigenValues(B);

        std::vector<double> eigen_values_2_ref = {4.0, -3.0};
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR( eigen_values_2, eigen_values_2_ref, 1e-6);

        // 2. Compute EigenValues using direct method
        Matrix C = CreateSymmetricTestMatrix3x3();
        noalias(eigen_values_3) = MPMMathUtilities<double>::EigenValuesDirectMethod(C);

        eigen_values_3_ref = {8.0, -1.0, -1.0};
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR( eigen_values_3, eigen_values_3_ref,1e-6);

        Matrix D = CreateSymmetricTest2Matrix3x3();
        noalias(eigen_values_3) = MPMMathUtilities<double>::EigenValuesDirectMethod(D);

        eigen_values_3_ref = {10, 6.0, 0.0};
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(eigen_values_3, eigen_values_3_ref, 1e-6);

        // 3. Compute EigenVectors and EigenValues of 3x3 symmetric matrices - using Gauss Seidel method
        MPMMathUtilities<double>::EigenVectors(C, eigen_vectors_3, eigen_values_3, comp_tolerance, num_iteration);

        eigen_values_3_ref = {-1.0, -1.0, 8.0};
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(eigen_values_3, eigen_values_3_ref, 1e-6);

        Matrix eigen_vectors_3_ref(3,3);
        eigen_vectors_3_ref(0,0) =  7.071068e-01;
        eigen_vectors_3_ref(0,1) =  0.000000e+00;
        eigen_vectors_3_ref(0,2) = -7.071068e-01;
        eigen_vectors_3_ref(1,0) = -2.357023e-01;
        eigen_vectors_3_ref(1,1) =  9.428090e-01;
        eigen_vectors_3_ref(1,2) = -2.357023e-01;
        eigen_vectors_3_ref(2,0) =  6.666667e-01;
        eigen_vectors_3_ref(2,1) =  3.333333e-01;
        eigen_vectors_3_ref(2,2) =  6.666667e-01;
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR( eigen_vectors_3, eigen_vectors_3_ref, 1e-6);

        MPMMathUtilities<double>::EigenVectors(D, eigen_vectors_3, eigen_values_3, comp_tolerance, num_iteration);

        eigen_values_3_ref = {0.0, 10.0, 6.0};
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR( eigen_values_3, eigen_values_3_ref, 1e-6);

        eigen_vectors_3_ref(0,0) =  8.944272e-01;
        eigen_vectors_3_ref(0,1) = -4.472136e-01;
        eigen_vectors_3_ref(0,2) =  0.000000e+00;
        eigen_vectors_3_ref(1,0) =  4.472136e-01;
        eigen_vectors_3_ref(1,1) =  8.944272e-01;
        eigen_vectors_3_ref(1,2) =  0.000000e+00;
        eigen_vectors_3_ref(2,0) =  0.000000e+00;
        eigen_vectors_3_ref(2,1) =  0.000000e+00;
        eigen_vectors_3_ref(2,2) =  1.000000e+00;
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR( eigen_vectors_3, eigen_vectors_3_ref, 1e-6);
    }

    /**
    * Check norm computation
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMMathUtilsNormCalculation, KratosMPMFastSuite)
    {
        Vector V3 = CreateTestVector3();
        Vector V6 = CreateTestVector6();
        Matrix M  = CreateSymmetricTest2Matrix3x3();

        MPMMathUtilities<double>::Normalize(V3);
        MPMMathUtilities<double>::Normalize(V6);
        const double norm_M = MPMMathUtilities<double>::NormTensor(M);

        for (unsigned int i = 0; i<3; ++i)
            KRATOS_EXPECT_NEAR(V3[i], 0.5773502692, 1e-6);

        for (unsigned int i = 0; i<6; ++i)
            KRATOS_EXPECT_NEAR(V6[i], 0.4082482905, 1e-6);

        KRATOS_EXPECT_NEAR(norm_M, 11.66190379, 1e-6);
    }

} // namespace Testing
} // namespace Kratos
