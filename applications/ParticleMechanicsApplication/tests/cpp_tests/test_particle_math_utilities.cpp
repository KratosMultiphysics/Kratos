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
#include "utilities/math_utils.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"

namespace Kratos
{
namespace Testing
{

    // Tolerance
    static constexpr double tolerance = 1.0e-6;

    // Number of Iteration and computation tolerance
    static constexpr unsigned int num_iteration = 100;
    static constexpr double comp_tolerance = 1.0e-9;

    // Generate 3x3 test matrix
    Matrix CreateTestMatrix3x3()
    {
        Matrix matrix = ZeroMatrix(3);
        matrix(0,0) = 3.0;
        matrix(0,1) = 0.0;
        matrix(0,2) = 2.0;
        matrix(1,0) = 2.0;
        matrix(1,1) = 0.0;
        matrix(1,2) =-2.0;
        matrix(2,0) = 0.0;
        matrix(2,1) = 1.0;
        matrix(2,2) = 1.0;
        return matrix;
    }

    // Generate 3x3 test matrix
    Matrix CreateTest2Matrix3x3()
    {
        Matrix matrix = ZeroMatrix(3);
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
        Matrix matrix = ZeroMatrix(2);
        matrix(0,0) =  2.0;
        matrix(0,1) =  2.0;
        matrix(1,0) =  5.0;
        matrix(1,1) = -1.0;
        return matrix;
    }

    // Generate 2x2 singular matrix
    Matrix CreateTestSingularMatrix2x2()
    {
        Matrix matrix = ZeroMatrix(2);
        matrix(0,0) =  4.0;
        matrix(0,1) =  2.0;
        matrix(1,0) =  2.0;
        matrix(1,1) =  1.0;
        return matrix;
    }

    // Generate 3x3 test matrix
    Matrix CreateSymmetricTestMatrix3x3()
    {
        Matrix matrix = ZeroMatrix(3);
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
        Matrix matrix = ZeroMatrix(3);
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
        Matrix matrix = ZeroMatrix(3);
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


    /**
    * Check whether the computation of eigenvalues and eigenvectors are performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleMathUtilsQRFactorizationCalculation, KratosParticleMechanicsFastSuite)
    {
        // Initialize vectors and matrices
        Matrix Q = ZeroMatrix(3,3);
        Matrix R = ZeroMatrix(3,3);

        // QR Factorization
        Matrix A = CreateSymmetricTest3Matrix3x3();
        ParticleMechanicsMathUtilities<double>::QRFactorization(A, Q, R);

        KRATOS_CHECK_LESS_EQUAL((-7.071068e-01 - Q(0,0))/Q(0,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 4.082483e-01 - Q(0,1))/Q(0,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL((-5.773503e-01 - Q(0,2))/Q(0,2), tolerance);
        KRATOS_CHECK_LESS_EQUAL((-7.071068e-01 - Q(1,0))/Q(1,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL((-4.082483e-01 - Q(1,1))/Q(1,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 5.773503e-01 - Q(1,2))/Q(1,2), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.000000e+00 - Q(2,0))/Q(2,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 8.164966e-01 - Q(2,1))/Q(2,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 5.773503e-01 - Q(2,2))/Q(2,2), tolerance);

        KRATOS_CHECK_LESS_EQUAL((-1.414214e+00 - R(0,0))/R(0,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL((-7.071068e-01 - R(0,1))/R(0,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL((-7.071068e-01 - R(0,2))/R(0,2), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.000000e+00 - R(1,0))/R(1,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 1.224745e+00 - R(1,1))/R(1,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 4.082483e-01 - R(1,2))/R(1,2), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.000000e+00 - R(2,0))/R(2,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.000000e+00 - R(2,1))/R(2,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 1.154701e+00 - R(2,2))/R(2,2), tolerance);

        // Check False
        A.resize(4,3,false);
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            ParticleMechanicsMathUtilities<double>::QRFactorization(A, Q, R),
            " GIVEN MATRIX IS NOT A SQUARE MATRIX: QRFactorization calculation");

    }


    /**
    * Check whether the computation of eigenvalues and eigenvectors are performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleMathUtilsEigenValueVectorsCalculation, KratosParticleMechanicsFastSuite)
    {
        // Initialize vectors and matrices
        Vector eigen_values_2  = ZeroVector(2);
        Vector eigen_values_3  = ZeroVector(3);
        Matrix eigen_vectors_3 = ZeroMatrix(3,3);

        // 1. Compute EigenValues
        Matrix A = CreateTest2Matrix3x3();
        noalias(eigen_values_3) = ParticleMechanicsMathUtilities<double>::EigenValues(A);

        KRATOS_CHECK_LESS_EQUAL(( -5.0 - eigen_values_3[0])/eigen_values_3[0], tolerance);
        KRATOS_CHECK_LESS_EQUAL((  3.0 - eigen_values_3[1])/eigen_values_3[1], tolerance);
        KRATOS_CHECK_LESS_EQUAL((  6.0 - eigen_values_3[2])/eigen_values_3[2], tolerance);

        Matrix B = CreateTestMatrix2x2();
        noalias(eigen_values_2) = ParticleMechanicsMathUtilities<double>::EigenValues(B);

        KRATOS_CHECK_LESS_EQUAL(( 4.0 - eigen_values_2[0])/eigen_values_2[0], tolerance);
        KRATOS_CHECK_LESS_EQUAL((-3.0 - eigen_values_2[1])/eigen_values_2[1], tolerance);

        // 2. Compute EigenValues using direct method
        Matrix C = CreateSymmetricTestMatrix3x3();
        noalias(eigen_values_3) = ParticleMechanicsMathUtilities<double>::EigenValuesDirectMethod(C);

        KRATOS_CHECK_LESS_EQUAL(( 8.0 - eigen_values_3[0])/eigen_values_3[0], tolerance);
        KRATOS_CHECK_LESS_EQUAL((-1.0 - eigen_values_3[1])/eigen_values_3[1], tolerance);
        KRATOS_CHECK_LESS_EQUAL((-1.0 - eigen_values_3[2])/eigen_values_3[2], tolerance);

        Matrix D = CreateSymmetricTest2Matrix3x3();
        noalias(eigen_values_3) = ParticleMechanicsMathUtilities<double>::EigenValuesDirectMethod(D);

        KRATOS_CHECK_LESS_EQUAL((10.0 - eigen_values_3[0])/eigen_values_3[0], tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 6.0 - eigen_values_3[1])/eigen_values_3[1], tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.0 - eigen_values_3[2])/eigen_values_3[2], tolerance);

        // 3. Compute EigenVectors and EigenValues of 3x3 symmetric matrices - using Gauss Seidel method
        ParticleMechanicsMathUtilities<double>::EigenVectors(C, eigen_vectors_3, eigen_values_3, comp_tolerance, num_iteration);

        KRATOS_CHECK_LESS_EQUAL((-1.0 - eigen_values_3[0])/eigen_values_3[0], tolerance);
        KRATOS_CHECK_LESS_EQUAL((-1.0 - eigen_values_3[1])/eigen_values_3[1], tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 8.0 - eigen_values_3[2])/eigen_values_3[2], tolerance);

        KRATOS_CHECK_LESS_EQUAL(( 7.071068e-01 - eigen_vectors_3(0,0))/eigen_vectors_3(0,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.000000e+00 - eigen_vectors_3(0,1))/eigen_vectors_3(0,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL((-7.071068e-01 - eigen_vectors_3(0,2))/eigen_vectors_3(0,2), tolerance);
        KRATOS_CHECK_LESS_EQUAL((-2.357023e-01 - eigen_vectors_3(1,0))/eigen_vectors_3(1,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 9.428090e-01 - eigen_vectors_3(1,1))/eigen_vectors_3(1,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL((-2.357023e-01 - eigen_vectors_3(1,2))/eigen_vectors_3(1,2), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 6.666667e-01 - eigen_vectors_3(2,0))/eigen_vectors_3(2,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 3.333333e-01 - eigen_vectors_3(2,1))/eigen_vectors_3(2,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 6.666667e-01 - eigen_vectors_3(2,2))/eigen_vectors_3(2,2), tolerance);

        ParticleMechanicsMathUtilities<double>::EigenVectors(D, eigen_vectors_3, eigen_values_3, comp_tolerance, num_iteration);

        KRATOS_CHECK_LESS_EQUAL(( 0.0 - eigen_values_3[0])/eigen_values_3[0], tolerance);
        KRATOS_CHECK_LESS_EQUAL((10.0 - eigen_values_3[1])/eigen_values_3[1], tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 6.0 - eigen_values_3[2])/eigen_values_3[2], tolerance);

        KRATOS_CHECK_LESS_EQUAL(( 8.944272e-01 - eigen_vectors_3(0,0))/eigen_vectors_3(0,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL((-4.472136e-01 - eigen_vectors_3(0,1))/eigen_vectors_3(0,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.000000e+00 - eigen_vectors_3(0,2))/eigen_vectors_3(0,2), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 4.472136e-01 - eigen_vectors_3(1,0))/eigen_vectors_3(1,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 8.944272e-01 - eigen_vectors_3(1,1))/eigen_vectors_3(1,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.000000e+00 - eigen_vectors_3(1,2))/eigen_vectors_3(1,2), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.000000e+00 - eigen_vectors_3(2,0))/eigen_vectors_3(2,0), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.000000e+00 - eigen_vectors_3(2,1))/eigen_vectors_3(2,1), tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 1.000000e+00 - eigen_vectors_3(2,2))/eigen_vectors_3(2,2), tolerance);

    }

    /**
    * Check inverse computation
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleMathUtilsInverseCalculation, KratosParticleMechanicsFastSuite)
    {
        // 1. Check 3x3 inverse
        Matrix A = CreateTestMatrix3x3();
        Matrix inv_A = ZeroMatrix(3);

        ParticleMechanicsMathUtilities<double>::InvertMatrix( A, inv_A);

        KRATOS_CHECK_NEAR(inv_A(0,0), 0.2,tolerance);
        KRATOS_CHECK_NEAR(inv_A(0,1), 0.2,tolerance);
        KRATOS_CHECK_NEAR(inv_A(0,2), 0.0,tolerance);
        KRATOS_CHECK_NEAR(inv_A(1,0),-0.2,tolerance);
        KRATOS_CHECK_NEAR(inv_A(1,1), 0.3,tolerance);
        KRATOS_CHECK_NEAR(inv_A(1,2), 1.0,tolerance);
        KRATOS_CHECK_NEAR(inv_A(2,0), 0.2,tolerance);
        KRATOS_CHECK_NEAR(inv_A(2,1),-0.3,tolerance);
        KRATOS_CHECK_NEAR(inv_A(2,2), 0.0,tolerance);

        // 2. Check 2x2 inverse
        Matrix B = CreateTestMatrix2x2();
        Matrix inv_B = ZeroMatrix(2);

        ParticleMechanicsMathUtilities<double>::InvertMatrix( B, inv_B);

        KRATOS_CHECK_NEAR(inv_B(0,0), 0.083333,tolerance);
        KRATOS_CHECK_NEAR(inv_B(0,1), 0.166667,tolerance);
        KRATOS_CHECK_NEAR(inv_B(1,0), 0.416667,tolerance);
        KRATOS_CHECK_NEAR(inv_B(1,1),-0.166667,tolerance);

        // 3. Check this will throw singular
        Matrix C = CreateTestSingularMatrix2x2();
        Matrix inv_C = ZeroMatrix(2);

        KRATOS_CHECK_NOT_EQUAL(ParticleMechanicsMathUtilities<double>::InvertMatrix(C, inv_C), 0);
    }



} // namespace Testing
} // namespace Kratos
