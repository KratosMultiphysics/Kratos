// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CheckHenckyStrainIncompressibleUniaxialElongation, KratosGeoMechanicsFastSuite)
{
    const double stretch = 2.0;
    Matrix Fmatrix = ZeroMatrix(3,3);
    Fmatrix(0,0) = stretch;
    Fmatrix(1,1) = 1.0/std::sqrt(stretch);
    Fmatrix(2,2) = Fmatrix(1,1);

    const Vector Evector = StressStrainUtilities::CalculateHenckyStrain(Fmatrix, 6);
    // E = 0.5 ln(C) = 0.5 ln(F^T F)
    KRATOS_EXPECT_DOUBLE_EQ(0.5*std::log(stretch*stretch), Evector[0]);
    KRATOS_EXPECT_DOUBLE_EQ(0.5*std::log(1./stretch), Evector[1]);
    KRATOS_EXPECT_DOUBLE_EQ(0.5*std::log(1./stretch), Evector[2]);
    KRATOS_EXPECT_DOUBLE_EQ(0.0, Evector[3]);
    KRATOS_EXPECT_DOUBLE_EQ(0.0, Evector[4]);
    KRATOS_EXPECT_DOUBLE_EQ(0.0, Evector[5]);
}

KRATOS_TEST_CASE_IN_SUITE(CheckHenckyStrainRigidRotation, KratosGeoMechanicsFastSuite)
{
    const double RotationAngle = 1.0;
    Matrix Fmatrix(2,2);
    Fmatrix(0,0) =  std::cos(RotationAngle); Fmatrix(0,1) =  std::sin(RotationAngle);
    Fmatrix(1,0) = -std::sin(RotationAngle); Fmatrix(1,1) =  std::cos(RotationAngle);

    // Plane strain has VoigtSize 4, rigid rotation results in all zero strain components
    const Vector Evector = StressStrainUtilities::CalculateHenckyStrain(Fmatrix, 4);
    KRATOS_EXPECT_DOUBLE_EQ(0.0, Evector[0]);
    KRATOS_EXPECT_DOUBLE_EQ(0.0, Evector[1]);
    KRATOS_EXPECT_DOUBLE_EQ(0.0, Evector[2]);
    KRATOS_EXPECT_DOUBLE_EQ(0.0, Evector[3]);
}

KRATOS_TEST_CASE_IN_SUITE(CheckHenckyStrainPureShear, KratosGeoMechanicsFastSuite)
{
    const double gamma = 0.5;
    Matrix Fmatrix = IdentityMatrix (2,2);
    Fmatrix(0,1) = gamma;
    Fmatrix(1,0) = gamma;
    const Vector Evector = StressStrainUtilities::CalculateHenckyStrain(Fmatrix, 4);

    // analytical eigenvalues of right Cauchy-Green deformation tensor
    // ( such that we can do the natural logarithm for the Hencky strain on a diagonal )
    const Matrix C = prod(trans(Fmatrix), Fmatrix); // Cauchy-Green
    const double T = C(0,0) + C(1,1); // Trace of C
    const double D = C(0,0)*C(1,1) - C(0,1)*C(1,0); // Determinant of C
    const double L1 = 0.5*T + std::sqrt(0.25*T*T - D); // First eigenvalue
    const double L2 = 0.5*T - std::sqrt(0.25*T*T - D); // Second eigenvalue
    // normalized eigen vectors
    Vector eigV1(2);
    eigV1(0) = L1 - C(1,1);
    eigV1(1) = C(1,0);
    eigV1 /= norm_2(eigV1);
    Vector eigV2(2);
    eigV2(0) = L2- C(1,1);
    eigV2(1) = C(1,0);
    eigV2 /= norm_2(eigV2);
    // natural logarithm on diagonal followed by triple product to rotate back to original coordinate system
    Matrix EHenckyPrincipal = zero_matrix(2,2);
    EHenckyPrincipal(0,0) = 0.5*std::log(L1);
    EHenckyPrincipal(1,1) = 0.5*std::log(L2);
    Matrix eigV12(2,2);
    eigV12(0,0) = eigV1(0);
    eigV12(1,0) = eigV1(1);
    eigV12(0,1) = eigV2(0);
    eigV12(1,1) = eigV2(1);

    const Matrix Eaux = prod(eigV12,EHenckyPrincipal);
    const Matrix EHenckyAnalytical = prod( Eaux, trans(eigV12));

    KRATOS_EXPECT_DOUBLE_EQ(EHenckyAnalytical(0,0), Evector(0));
    KRATOS_EXPECT_DOUBLE_EQ(EHenckyAnalytical(1,1), Evector(1));
    KRATOS_EXPECT_DOUBLE_EQ(0.0, Evector(2));
    KRATOS_EXPECT_DOUBLE_EQ(2.0*EHenckyAnalytical(1,0), Evector(3));
}

}