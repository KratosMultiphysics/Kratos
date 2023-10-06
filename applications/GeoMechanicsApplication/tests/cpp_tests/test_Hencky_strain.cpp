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

#include "testing/testing.h"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.hpp"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CheckHenckyStrainIncompressibleUniaxialElongation, KratosGeoMechanicsFastSuite)
{
    double stretch = 2.0;
    Matrix Fmatrix = ZeroMatrix(3,3);
    Fmatrix(0,0) = stretch;
    Fmatrix(1,1) = 1.0/std::sqrt(stretch);
    Fmatrix(2,2) = Fmatrix(1,1);

    Vector Evector = StressStrainUtilities::CalculateHenckyStrain(Fmatrix, 6);
    KRATOS_EXPECT_DOUBLE_EQ(0.5*std::log(stretch*stretch), Evector[0]);
    KRATOS_EXPECT_DOUBLE_EQ(0.5*std::log(1./stretch), Evector[1]);
    KRATOS_EXPECT_DOUBLE_EQ(0.5*std::log(1./stretch), Evector[2]);
    KRATOS_EXPECT_DOUBLE_EQ(0., Evector[3]);
}

KRATOS_TEST_CASE_IN_SUITE(CheckHenckyStrainRigidRotation, KratosGeoMechanicsFastSuite)
{
    double RotationAngle = 1.;
    Matrix Fmatrix = ZeroMatrix(2,2);
    Fmatrix(0,0) =  std::cos(RotationAngle); Fmatrix(0,1) =  std::sin(RotationAngle);
    Fmatrix(1,0) = -std::sin(RotationAngle); Fmatrix(1,1) =  std::cos(RotationAngle);

    // Plane strain has VoigtSize 4, rigid rotation results in all zero strain components
    Vector Evector = StressStrainUtilities::CalculateHenckyStrain(Fmatrix, 4);
    KRATOS_EXPECT_DOUBLE_EQ(0., Evector[0]);
    KRATOS_EXPECT_DOUBLE_EQ(0., Evector[1]);
    KRATOS_EXPECT_DOUBLE_EQ(0., Evector[2]);
    KRATOS_EXPECT_DOUBLE_EQ(0., Evector[3]);
}

KRATOS_TEST_CASE_IN_SUITE(CheckHenckyStrainPureShear, KratosGeoMechanicsFastSuite)
{
    double gamma = 0.5;
    Matrix Fmatrix = IdentityMatrix (2,2);
    Fmatrix(0,1) = gamma;
    Fmatrix(1,0) = gamma;
    Vector Evector = StressStrainUtilities::CalculateHenckyStrain(Fmatrix, 4);

    // analytical eigenvalues of right Cauchy-Green deformation tensor
    // ( such that we can do the natural logarithm for the Hencky strain on a diagonal )
    Matrix C = prod(trans(Fmatrix), Fmatrix);
    double T = C(0,0) + C(1,1);
    double D = C(0,0)*C(1,1) - C(0,1)*C(1,0);
    double L1 = 0.5*T + std::sqrt(0.25*T*T - D);
    double L2 = 0.5*T - std::sqrt(0.25*T*T - D);
    // normalized eigen vectors
    Vector V1(2);
    V1(0) = L1 - C(1,1);
    V1(1) = C(1,0);
    V1 /= norm_2(V1);
    Vector V2(2);
    V2(0) = L2- C(1,1);
    V2(1) = C(1,0);
    V2 /= norm_2(V2);
    // natural logarithm on diagonal followed by triple product to rotate back to original coordinate system
    Matrix Eprin = zero_matrix(2,2);
    Eprin(0,0) = 0.5*std::log(L1);
    Eprin(1,1) = 0.5*std::log(L2);
    Matrix V12(2,2);
    V12(0,0) = V1(0);
    V12(1,0) = V1(1);
    V12(0,1) = V2(0);
    V12(1,1) = V2(1);

    Matrix Eaux = prod(V12,Eprin);
    Matrix Eana = prod( Eaux, trans(V12));

    KRATOS_EXPECT_DOUBLE_EQ(Eana(0,0),Evector(0));
    KRATOS_EXPECT_DOUBLE_EQ(Eana(1,1),Evector(1));
    KRATOS_EXPECT_DOUBLE_EQ(0.,Evector(2));
    KRATOS_EXPECT_DOUBLE_EQ(2.*Eana(1,0),Evector(3));
}

}