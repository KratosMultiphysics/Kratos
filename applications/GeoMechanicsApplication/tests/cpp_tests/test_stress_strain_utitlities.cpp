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

#include "custom_utilities/math_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "testing/testing.h"
#include "utilities/math_utils.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateTrace, KratosGeoMechanicsFastSuite)
{
    const Vector eye = MathUtils<double>::StressTensorToVector(IdentityMatrix(3));
    KRATOS_EXPECT_DOUBLE_EQ(3.0, StressStrainUtilities::CalculateTrace(eye));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMeanStress, KratosGeoMechanicsFastSuite)
{
    const Vector eye = MathUtils<double>::StressTensorToVector(IdentityMatrix(3));
    KRATOS_EXPECT_DOUBLE_EQ(1.0, StressStrainUtilities::CalculateMeanStress(eye));
}

KRATOS_TEST_CASE_IN_SUITE(CheckGreenLagrangeStrainTensor, KratosGeoMechanicsFastSuite)
{
    const Matrix eye =
        StressStrainUtilities::CalculateGreenLagrangeStrainTensor(std::sqrt(3.) * IdentityMatrix(3));
    KRATOS_EXPECT_MATRIX_NEAR(eye, IdentityMatrix(3), 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateVonMisesStressHydrostatic, KratosGeoMechanicsFastSuite)
{
    Vector stress_vector(4);
    stress_vector <<= -2.0, -2.0, -2.0, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateVonMisesStress(stress_vector));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateVonMisesStressPureShear, KratosGeoMechanicsFastSuite)
{
    Vector stress_vector(4);
    stress_vector <<= 0.0, 0.0, 0.0, 2.0;
    KRATOS_EXPECT_DOUBLE_EQ(2. * std::sqrt(3.), StressStrainUtilities::CalculateVonMisesStress(stress_vector));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateLodeAngle, KratosGeoMechanicsFastSuite)
{
    Vector stress_vector(4);
    stress_vector <<= 0.5, -1.0, 0.5, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(30.),
                            StressStrainUtilities::CalculateLodeAngle(stress_vector));
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(-30.),
                            StressStrainUtilities::CalculateLodeAngle(-1. * stress_vector));
    Vector stress_vector2(4);
    stress_vector2 <<= -1.0, 0.0, 1.0, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(0.),
                            StressStrainUtilities::CalculateLodeAngle(stress_vector2));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMohrCoulombShearCapacityZeroStress, KratosGeoMechanicsFastSuite)
{
    const Vector stress_vector  = ZeroVector(4);
    const double cohesion       = 2.;
    const double friction_angle = 0.;
    KRATOS_EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateMohrCoulombShearCapacity(
                                     stress_vector, cohesion, friction_angle));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMohrCoulombShearCapacityHydrostatic, KratosGeoMechanicsFastSuite)
{
    Vector stress_vector(4);
    stress_vector <<= -2.0, -2.0, -2.0, 0.0;
    const double cohesion       = 0.;
    const double friction_angle = MathUtils<>::DegreesToRadians(90.);
    KRATOS_EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateMohrCoulombShearCapacity(
                                     stress_vector, cohesion, friction_angle));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMohrCoulombShearCapacityShearOnly, KratosGeoMechanicsFastSuite)
{
    Vector stress_vector(4);
    stress_vector <<= -2.0, 0.0, 2.0, 0.0;
    const double cohesion       = 2.;
    const double friction_angle = 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(1.0, StressStrainUtilities::CalculateMohrCoulombShearCapacity(
                                     stress_vector, cohesion, friction_angle));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMohrCoulombPressureCapacityZeroStress, KratosGeoMechanicsFastSuite)
{
    Vector       stress_vector  = ZeroVector(4);
    const double cohesion       = 2.;
    const double friction_angle = 0.;
    KRATOS_EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateMohrCoulombPressureCapacity(
                                     stress_vector, cohesion, friction_angle));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMohrCoulombPressureCapacityShearOnly, KratosGeoMechanicsFastSuite)
{
    Vector stress_vector(4);
    stress_vector <<= -2.0, 0.0, 2.0, 0.0;
    const double cohesion       = std::sqrt(2. * 2. * 3.);
    const double friction_angle = MathUtils<>::DegreesToRadians(30.0);
    KRATOS_EXPECT_NEAR(3. * std::sin(friction_angle),
                       StressStrainUtilities::CalculateMohrCoulombPressureCapacity(
                           stress_vector, cohesion, friction_angle),
                       1.E-10);
}

} // namespace Kratos::Testing