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

#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "testing/testing.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateTrace, KratosGeoMechanicsFastSuite)
{
    Vector eye = MathUtils<double>::StressTensorToVector(IdentityMatrix(3));
    KRATOS_EXPECT_DOUBLE_EQ(3.0, StressStrainUtilities::CalculateTrace(eye));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMeanStress, KratosGeoMechanicsFastSuite)
{
    Vector eye = MathUtils<double>::StressTensorToVector(IdentityMatrix(3));
    KRATOS_EXPECT_DOUBLE_EQ(1.0, StressStrainUtilities::CalculateMeanStress(eye));
}

KRATOS_TEST_CASE_IN_SUITE(CheckGreenLagrangeStrainTensor, KratosGeoMechanicsFastSuite)
{
    Matrix eye = StressStrainUtilities::CalculateGreenLagrangeStrainTensor(sqrt(3.) * IdentityMatrix(3));
    KRATOS_EXPECT_DOUBLE_EQ(1.0, eye(0,0));
    KRATOS_EXPECT_DOUBLE_EQ(1.0, eye(1,1));
    KRATOS_EXPECT_DOUBLE_EQ(1.0, eye(2,2));
    KRATOS_EXPECT_DOUBLE_EQ(0.0, eye(0,1));
    KRATOS_EXPECT_DOUBLE_EQ(0.0, eye(0,2));
    KRATOS_EXPECT_DOUBLE_EQ(0.0, eye(1,2));
}


KRATOS_TEST_CASE_IN_SUITE(CheckCalculateVonMisesStress, KratosGeoMechanicsFastSuite)
{
    Vector stress_vector(4);
    stress_vector <<= -2.0, -2.0, -2.0, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateVonMisesStress(stress_vector));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateLodeAngle, KratosGeoMechanicsFastSuite)
{
    Vector stress_vector(4);
    stress_vector <<= 0.5, -1.0, 0.5, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(30.), StressStrainUtilities::CalculateLodeAngle(stress_vector));
    KRATOS_EXPECT_DOUBLE_EQ(-MathUtils<>::DegreesToRadians(30.), StressStrainUtilities::CalculateLodeAngle(-1.*stress_vector));
    Vector stress_vector2(4);
    stress_vector <<= -1.0, 0.0, 1.0, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(0.), StressStrainUtilities::CalculateLodeAngle(stress_vector2));
}

}