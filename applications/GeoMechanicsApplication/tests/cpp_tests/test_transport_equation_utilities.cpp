// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_utilities/transport_equation_utilities.hpp"
#include "testing/testing.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateBiotModulusInverse_GivesExpectedResult, KratosGeoMechanicsFastSuite)
{
    Properties properties;
    properties[IGNORE_UNDRAINED]          = false;
    properties[POROSITY]                  = 0.5;
    properties[BULK_MODULUS_SOLID]        = 1.0e9;
    properties[BULK_MODULUS_FLUID]        = 2.0e6;
    const double biot_coefficient         = 1.0;
    const double degree_of_saturation     = 0.3;
    const double derivative_of_saturation = 0.2;

    const double expected_value = -0.09999992485;
    KRATOS_EXPECT_DOUBLE_EQ(GeoTransportEquationUtilities::CalculateBiotModulusInverse(
                                biot_coefficient, degree_of_saturation, derivative_of_saturation, properties),
                            expected_value);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateBiotModulusInverse_ReturnsLargeNumber_WhenIgnoreUndrained, KratosGeoMechanicsFastSuite)
{
    Properties properties;
    properties[IGNORE_UNDRAINED]          = true;
    properties[POROSITY]                  = 0.5;
    properties[BULK_MODULUS_SOLID]        = 1.0e9;
    properties[BULK_MODULUS_FLUID]        = 2.0e6;
    const double biot_coefficient         = 1.0;
    const double degree_of_saturation     = 0.3;
    const double derivative_of_saturation = 0.2;

    const auto large_number = 1e10;
    KRATOS_EXPECT_TRUE(GeoTransportEquationUtilities::CalculateBiotModulusInverse(
                           biot_coefficient, degree_of_saturation, derivative_of_saturation, properties) > large_number)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateBiotModulusInverse_DoesNotThrow_ForEmptyProperties, KratosGeoMechanicsFastSuite)
{
    Properties   properties;
    const double biot_coefficient         = 1.0;
    const double degree_of_saturation     = 0.3;
    const double derivative_of_saturation = 0.2;

    KRATOS_EXPECT_TRUE(std::isnan(GeoTransportEquationUtilities::CalculateBiotModulusInverse(
        biot_coefficient, degree_of_saturation, derivative_of_saturation, properties)))
}

} // namespace Kratos::Testing