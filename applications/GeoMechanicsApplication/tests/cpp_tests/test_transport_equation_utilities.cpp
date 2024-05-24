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
#include <boost/numeric/ublas/assignment.hpp>

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
    KRATOS_EXPECT_DOUBLE_EQ(GeoTransportEquationUtilities::CalculateInverseBiotModulus(
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
    KRATOS_EXPECT_GT(GeoTransportEquationUtilities::CalculateInverseBiotModulus(
                         biot_coefficient, degree_of_saturation, derivative_of_saturation, properties),
                     large_number);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateBiotModulusInverse_DoesNotThrow_ForEmptyProperties, KratosGeoMechanicsFastSuite)
{
    const Properties properties;
    const double     biot_coefficient         = 1.0;
    const double     degree_of_saturation     = 0.3;
    const double     derivative_of_saturation = 0.2;

    KRATOS_EXPECT_TRUE(std::isnan(GeoTransportEquationUtilities::CalculateInverseBiotModulus(
        biot_coefficient, degree_of_saturation, derivative_of_saturation, properties)))
}

KRATOS_TEST_CASE_IN_SUITE(CalculateBulkModulus_ReturnsZeroForZeroConstitutiveMatrix, KratosGeoMechanicsFastSuite)
{
    ZeroMatrix constitutive_matrix(3, 3);
    KRATOS_EXPECT_DOUBLE_EQ(GeoTransportEquationUtilities::CalculateBulkModulus(constitutive_matrix), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateBulkModulus_Throws_WhenConstitutiveMatrixIsEmpty, KratosGeoMechanicsFastSuite)
{
    const Matrix constitutive_matrix;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto bulk_modulus =
            GeoTransportEquationUtilities::CalculateBulkModulus(constitutive_matrix),
        "Constitutive matrix is empty, aborting bulk modulus calculation.")
}

KRATOS_TEST_CASE_IN_SUITE(CalculateBulkModulus_GivesExpectedResult_ForFilledConstitutiveMatrix, KratosGeoMechanicsFastSuite)
{
    Matrix constitutive_matrix(3, 3);
    constitutive_matrix <<= 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;

    KRATOS_EXPECT_DOUBLE_EQ(GeoTransportEquationUtilities::CalculateBulkModulus(constitutive_matrix), -11.0);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateBiotCoefficients_GivesExpectedResults_WhenPropertyHasBiotCoefficient,
                          KratosGeoMechanicsFastSuite)
{
    const std::vector<Matrix> constitutive_matrices(2, ZeroMatrix(3, 3));

    const double biot_coefficient = 1.2;
    Properties   properties;
    properties[BIOT_COEFFICIENT] = biot_coefficient;

    const auto actual_values =
        GeoTransportEquationUtilities::CalculateBiotCoefficients(constitutive_matrices, properties);

    const std::vector<double> expected_values(2, biot_coefficient);
    KRATOS_CHECK_VECTOR_NEAR(expected_values, actual_values, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateBiotCoefficients_GivesExpectedResults, KratosGeoMechanicsFastSuite)
{
    std::vector<Matrix> constitutive_matrices;
    constitutive_matrices.emplace_back(ScalarMatrix(3, 3, 1.0));
    constitutive_matrices.emplace_back(ScalarMatrix(3, 3, 3.0));

    Properties properties;
    properties[BULK_MODULUS_SOLID] = 1.0;

    const auto actual_values =
        GeoTransportEquationUtilities::CalculateBiotCoefficients(constitutive_matrices, properties);

    const std::vector<double> expected_values = {4.0 / 3.0, 2.0};
    KRATOS_CHECK_VECTOR_NEAR(expected_values, actual_values, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateBiotCoefficients_GivesInfResults_ForZeroBulkModulus, KratosGeoMechanicsFastSuite)
{
    std::vector<Matrix> constitutive_matrices;
    constitutive_matrices.emplace_back(ScalarMatrix(3, 3, 1.0));
    constitutive_matrices.emplace_back(ScalarMatrix(3, 3, 3.0));

    Properties properties;
    properties[BULK_MODULUS_SOLID] = 0.0;

    const auto actual_values =
        GeoTransportEquationUtilities::CalculateBiotCoefficients(constitutive_matrices, Properties());

    // Due to a division by zero, we end up with inf
    KRATOS_EXPECT_TRUE(std::all_of(actual_values.begin(), actual_values.end(),
                                   [](const double value) { return std::isinf(value); }))
}

} // namespace Kratos::Testing