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
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include <boost/numeric/ublas/assignment.hpp>

namespace
{

constexpr auto tolerance = 1.0e-12;

}

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateBiotModulusInverse_GivesExpectedResult)
{
    Properties properties;
    properties[IGNORE_UNDRAINED]                       = false;
    properties[POROSITY]                               = 0.5;
    properties[BULK_MODULUS_SOLID]                     = 1.0e9;
    properties[BULK_MODULUS_FLUID]                     = 2.0e6;
    const std::vector<double> biot_coefficient         = {1.0};
    const std::vector<double> degree_of_saturation     = {0.3};
    const std::vector<double> derivative_of_saturation = {0.2};

    const double expected_value = -0.09999992485;
    KRATOS_EXPECT_DOUBLE_EQ(
        GeoTransportEquationUtilities::CalculateInverseBiotModuli(
            biot_coefficient, degree_of_saturation, derivative_of_saturation, properties)[0],
        expected_value);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateBiotModulusInverse_ReturnsLargeNumber_WhenIgnoreUndrained)
{
    Properties properties;
    properties[IGNORE_UNDRAINED]                       = true;
    properties[POROSITY]                               = 0.5;
    properties[BULK_MODULUS_SOLID]                     = 1.0e9;
    properties[BULK_MODULUS_FLUID]                     = 2.0e6;
    const std::vector<double> biot_coefficient         = {1.0};
    const std::vector<double> degree_of_saturation     = {0.3};
    const std::vector<double> derivative_of_saturation = {0.2};

    const auto large_number = 1e10;
    KRATOS_EXPECT_GT(GeoTransportEquationUtilities::CalculateInverseBiotModuli(
                         biot_coefficient, degree_of_saturation, derivative_of_saturation, properties)[0],
                     large_number);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateBiotModulusInverse_DoesNotThrow_ForEmptyProperties)
{
    const Properties          properties;
    const std::vector<double> biot_coefficient         = {1.0};
    const std::vector<double> degree_of_saturation     = {0.3};
    const std::vector<double> derivative_of_saturation = {0.2};

    KRATOS_EXPECT_TRUE(std::isnan(GeoTransportEquationUtilities::CalculateInverseBiotModuli(
        biot_coefficient, degree_of_saturation, derivative_of_saturation, properties)[0]))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateBulkModulus_ReturnsZeroForZeroConstitutiveMatrix)
{
    ZeroMatrix constitutive_matrix(3, 3);
    KRATOS_EXPECT_DOUBLE_EQ(GeoTransportEquationUtilities::CalculateBulkModulus(constitutive_matrix), 0.0);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateBulkModulus_Throws_WhenConstitutiveMatrixIsEmpty)
{
    const Matrix constitutive_matrix;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto bulk_modulus =
            GeoTransportEquationUtilities::CalculateBulkModulus(constitutive_matrix),
        "Constitutive matrix is empty, aborting bulk modulus calculation.")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateBulkModulus_GivesExpectedResult_ForFilledConstitutiveMatrix)
{
    Matrix constitutive_matrix(3, 3);
    constitutive_matrix <<= 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;

    KRATOS_EXPECT_DOUBLE_EQ(GeoTransportEquationUtilities::CalculateBulkModulus(constitutive_matrix), -11.0);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateBiotCoefficients_GivesExpectedResults_WhenPropertyHasBiotCoefficient)
{
    const std::vector<Matrix> constitutive_matrices(2, ZeroMatrix(3, 3));

    const double biot_coefficient = 1.2;
    Properties   properties;
    properties[BIOT_COEFFICIENT] = biot_coefficient;

    const auto actual_values =
        GeoTransportEquationUtilities::CalculateBiotCoefficients(constitutive_matrices, properties);

    const std::vector<double> expected_values(2, biot_coefficient);
    KRATOS_EXPECT_VECTOR_NEAR(expected_values, actual_values, 1e-12)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateBiotCoefficients_GivesExpectedResults)
{
    std::vector<Matrix> constitutive_matrices;
    constitutive_matrices.emplace_back(ScalarMatrix(3, 3, 1.0));
    constitutive_matrices.emplace_back(ScalarMatrix(3, 3, 3.0));

    Properties properties;
    properties[BULK_MODULUS_SOLID] = 1.0;

    const auto actual_values =
        GeoTransportEquationUtilities::CalculateBiotCoefficients(constitutive_matrices, properties);

    const std::vector<double> expected_values = {4.0 / 3.0, 2.0};
    KRATOS_EXPECT_VECTOR_NEAR(expected_values, actual_values, 1e-12)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateBiotCoefficients_GivesInfResults_ForZeroBulkModulus)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, EachFluidPressureIsTheInnerProductOfShapeFunctionsAndPressure)
{
    auto shape_function_values = Matrix{3, 3, 0.0};
    // clang-format off
    shape_function_values <<= 0.8, 0.2, 0.3,
                              0.1, 0.7, 0.4,
                              0.2, 0.5, 0.6;
    // clang-format on

    auto pore_water_pressures = Vector(3);
    pore_water_pressures <<= 2.0, 3.0, 4.0;

    const auto expected_fluid_pressures = std::vector<double>{3.4, 3.9, 4.3};
    KRATOS_EXPECT_VECTOR_NEAR(GeoTransportEquationUtilities::CalculateFluidPressures(
                                  shape_function_values, pore_water_pressures),
                              expected_fluid_pressures, tolerance)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PermeabilityUpdateFactorEqualsOneWhenChangeInverseFactorIsNotGiven)
{
    const auto unused_strain_vectors = std::vector<Vector>(3, Vector{});
    auto       properties            = Properties{};

    const auto expected_factors = std::vector<double>(unused_strain_vectors.size(), 1.0);
    KRATOS_EXPECT_VECTOR_NEAR(GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(
                                  unused_strain_vectors, properties),
                              expected_factors, tolerance)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PermeabilityUpdateFactorEqualsOneWhenChangeInverseFactorIsNonPositive)
{
    const auto unused_strain_vectors               = std::vector<Vector>(3, Vector{});
    auto       properties                          = Properties{};
    properties[PERMEABILITY_CHANGE_INVERSE_FACTOR] = -1.0;

    const auto expected_factors = std::vector<double>(unused_strain_vectors.size(), 1.0);
    KRATOS_EXPECT_VECTOR_NEAR(GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(
                                  unused_strain_vectors, properties),
                              expected_factors, tolerance)

    properties[PERMEABILITY_CHANGE_INVERSE_FACTOR] = 0.0;

    KRATOS_EXPECT_VECTOR_NEAR(GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(
                                  unused_strain_vectors, properties),
                              expected_factors, tolerance)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       PermeabilityUpdateFactorIsComputedFromStrainsAndPropertiesWhenChangeInverseFactorIsPositive)
{
    auto test_strains = Vector{3};
    test_strains <<= 0.001, 0.002, 0.0;
    auto strain_vectors = std::vector<Vector>{test_strains, 2.0 * test_strains, 4.0 * test_strains};

    auto properties                                = Properties{};
    properties[PERMEABILITY_CHANGE_INVERSE_FACTOR] = 0.5;
    properties[POROSITY]                           = 0.2;

    const auto expected_factors = std::vector<double>{1.00433, 1.0087, 1.01753};
    KRATOS_EXPECT_VECTOR_NEAR(
        GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(strain_vectors, properties),
        expected_factors, 1e-5)
}

} // namespace Kratos::Testing
