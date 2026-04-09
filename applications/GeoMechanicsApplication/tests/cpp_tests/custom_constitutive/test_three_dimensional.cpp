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
#include "custom_constitutive/three_dimensional.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensional_GetStrainSizeReturnsSix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_three_dimensional = std::make_unique<ThreeDimensional>();

    // Act & Assert
    KRATOS_EXPECT_EQ(p_three_dimensional->GetStrainSize(), VOIGT_SIZE_3D);
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensional_GetDimensionReturnsThree, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_three_dimensional = std::make_unique<ThreeDimensional>();

    // Act & Assert
    KRATOS_EXPECT_EQ(p_three_dimensional->GetDimension(), N_DIM_3D);
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensional_GetNumberOfNormalComponentsReturnsThree, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_three_dimensional = std::make_unique<ThreeDimensional>();

    // Act & Assert
    KRATOS_EXPECT_EQ(p_three_dimensional->GetNumberOfNormalComponents(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensional_GetSpatialTypeReturnsThreeDimensionalLaw, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_three_dimensional = std::make_unique<ThreeDimensional>();

    // Act & Assert
    KRATOS_EXPECT_EQ(p_three_dimensional->GetSpatialType(), ConstitutiveLaw::THREE_DIMENSIONAL_LAW);
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensional_CalculateElasticMatrixReturnsConstitutiveTensor,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_three_dimensional = std::make_unique<ThreeDimensional>();
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0);
    properties.SetValue(POISSON_RATIO, 0.0);

    // Act
    auto constitutive_tensor = p_three_dimensional->CalculateElasticConstitutiveTensor(properties);

    // Assert
    auto expected_tensor = UblasUtilities::CreateMatrix({{1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                                         {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                                                         {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                                                         {0.0, 0.0, 0.0, 0.5, 0.0, 0.0},
                                                         {0.0, 0.0, 0.0, 0.0, 0.5, 0.0},
                                                         {0.0, 0.0, 0.0, 0.0, 0.0, 0.5}});
    KRATOS_EXPECT_MATRIX_EQ(constitutive_tensor, expected_tensor);

    // Arrange
    properties.SetValue(GEO_DRAINAGE_TYPE, static_cast<int>(DrainageType::UNDRAINED));
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.5);

    // Act
    constitutive_tensor = p_three_dimensional->CalculateElasticConstitutiveTensor(properties);

    // Assert
    expected_tensor = UblasUtilities::CreateMatrix({{4.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0},
                                                    {1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0},
                                                    {1.0 / 3.0, 1.0 / 3.0, 4.0 / 3.0, 0.0, 0.0, 0.0},
                                                    {0.0, 0.0, 0.0, 0.5, 0.0, 0.0},
                                                    {0.0, 0.0, 0.0, 0.0, 0.5, 0.0},
                                                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.5}});
    KRATOS_EXPECT_MATRIX_EQ(constitutive_tensor, expected_tensor);

    // Arrange
    properties.Erase(BULK_MODULUS_FLUID);
    properties.Erase(BULK_MODULUS_SOLID);
    properties.Erase(POROSITY);
    properties.SetValue(SKEMPTON_B, 0.5);

    // Act
    constitutive_tensor = p_three_dimensional->CalculateElasticConstitutiveTensor(properties);

    // Assert
    KRATOS_EXPECT_MATRIX_EQ(constitutive_tensor, expected_tensor);

    // Arrange
    properties.Erase(SKEMPTON_B);
    properties.Erase(BIOT_COEFFICIENT);
    properties.SetValue(POISSON_UNDRAINED, 0.2);

    // Act
    constitutive_tensor = p_three_dimensional->CalculateElasticConstitutiveTensor(properties);

    // Assert
    KRATOS_EXPECT_MATRIX_EQ(constitutive_tensor, expected_tensor);
}

} // namespace Kratos::Testing