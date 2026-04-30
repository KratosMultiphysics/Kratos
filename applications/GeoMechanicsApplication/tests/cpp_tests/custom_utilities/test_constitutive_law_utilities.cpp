// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/checks.h"
#include "tests/cpp_tests/custom_constitutive/mock_constitutive_law.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <numbers>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SetSixConstitutiveParametersCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    ConstitutiveLaw::Parameters ConstitutiveParameters;

    auto             strain_vector               = UblasUtilities::CreateVector({1.0, 2.0, 3.0});
    Matrix           constitutive_matrix         = IdentityMatrix(5, 5);
    const auto       N                           = UblasUtilities::CreateVector({0.1, 0.2, 0.5});
    const Matrix     shape_functions_derivatives = ScalarMatrix(3, 3, 5.0);
    const Matrix     deformation_gradient_F      = ScalarMatrix(3, 3, 10.0);
    constexpr double determinant_of_F            = 10.0;
    ConstitutiveLawUtilities::SetConstitutiveParameters(
        ConstitutiveParameters, strain_vector, constitutive_matrix, N, shape_functions_derivatives,
        deformation_gradient_F, determinant_of_F);

    KRATOS_CHECK_VECTOR_NEAR(ConstitutiveParameters.GetStrainVector(), strain_vector, 1e-12)

    KRATOS_CHECK_MATRIX_NEAR(ConstitutiveParameters.GetConstitutiveMatrix(), constitutive_matrix, 1e-12)

    KRATOS_CHECK_VECTOR_NEAR(ConstitutiveParameters.GetShapeFunctionsValues(), N, 1e-12)

    KRATOS_CHECK_MATRIX_NEAR(ConstitutiveParameters.GetShapeFunctionsDerivatives(),
                             shape_functions_derivatives, 1e-12)

    KRATOS_CHECK_MATRIX_NEAR(ConstitutiveParameters.GetDeformationGradientF(), deformation_gradient_F, 1e-12)

    KRATOS_CHECK_NEAR(ConstitutiveParameters.GetDeterminantF(), determinant_of_F, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(CohesionCanBeFetchedFromGeoCohesionProperty, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    properties.SetValue(GEO_COHESION, 2.0);

    KRATOS_EXPECT_DOUBLE_EQ(ConstitutiveLawUtilities::GetCohesion(properties), 2.0);
}

KRATOS_TEST_CASE_IN_SUITE(CohesionCanBeFetchedFromUMatParameters, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       properties      = Properties{};
    const auto umat_parameters = UblasUtilities::CreateVector({2.0, 30.0});
    properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 1);

    KRATOS_EXPECT_DOUBLE_EQ(ConstitutiveLawUtilities::GetCohesion(properties), 2.0);

    properties.Erase(INDEX_OF_UMAT_C_PARAMETER);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ConstitutiveLawUtilities::GetCohesion(properties),
        "Error: ConstitutiveLawUtilities::GetCohesion failed. There is no GEO_COHESION available "
        "and attempting to get the cohesion from UMAT parameters resulted in the following Error: "
        "There is no INDEX_OF_UMAT_C_PARAMETER for material 0.");

    properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 1);
    properties.Erase(UMAT_PARAMETERS);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ConstitutiveLawUtilities::GetCohesion(properties),
        "ConstitutiveLawUtilities::GetCohesion failed. There is no GEO_COHESION available and "
        "attempting to get the cohesion from UMAT parameters resulted in the following Error: "
        "There is no UMAT_PARAMETERS for material 0.");
}

KRATOS_TEST_CASE_IN_SUITE(FrictionAngleCanBeFetchedFromGeoFrictionAngleProperty, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);

    KRATOS_EXPECT_DOUBLE_EQ(ConstitutiveLawUtilities::GetFrictionAngleInDegrees(properties), 30.0);
}

KRATOS_TEST_CASE_IN_SUITE(FrictionAngleCanBeFetchedFromUMatParameters, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       properties      = Properties{};
    const auto umat_parameters = UblasUtilities::CreateVector({2.0, 30.0});
    properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 2);

    KRATOS_EXPECT_DOUBLE_EQ(ConstitutiveLawUtilities::GetFrictionAngleInDegrees(properties), 30.0);

    properties.Erase(INDEX_OF_UMAT_PHI_PARAMETER);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ConstitutiveLawUtilities::GetFrictionAngleInDegrees(properties), "ConstitutiveLawUtilities::GetFrictionAngleInDegrees failed. There is no GEO_FRICTION_ANGLE available and attempting to get the friction angle from UMAT parameters resulted in the following Error: There is no INDEX_OF_UMAT_PHI_PARAMETER for material 0.");

    properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 2);
    properties.Erase(UMAT_PARAMETERS);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ConstitutiveLawUtilities::GetFrictionAngleInDegrees(properties), "ConstitutiveLawUtilities::GetFrictionAngleInDegrees failed. There is no GEO_FRICTION_ANGLE available and attempting to get the friction angle from UMAT parameters resulted in the following Error: There is no UMAT_PARAMETERS for material 0.");
}

KRATOS_TEST_CASE_IN_SUITE(RaiseADebugErrorWhenIndexInUMatParametersIsOutOfBounds, KratosGeoMechanicsFastSuiteWithoutKernel)
{
#ifndef KRATOS_DEBUG
    GTEST_SKIP() << "This test requires a debug build";
#endif

    auto       properties      = Properties{};
    const auto umat_parameters = UblasUtilities::CreateVector({2.0, 30.0});
    properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 0); // 1-based index

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ConstitutiveLawUtilities::GetCohesion(properties),
        "ConstitutiveLawUtilities::GetCohesion failed. There is no GEO_COHESION available and "
        "attempting to get the cohesion from UMAT parameters resulted in the following Error: Got "
        "out-of-bounds INDEX_OF_UMAT_C_PARAMETER (material ID: 0): 0 is not in range [1, 2].");

    properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 3); // 1-based index
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ConstitutiveLawUtilities::GetCohesion(properties),
        "ConstitutiveLawUtilities::GetCohesion failed. There is no GEO_COHESION available and "
        "attempting to get the cohesion from UMAT parameters resulted in the following Error: Got "
        "out-of-bounds INDEX_OF_UMAT_C_PARAMETER (material ID: 0): 3 is not in range [1, 2].");
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_GetStateVariableIndex, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Act and Assert
    constexpr auto expected_index = 1;
    KRATOS_EXPECT_EQ(ConstitutiveLawUtilities::GetStateVariableIndex(STATE_VARIABLE_2), expected_index);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_CheckStrainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties                       = Properties{};
    properties.GetValue(CONSTITUTIVE_LAW) = Kratos::make_shared<MockConstitutiveLaw>();

    // Act and Assert
    std::size_t           expected_size = 2;
    constexpr std::size_t element_id    = 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ConstitutiveLawUtilities::CheckStrainSize(properties, expected_size, element_id), "Wrong constitutive law is used: strain size is 4 when it is expected to be 2 at element Id = 1.");

    expected_size = properties[CONSTITUTIVE_LAW]->GetStrainSize();
    EXPECT_NO_THROW(ConstitutiveLawUtilities::CheckStrainSize(properties, expected_size, element_id));
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_CheckStrainMeasures, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties                       = Properties{};
    auto constitutive_law                 = Kratos::make_shared<MockConstitutiveLaw>();
    properties.GetValue(CONSTITUTIVE_LAW) = constitutive_law;

    // Act and Assert
    constexpr std::size_t element_id = 1;
    constitutive_law->AddStrainMeasure_Infinitesimal(false);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ConstitutiveLawUtilities::CheckHasStrainMeasure_Infinitesimal(properties, element_id), " Constitutive law is not compatible with the strain type StrainMeasure_Infinitesimal at element 1.");

    constitutive_law->AddStrainMeasure_Infinitesimal(true);
    EXPECT_NO_THROW(ConstitutiveLawUtilities::CheckHasStrainMeasure_Infinitesimal(properties, element_id));
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_CalculateK0NCFromFrictionAngleInDegrees,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    EXPECT_NEAR(ConstitutiveLawUtilities::CalculateK0NCFromFrictionAngleInRadians(MathUtils<>::DegreesToRadians(30.0)),
                0.5, Defaults::absolute_tolerance);
    EXPECT_NEAR(ConstitutiveLawUtilities::CalculateK0NCFromFrictionAngleInRadians(MathUtils<>::DegreesToRadians(60.0)),
                1.0 - 0.5 * std::numbers::sqrt3, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_HasFrictionAngle, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // 1) GEO_FRICTION_ANGLE provided
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);
    KRATOS_EXPECT_EQ(ConstitutiveLawUtilities::HasFrictionAngle(properties), true);

    // 2) UMAT_PARAMETERS + INDEX_OF_UMAT_PHI_PARAMETER provided
    properties           = Properties{};
    auto umat_parameters = UblasUtilities::CreateVector({2.0, 30.0});
    properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 2);
    KRATOS_EXPECT_EQ(ConstitutiveLawUtilities::HasFrictionAngle(properties), true);

    // 3) only INDEX_OF_UMAT_PHI_PARAMETER provided -> should be false
    properties = Properties{};
    properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 1);
    KRATOS_EXPECT_EQ(ConstitutiveLawUtilities::HasFrictionAngle(properties), false);

    // 4) only UMAT_PARAMETERS provided -> should be false
    properties      = Properties{};
    umat_parameters = UblasUtilities::CreateVector({2.0, 30.0});
    properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    KRATOS_EXPECT_EQ(ConstitutiveLawUtilities::HasFrictionAngle(properties), false);

    // 5) neither provided -> false
    properties = Properties{};
    KRATOS_EXPECT_EQ(ConstitutiveLawUtilities::HasFrictionAngle(properties), false);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_ValidateFrictionAngle, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr std::size_t element_id = 1;

    // Valid: GEO_FRICTION_ANGLE provided
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);
    EXPECT_NO_THROW(ConstitutiveLawUtilities::ValidateFrictionAngle(properties, element_id));

    // Valid: UMAT_PARAMETERS + INDEX_OF_UMAT_PHI_PARAMETER provided
    properties           = Properties{};
    auto umat_parameters = UblasUtilities::CreateVector({2.0, 30.0});
    properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 2);
    EXPECT_NO_THROW(ConstitutiveLawUtilities::ValidateFrictionAngle(properties, element_id));

    // Missing both -> error
    properties = Properties{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ConstitutiveLawUtilities::ValidateFrictionAngle(properties, element_id), "Properties ( 0) of element ( 1) does not have GEO_FRICTION_ANGLE nor INDEX_OF_UMAT_PHI_PARAMETER.");

    // GEO_FRICTION_ANGLE out of range -> error
    properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ConstitutiveLawUtilities::ValidateFrictionAngle(properties, element_id), " Properties ( 0) of element ( 1): GEO_FRICTION_ANGLE (-1 degrees) should be between 0 and 90 degrees.");

    // UMAT phi value out of range -> error
    properties      = Properties{};
    umat_parameters = UblasUtilities::CreateVector({2.0, -5.0});
    properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ConstitutiveLawUtilities::ValidateFrictionAngle(properties, element_id),
        " Properties ( 0) of element ( 1): Phi (-5 degrees) should be between 0 and 90 degrees.");

    // INDEX_OF_UMAT_PHI_PARAMETER out of bounds -> error
    properties      = Properties{};
    umat_parameters = UblasUtilities::CreateVector({2.0, 30.0});
    properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ConstitutiveLawUtilities::ValidateFrictionAngle(properties, element_id), "Properties ( 0) of element ( 1): INDEX_OF_UMAT_PHI_PARAMETER (3) is not in range [1, size of UMAT_PARAMETERS].");
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_CalculateUndrainedYoungsModulus_ReturnsExpectedValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0);
    properties.SetValue(POISSON_RATIO, 0.2);

    // Act
    const auto undrained_youngs_modulus =
        ConstitutiveLawUtilities::CalculateUndrainedYoungsModulus(properties, 0.4);

    // Assert
    KRATOS_EXPECT_EQ(undrained_youngs_modulus, 7.0 / 6.0);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_CalculateElasticProperties, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange: plain properties
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 2.5);
    properties.SetValue(POISSON_RATIO, 0.25);

    // Act
    const auto [E_plain, nu_plain] = ConstitutiveLawUtilities::GetOrCalculateElasticProperties(properties);

    // Assert
    KRATOS_EXPECT_EQ(E_plain, 2.5);
    KRATOS_EXPECT_EQ(nu_plain, 0.25);

    // Arrange: properties for undrained computation
    properties = Properties{};
    properties.SetValue(YOUNG_MODULUS, 1.0);
    properties.SetValue(POISSON_RATIO, 0.2);
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.5);

    // Act: undrained
    const auto [E_undrained, nu_undrained] =
        ConstitutiveLawUtilities::GetOrCalculateElasticProperties(properties, true);

    // Assert
    const auto expected_E  = 10.0 / 9.0;
    const auto expected_nu = 1.0 / 3.0;
    KRATOS_EXPECT_NEAR(E_undrained, expected_E, Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(nu_undrained, expected_nu, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_GetOrCalculatedSkemptonB_ReturnsExpectedValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0);
    properties.SetValue(POISSON_RATIO, 0.2);
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.5);

    // Act
    auto skempton_b = ConstitutiveLawUtilities::GetOrCalculateSkemptonB(properties);

    // Assert
    KRATOS_EXPECT_EQ(skempton_b, 0.5);

    // Arrange
    properties.SetValue(BIOT_COEFFICIENT, 1.0e-8);
    properties.SetValue(POROSITY, 1.0e-8);

    // Act
    skempton_b = ConstitutiveLawUtilities::GetOrCalculateSkemptonB(properties);

    // Assert
    constexpr auto tolerance = 1.0e6 * Defaults::absolute_tolerance;
    KRATOS_EXPECT_NEAR(skempton_b, 0.5, tolerance);

    // Arrange
    properties.SetValue(BIOT_COEFFICIENT, 0.0);

    // Act
    skempton_b = ConstitutiveLawUtilities::GetOrCalculateSkemptonB(properties);

    // Assert
    KRATOS_EXPECT_EQ(skempton_b, 0.0);

    // Arrange
    properties.SetValue(GEO_SKEMPTON_B, 0.8);
    // Act
    skempton_b = ConstitutiveLawUtilities::GetOrCalculateSkemptonB(properties);

    // Assert
    KRATOS_EXPECT_EQ(skempton_b, 0.8);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_GetOrCalculatedSkemptonB_ThrowsError,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(BIOT_COEFFICIENT, 0.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.0);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((void)ConstitutiveLawUtilities::GetOrCalculateSkemptonB(properties),
                                      "Non-physical values: denominator < epsilon.");

    // Arrange
    properties.SetValue(BIOT_COEFFICIENT, 0.5);
    properties.SetValue(BULK_MODULUS_SOLID, 1.E2);
    properties.SetValue(POROSITY, 0.1);
    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((void)ConstitutiveLawUtilities::GetOrCalculateSkemptonB(properties),
                                      "Calculated Skempton B (1.08696) is out of range [0,1].");

    // Arrange
    properties.SetValue(BIOT_COEFFICIENT, -0.0001);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.1);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((void)ConstitutiveLawUtilities::GetOrCalculateSkemptonB(properties),
                                      "Calculated Skempton B (-0.0010011) is out of range [0,1].");
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_CalculateUndrainedPoissonsRatio_ReturnsExpectedValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(POISSON_RATIO, 0.2);
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.5);

    // Act
    auto undrained_poissons_ratio = ConstitutiveLawUtilities::CalculateUndrainedPoissonsRatio(properties);

    // Assert
    KRATOS_EXPECT_NEAR(undrained_poissons_ratio, 1.0 / 3.0, Defaults::absolute_tolerance);

    // Arrange
    properties.SetValue(GEO_SKEMPTON_B, 0.2);

    // Act
    undrained_poissons_ratio = ConstitutiveLawUtilities::CalculateUndrainedPoissonsRatio(properties);

    // Assert
    KRATOS_EXPECT_NEAR(undrained_poissons_ratio, 0.25, Defaults::absolute_tolerance);

    // Arrange
    properties.SetValue(POISSON_RATIO, 0.3);
    properties.SetValue(GEO_SKEMPTON_B, 1.0);

    // Act
    undrained_poissons_ratio = ConstitutiveLawUtilities::CalculateUndrainedPoissonsRatio(properties);

    // Assert
    KRATOS_EXPECT_NEAR(undrained_poissons_ratio, 0.5, Defaults::absolute_tolerance);

    // Arrange
    properties.SetValue(BIOT_COEFFICIENT, 0.0);

    // Act
    undrained_poissons_ratio = ConstitutiveLawUtilities::CalculateUndrainedPoissonsRatio(properties);

    // Assert
    KRATOS_EXPECT_NEAR(undrained_poissons_ratio, 0.3, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_GetOrCalculateUndrainedPoissonsRatio_ReturnsExpectedValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(POISSON_RATIO, 0.2);
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.5);

    // Act
    auto undrained_poissons_ratio = ConstitutiveLawUtilities::GetOrCalculateUndrainedPoissonsRatio(properties);

    // Assert
    const auto expected_undrained_poissons_ratio = 1.0 / 3.0;
    KRATOS_EXPECT_NEAR(undrained_poissons_ratio, expected_undrained_poissons_ratio, Defaults::absolute_tolerance);

    // When explicit GEO_POISSON_UNDRAINED is provided it should be returned (and clamped if needed)
    properties.SetValue(GEO_POISSON_UNDRAINED, 0.4);
    undrained_poissons_ratio = ConstitutiveLawUtilities::GetOrCalculateUndrainedPoissonsRatio(properties);
    KRATOS_EXPECT_NEAR(undrained_poissons_ratio, 0.4, Defaults::absolute_tolerance);

    properties.SetValue(GEO_POISSON_UNDRAINED, 0.6);
    undrained_poissons_ratio = ConstitutiveLawUtilities::GetOrCalculateUndrainedPoissonsRatio(properties);
    KRATOS_EXPECT_NEAR(undrained_poissons_ratio, 0.495, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_CalculateUndrainedPoissonsRatio_ThrowsError,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(POISSON_RATIO, -0.99999999999999995);
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(GEO_SKEMPTON_B, 1.0);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((void)ConstitutiveLawUtilities::CalculateUndrainedPoissonsRatio(properties),
                                      "Non-physical values: denominator < epsilon.");
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitilties_MakeContinuumConstitutiveTensorReturnsConstitutiveTensor,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange & Act
    constexpr auto youngs_modulus              = 1.0;
    constexpr auto poissons_ratio              = 0.25;
    constexpr auto strain_size                 = 4;
    constexpr auto number_of_normal_components = 2;
    const auto constitutive_tensor = ConstitutiveLawUtilities::MakeContinuumElasticConstitutiveTensor(
        youngs_modulus, poissons_ratio, strain_size, number_of_normal_components);

    // Assert
    const auto expected_tensor = UblasUtilities::CreateMatrix(
        {{1.2, 0.4, 0.0, 0.0}, {0.4, 1.2, 0.0, 0.0}, {0.0, 0.0, 0.4, 0.0}, {0.0, 0.0, 0.0, 0.4}});
    KRATOS_EXPECT_MATRIX_EQ(constitutive_tensor, expected_tensor);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_MakeInterfaceElasticConstitutiveTensorReturnsConstitutiveTensor,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange & Act
    constexpr auto youngs_modulus              = 2.0;
    constexpr auto poissons_ratio              = 1.0;
    constexpr auto strain_size                 = 4;
    constexpr auto number_of_normal_components = 2;
    const auto constitutive_tensor = ConstitutiveLawUtilities::MakeInterfaceElasticConstitutiveTensor(
        youngs_modulus, poissons_ratio, strain_size, number_of_normal_components);

    // Assert
    const auto expected_tensor = UblasUtilities::CreateMatrix(
        {{2.0, 0.0, 0.0, 0.0}, {0.0, 2.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}});
    KRATOS_EXPECT_MATRIX_EQ(constitutive_tensor, expected_tensor);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_CalculateK0NCFromFrictionAngleGivesK0NC,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange & Act
    const auto calculated_k0_nc =
        ConstitutiveLawUtilities::CalculateK0NCFromFrictionAngleInRadians(30.0 * Globals::Pi / 180.0);

    // Assert
    KRATOS_EXPECT_EQ(calculated_k0_nc, 0.5);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_ReplaceIgnoreUndrainedByDrainageTypeEndsWithoutIgnoreUndrained,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model               my_model;
    ModelPart&          my_model_part = my_model.CreateModelPart("Main");
    Properties::Pointer p_properties  = my_model_part.CreateNewProperties(0);
    p_properties->SetValue(IGNORE_UNDRAINED, true);
    // Act
    ConstitutiveLawUtilities::ReplaceIgnoreUndrainedByDrainageType(*p_properties);
    // Assert
    EXPECT_FALSE(p_properties->Has(IGNORE_UNDRAINED));
    EXPECT_TRUE(p_properties->Has(GEO_DRAINAGE_TYPE));
    KRATOS_EXPECT_EQ(ConstitutiveLawUtilities::StringToDrainageType(p_properties->GetValue(GEO_DRAINAGE_TYPE)),
                     DrainageType::CONSTANT_WATER_PRESSURE);

    // Arrange
    p_properties->Erase(GEO_DRAINAGE_TYPE);
    p_properties->SetValue(IGNORE_UNDRAINED, false);
    // Act
    ConstitutiveLawUtilities::ReplaceIgnoreUndrainedByDrainageType(*p_properties);
    // Assert
    EXPECT_FALSE(p_properties->Has(IGNORE_UNDRAINED));
    EXPECT_TRUE(p_properties->Has(GEO_DRAINAGE_TYPE));
    KRATOS_EXPECT_EQ(ConstitutiveLawUtilities::StringToDrainageType(p_properties->GetValue(GEO_DRAINAGE_TYPE)),
                     DrainageType::FULLY_COUPLED);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_CalculateExcessPorePressureIncrementGivesDeltaPw,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.5);
    constexpr auto volumetric_strain_increment = 1.0;

    // Act
    auto excess_pore_pressure_increment = ConstitutiveLawUtilities::CalculateExcessPorePressureIncrement(
        properties, volumetric_strain_increment);

    // Assert
    KRATOS_EXPECT_NEAR(excess_pore_pressure_increment, 4.0e3 / 3.0, Defaults::absolute_tolerance);

    // Arrange
    properties.SetValue(POROSITY, 0.0);

    // Act
    excess_pore_pressure_increment = ConstitutiveLawUtilities::CalculateExcessPorePressureIncrement(
        properties, volumetric_strain_increment);

    // Assert
    KRATOS_EXPECT_NEAR(excess_pore_pressure_increment, properties.GetValue(BULK_MODULUS_SOLID),
                       Defaults::absolute_tolerance);

    // Arrange
    properties.SetValue(BIOT_COEFFICIENT, 0.0);
    properties.SetValue(POROSITY, 0.5);

    // Act
    excess_pore_pressure_increment = ConstitutiveLawUtilities::CalculateExcessPorePressureIncrement(
        properties, volumetric_strain_increment);

    // Assert
    KRATOS_EXPECT_NEAR(excess_pore_pressure_increment, 0.0, Defaults::absolute_tolerance);

    // Arrange
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(POROSITY, 1.0);

    // Act
    excess_pore_pressure_increment = ConstitutiveLawUtilities::CalculateExcessPorePressureIncrement(
        properties, volumetric_strain_increment);

    // Assert
    KRATOS_EXPECT_NEAR(excess_pore_pressure_increment, properties.GetValue(BULK_MODULUS_FLUID),
                       Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_CalculateExcessPorePressureIncrementThrowsErrorWhenDenominatorEqualsZero,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(BIOT_COEFFICIENT, 0.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.0);
    constexpr auto volumetric_strain_increment = 1.0;

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (void)ConstitutiveLawUtilities::CalculateExcessPorePressureIncrement(properties, volumetric_strain_increment),
        "Non-physical values: denominator < epsilon for property Id of 0.");
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_GetNumberOfNormalStrainComponents_ThreeDimensional,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    auto p_law      = Kratos::make_shared<MockConstitutiveLaw>();
    p_law->SetThreeDimensionalLaw(true);
    p_law->SetStrainSize(6);
    properties.SetValue(CONSTITUTIVE_LAW, p_law);

    KRATOS_EXPECT_EQ(ConstitutiveLawUtilities::GetNumberOfNormalStrainComponents(properties), 3);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_GetNumberOfNormalStrainComponents_PlaneStrain,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    auto p_law      = Kratos::make_shared<MockConstitutiveLaw>();
    p_law->SetPlaneStrainLaw(true);
    p_law->SetStrainSize(4);
    properties.SetValue(CONSTITUTIVE_LAW, p_law);

    KRATOS_EXPECT_EQ(ConstitutiveLawUtilities::GetNumberOfNormalStrainComponents(properties), 3);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_CalculateVolumetricStrain_UsesNormalComponentsOnly,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    auto p_law      = Kratos::make_shared<MockConstitutiveLaw>();
    p_law->SetThreeDimensionalLaw(true);
    p_law->SetStrainSize(6);
    properties.SetValue(CONSTITUTIVE_LAW, p_law);

    const auto strain_vector = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 100.0, 200.0, 300.0});
    KRATOS_EXPECT_DOUBLE_EQ(ConstitutiveLawUtilities::CalculateVolumetricStrain(strain_vector, properties), 6.0);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_CalculateVolumetricStrain_ThrowsForInvalidStrainSize,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    auto p_law      = Kratos::make_shared<MockConstitutiveLaw>();
    p_law->SetThreeDimensionalLaw(true);
    p_law->SetStrainSize(6);
    properties.SetValue(CONSTITUTIVE_LAW, p_law);

    const auto too_small_strain_vector = UblasUtilities::CreateVector({1.0, 2.0});
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (void)ConstitutiveLawUtilities::CalculateVolumetricStrain(too_small_strain_vector, properties),
        "NumberOfNormalComponents (3) exceeds strain vector size (2).");
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_CalculateExcessPorePressureForce_ReturnsExpectedValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    auto p_law      = Kratos::make_shared<MockConstitutiveLaw>();
    p_law->SetPlaneStrainLaw(true);
    p_law->SetStrainSize(4);
    properties.SetValue(CONSTITUTIVE_LAW, p_law);

    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.0e3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.0e3);
    properties.SetValue(POROSITY, 0.5);

    const auto strain_vector = UblasUtilities::CreateVector({0.3, 0.1, 0.1, 99.0});
    const auto voigt_vector  = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    const auto B =
        UblasUtilities::CreateMatrix({{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}, {0.0, 0.0, 0.0}});
    const auto     excess_pore_pressure_previous = UblasUtilities::CreateVector({0.05});
    constexpr auto integration_coefficient       = 2.0;

    const auto force = ConstitutiveLawUtilities::CalculateExcessPorePressureForce(
        properties, strain_vector, B, voigt_vector, integration_coefficient, 0, excess_pore_pressure_previous);

    const auto expected_force = UblasUtilities::CreateVector({1200.0, 2400.0, 3600.0});
    KRATOS_EXPECT_VECTOR_NEAR(force, expected_force, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_CalculateExcessPorePressureForce_ThrowsForInvalidIntegrationPoint,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    auto p_law      = Kratos::make_shared<MockConstitutiveLaw>();
    p_law->SetPlaneStrainLaw(true);
    p_law->SetStrainSize(4);
    properties.SetValue(CONSTITUTIVE_LAW, p_law);

    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.0e3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.0e3);
    properties.SetValue(POROSITY, 0.5);

    const auto strain_vector = UblasUtilities::CreateVector({0.3, 0.1, 0.1, 0.0});
    const auto voigt_vector  = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    const auto B =
        UblasUtilities::CreateMatrix({{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}, {0.0, 0.0, 0.0}});
    const auto excess_pore_pressure_previous = UblasUtilities::CreateVector({0.05});

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (void)ConstitutiveLawUtilities::CalculateExcessPorePressureForce(
            properties, strain_vector, B, voigt_vector, 1.0, 1, excess_pore_pressure_previous),
        "Integration point index (1) exceeds cached previous volumetric strain size (1).");
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_AssembleExcessPorePressureForces_WithUndrainedMaterial,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    auto p_law      = Kratos::make_shared<MockConstitutiveLaw>();
    p_law->SetPlaneStrainLaw(true);
    p_law->SetStrainSize(4);
    properties.SetValue(CONSTITUTIVE_LAW, p_law);
    properties.SetValue(GEO_DRAINAGE_TYPE, "Undrained");

    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.0e3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.0e3);
    properties.SetValue(POROSITY, 0.5);

    // Two integration points with different strains
    auto strain_vectors = std::vector<Vector>{UblasUtilities::CreateVector({0.3, 0.1, 0.1, 99.0}),
                                              UblasUtilities::CreateVector({0.2, 0.15, 0.05, 99.0})};

    auto B_matrices = std::vector<Matrix>{
        UblasUtilities::CreateMatrix({{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}, {0.0, 0.0, 0.0}}),
        UblasUtilities::CreateMatrix({{0.5, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.5}, {0.0, 0.0, 0.0}})};

    const auto voigt_vector                  = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    const auto integration_coefficients      = std::vector<double>{2.0, 3.0};
    const auto excess_pore_pressure_previous = UblasUtilities::CreateVector({0.05, 0.03});

    auto result_vector = Vector(6, 0.0);

    ConstitutiveLawUtilities::AssembleExcessPorePressureForces(
        result_vector, properties, strain_vectors, B_matrices, voigt_vector,
        integration_coefficients, excess_pore_pressure_previous);

    // GP 0: volumetric_strain = 0.5, increment = 0.45, delta_p = (1.0 * 0.45 / (0.5/1.0e3 + (1.0-0.5)/2.0e3)) = 600.0
    //       force = [1.0*600.0, 2.0*600.0, 3.0*600.0] * 2.0 = [1200, 2400, 3600]
    // GP 1: volumetric_strain = 0.4, increment = 0.37, delta_p = (1.0 * 0.37 / (0.5/1.0e3 + 0.5/2.0e3)) = 493.33...
    //       force = [0.5*493.33, 1.0*493.33, 1.5*493.33] * 3.0 = [740, 1480, 2220]
    const auto expected =
        UblasUtilities::CreateVector({1200.0 + 740.0, 2400.0 + 1480.0, 3600.0 + 2220.0, 0.0, 0.0, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(result_vector, expected, 5.0); // 5 Pa tolerance for rounding
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_AssembleExcessPorePressureForces_SkipsWhenNotUndrained,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    auto p_law      = Kratos::make_shared<MockConstitutiveLaw>();
    p_law->SetPlaneStrainLaw(true);
    p_law->SetStrainSize(4);
    properties.SetValue(CONSTITUTIVE_LAW, p_law);
    properties.SetValue(GEO_DRAINAGE_TYPE, "Drained"); // Not undrained

    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.0e3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.0e3);
    properties.SetValue(POROSITY, 0.5);

    auto strain_vectors = std::vector<Vector>{UblasUtilities::CreateVector({0.3, 0.1, 0.1, 99.0})};
    auto B_matrices     = std::vector<Matrix>{UblasUtilities::CreateMatrix(
        {{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}, {0.0, 0.0, 0.0}})};

    const auto voigt_vector                  = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    const auto integration_coefficients      = std::vector<double>{2.0};
    const auto excess_pore_pressure_previous = UblasUtilities::CreateVector({0.05});

    auto result_vector = UblasUtilities::CreateVector({0.0, 0.0, 0.0});

    // Should return early without assembling anything
    ConstitutiveLawUtilities::AssembleExcessPorePressureForces(
        result_vector, properties, strain_vectors, B_matrices, voigt_vector,
        integration_coefficients, excess_pore_pressure_previous);

    KRATOS_EXPECT_VECTOR_NEAR(result_vector, UblasUtilities::CreateVector({0.0, 0.0, 0.0}), 1e-12);
}

} // namespace Kratos::Testing
