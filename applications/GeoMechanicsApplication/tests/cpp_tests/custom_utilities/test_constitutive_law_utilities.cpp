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
#include "geo_mechanics_application_variables.h"
#include "includes/checks.h"
#include "tests/cpp_tests/custom_constitutive/mock_constitutive_law.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <numbers>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SetSixConstitutiveParametersCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    ConstitutiveLaw::Parameters ConstitutiveParameters;

    Vector strain_vector(3);
    strain_vector <<= 1.0, 2.0, 3.0;
    Matrix constitutive_matrix = IdentityMatrix(5, 5);
    Vector N(3);
    N <<= 0.1, 0.2, 0.5;
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
    auto properties      = Properties{};
    auto umat_parameters = Vector{2};
    umat_parameters <<= 2.0, 30.0;
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
    auto properties      = Properties{};
    auto umat_parameters = Vector{2};
    umat_parameters <<= 2.0, 30.0;
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

    auto properties      = Properties{};
    auto umat_parameters = Vector{2};
    umat_parameters <<= 2.0, 30.0;
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
    using namespace std::numbers;

    EXPECT_NEAR(ConstitutiveLawUtilities::CalculateK0NCFromFrictionAngleInDegrees(30.0), 0.5,
                Defaults::absolute_tolerance);
    EXPECT_NEAR(ConstitutiveLawUtilities::CalculateK0NCFromFrictionAngleInDegrees(60.0),
                1.0 - 0.5 * sqrt3, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
