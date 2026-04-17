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

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilities_ReplaceIgnoreUndrainedByDrainageTypeEndsWithoutIgnoreUndrained,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties* pProperties = nullptr;
    EXPECT_NO_THROW(ConstitutiveLawUtilities::ReplaceIgnoreUndrainedByDrainageType(pProperties));

    Model               my_model;
    ModelPart&          my_model_part = my_model.CreateModelPart("Main");
    Properties::Pointer p_properties  = my_model_part.CreateNewProperties(0);
    p_properties->SetValue(IGNORE_UNDRAINED, true);
    ConstitutiveLawUtilities::ReplaceIgnoreUndrainedByDrainageType(p_properties.get());
    EXPECT_FALSE(p_properties->Has(IGNORE_UNDRAINED));
    EXPECT_TRUE(p_properties->Has(GEO_DRAINAGE_TYPE));

    // would be nice to test the warning and the value of GEO_DRAINAGE_TYPE
}

} // namespace Kratos::Testing
