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
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

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
        "Material 0 does not have UMAT_PARAMETERS and/or INDEX_OF_UMAT_C_PARAMETER");

    properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 1);
    properties.Erase(UMAT_PARAMETERS);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ConstitutiveLawUtilities::GetCohesion(properties),
        "Material 0 does not have UMAT_PARAMETERS and/or INDEX_OF_UMAT_C_PARAMETER");
}

KRATOS_TEST_CASE_IN_SUITE(FrictionAngleCanBeFetchedFromGeoFrictionAngleProperty, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);

    KRATOS_EXPECT_DOUBLE_EQ(ConstitutiveLawUtilities::GetFrictionAngle(properties), 30.0);
}

KRATOS_TEST_CASE_IN_SUITE(FrictionAngleCanBeFetchedFromUMatParameters, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties      = Properties{};
    auto umat_parameters = Vector{2};
    umat_parameters <<= 2.0, 30.0;
    properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 2);

    KRATOS_EXPECT_DOUBLE_EQ(ConstitutiveLawUtilities::GetFrictionAngle(properties), 30.0);

    properties.Erase(INDEX_OF_UMAT_PHI_PARAMETER);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ConstitutiveLawUtilities::GetFrictionAngle(properties),
        "Material 0 does not have UMAT_PARAMETERS and/or INDEX_OF_UMAT_PHI_PARAMETER");

    properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 2);
    properties.Erase(UMAT_PARAMETERS);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ConstitutiveLawUtilities::GetFrictionAngle(properties),
        "Material 0 does not have UMAT_PARAMETERS and/or INDEX_OF_UMAT_PHI_PARAMETER");
}

// The following test only raises errors when using debug builds
#ifdef KRATOS_DEBUG

KRATOS_TEST_CASE_IN_SUITE(ThrowExceptionWhenIndexInUMatParametersIsOutOfBounds, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties      = Properties{};
    auto umat_parameters = Vector{2};
    umat_parameters <<= 2.0, 30.0;
    properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 0); // 1-based index

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ConstitutiveLawUtilities::GetCohesion(properties),
        "Got out-of-bounds INDEX_OF_UMAT_C_PARAMETER (material ID: 0): 0 is not in range [1, 2]");

    properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 3); // 1-based index
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ConstitutiveLawUtilities::GetCohesion(properties),
        "Got out-of-bounds INDEX_OF_UMAT_C_PARAMETER (material ID: 0): 3 is not in range [1, 2]");
}

#endif

} // namespace Kratos::Testing
