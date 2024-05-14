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

#include "includes/checks.h"
#include "includes/constitutive_law.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SetSixConstitutiveParametersCorrectResults, KratosGeoMechanicsFastSuite)
{
    ConstitutiveLaw::Parameters ConstitutiveParameters;

    Vector strain_vector(3);
    strain_vector(0)           = 1.0;
    strain_vector(1)           = 2.0;
    strain_vector(2)           = 3.0;
    Matrix constitutive_matrix = IdentityMatrix(5, 5);
    Vector N(3);
    N(0)                                         = 0.1;
    N(1)                                         = 0.2;
    N(2)                                         = 0.5;
    const Matrix     shape_functions_derivatives = ScalarMatrix(3, 3, 5.0);
    const Matrix     deformation_gradient_F      = ScalarMatrix(3, 3, 10.0);
    constexpr double determinant_of_F            = 10.0;
    ConstitutiveLaw::Parameters::SetSixConstitutiveParameters(
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

} // namespace Kratos::Testing