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
//                   Gennady Markelov
//

#include "custom_utilities/math_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"
#include "utilities/math_utils.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateTrace, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const Vector eye = MathUtils<double>::StressTensorToVector(IdentityMatrix(3));
    KRATOS_EXPECT_DOUBLE_EQ(3.0, StressStrainUtilities::CalculateTrace(eye));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMeanStress, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const Vector eye = MathUtils<double>::StressTensorToVector(IdentityMatrix(3));
    KRATOS_EXPECT_DOUBLE_EQ(1.0, StressStrainUtilities::CalculateMeanStress(eye));
}

KRATOS_TEST_CASE_IN_SUITE(CheckGreenLagrangeStrainTensor, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const Matrix eye =
        StressStrainUtilities::CalculateGreenLagrangeStrainTensor(std::sqrt(3.) * IdentityMatrix(3));
    KRATOS_EXPECT_MATRIX_NEAR(eye, IdentityMatrix(3), Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateVonMisesStressHydrostatic, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector stress_vector(4);
    stress_vector <<= -2.0, -2.0, -2.0, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateVonMisesStress(stress_vector));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateVonMisesStressPureShear, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector stress_vector(4);
    stress_vector <<= 0.0, 0.0, 0.0, 2.0;
    KRATOS_EXPECT_DOUBLE_EQ(2. * std::sqrt(3.), StressStrainUtilities::CalculateVonMisesStress(stress_vector));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateLodeAngle, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector stress_vector(4);
    stress_vector <<= 0.5, -1.0, 0.5, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(30.),
                            StressStrainUtilities::CalculateLodeAngle(stress_vector));
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(-30.),
                            StressStrainUtilities::CalculateLodeAngle(-1. * stress_vector));
    Vector stress_vector2(4);
    stress_vector2 <<= -1.0, 0.0, 1.0, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(0.),
                            StressStrainUtilities::CalculateLodeAngle(stress_vector2));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMohrCoulombShearCapacityZeroStress, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const Vector stress_vector  = ZeroVector(4);
    const double cohesion       = 2.;
    const double friction_angle = 0.;
    KRATOS_EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateMohrCoulombShearCapacity(
                                     stress_vector, cohesion, friction_angle));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMohrCoulombShearCapacityHydrostatic, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector stress_vector(4);
    stress_vector <<= -2.0, -2.0, -2.0, 0.0;
    const double cohesion       = 0.;
    const double friction_angle = MathUtils<>::DegreesToRadians(90.);
    KRATOS_EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateMohrCoulombShearCapacity(
                                     stress_vector, cohesion, friction_angle));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMohrCoulombShearCapacityShearOnly, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector stress_vector(4);
    stress_vector <<= -2.0, 0.0, 2.0, 0.0;
    const double cohesion       = 2.;
    const double friction_angle = 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(1.0, StressStrainUtilities::CalculateMohrCoulombShearCapacity(
                                     stress_vector, cohesion, friction_angle));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMohrCoulombPressureCapacityZeroStress, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector       stress_vector  = ZeroVector(4);
    const double cohesion       = 2.;
    const double friction_angle = 0.;
    KRATOS_EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateMohrCoulombPressureCapacity(
                                     stress_vector, cohesion, friction_angle));
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateMohrCoulombPressureCapacityShearOnly, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector stress_vector(4);
    stress_vector <<= -2.0, 0.0, 2.0, 0.0;
    const double cohesion       = std::sqrt(2. * 2. * 3.);
    const double friction_angle = MathUtils<>::DegreesToRadians(30.0);
    KRATOS_EXPECT_NEAR(3. * std::sin(friction_angle),
                       StressStrainUtilities::CalculateMohrCoulombPressureCapacity(
                           stress_vector, cohesion, friction_angle),
                       Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateCauchyStrain, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Matrix B(5, 5);
    // clang-format off
    B <<=1.0, 2.0, 3.0, 4.0, 5.0,
         1.1, 1.2, 1.3, 1.4, 1.5,
         1.0, 0.9, 0.8, 0.7, 0.6,
         0.0, 1.0, 2.0, 3.0, 4.0,
         5.0, 4.0, 3.0, 2.0, 1.0;
    // clang-format on

    Vector displacements(5);
    displacements <<= 0.01, 0.02, 0.03, 0.04, 0.05;

    const auto strain = StressStrainUtilities::CalculateCauchyStrain(B, displacements);

    Vector expected_strain(5);
    expected_strain <<= 0.55, 0.205, 0.11, 0.4, 0.35;

    KRATOS_EXPECT_VECTOR_NEAR(strain, expected_strain, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculateStrains, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Matrix B(4, 4);
    // clang-format off
    B <<=1.0, 2.0, 3.0, 4.0,
         1.1, 1.2, 1.3, 1.4,
         1.0, 0.9, 0.8, 0.7,
         0.0, 1.0, 2.0, 3.0;
    // clang-format on

    Vector displacements(4);
    displacements <<= 0.01, 0.02, 0.03, 0.04;

    bool        use_hencky_strain = false;
    std::size_t voigt_size        = 4;

    std::vector<Matrix> Bs;
    Bs.push_back(B);
    Bs.push_back(B);

    std::vector<Matrix> deformation_gradients;
    Matrix              deformation_gradient = std::sqrt(3.) * IdentityMatrix(3);
    deformation_gradients.push_back(deformation_gradient);
    deformation_gradients.push_back(deformation_gradient);

    auto strains = StressStrainUtilities::CalculateStrains(deformation_gradients, Bs, displacements,
                                                           use_hencky_strain, voigt_size);
    Vector expected_strain(4);
    expected_strain <<= 0.3, 0.13, 0.08, 0.2;
    std::vector<Vector> expected_strains;
    expected_strains.push_back(expected_strain);
    expected_strains.push_back(expected_strain);

    for (size_t i = 0; i < strains.size(); ++i)
        KRATOS_EXPECT_VECTOR_NEAR(strains[i], expected_strains[i], Defaults::absolute_tolerance)

    use_hencky_strain = true;
    strains = StressStrainUtilities::CalculateStrains(deformation_gradients, Bs, displacements,
                                                      use_hencky_strain, voigt_size);

    expected_strain <<= 0.549306, 0.549306, 0.549306, 0;
    expected_strains.clear();
    expected_strains.push_back(expected_strain);
    expected_strains.push_back(expected_strain);

    for (size_t i = 0; i < strains.size(); ++i)
        KRATOS_EXPECT_VECTOR_NEAR(strains[i], expected_strains[i], 1.E-6)
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculatePrincipalStresses, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector cauchy_stresses = ZeroVector(6);
    cauchy_stresses <<= 80.0, 50.0, 20.0, 40.0, 35.0, 45.0;

    Vector principal_stresses;
    Matrix rotation_matrix;
    StressStrainUtilities::CalculatePrincipalStresses(cauchy_stresses, principal_stresses, rotation_matrix);

    Vector expected_solution_vector = ZeroVector(3);
    expected_solution_vector <<= 135.736961146391, 22.5224297324582, -8.25939087884923;
    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses, expected_solution_vector, Defaults::absolute_tolerance)

    Matrix expected_solution_matrix = ZeroMatrix(3, 3);
    expected_solution_matrix <<= 0.7304344631839778, -0.6094259674898185, -0.3083269127764114,
        0.5210391553001751, 0.7890961797761399, -0.3253389274384210, 0.4415695796102967,
        0.0769883706270019, 0.8939178357941993;
    KRATOS_EXPECT_MATRIX_NEAR(rotation_matrix, expected_solution_matrix, Defaults::absolute_tolerance)
}

} // namespace Kratos::Testing
