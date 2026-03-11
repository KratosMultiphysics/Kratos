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

#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "custom_utilities/ublas_utilities.h"
#include "tests/cpp_tests/test_utilities.h"
#include "utilities/math_utils.h"

using namespace Kratos;

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateTrace)
{
    const Vector eye = MathUtils<double>::StressTensorToVector(IdentityMatrix(3));
    EXPECT_DOUBLE_EQ(3.0, StressStrainUtilities::CalculateTrace(eye));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateMeanStress)
{
    const Vector eye = MathUtils<double>::StressTensorToVector(IdentityMatrix(3));
    EXPECT_DOUBLE_EQ(1.0, StressStrainUtilities::CalculateMeanStress(eye));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckGreenLagrangeStrainTensor)
{
    const Matrix eye =
        StressStrainUtilities::CalculateGreenLagrangeStrainTensor(std::sqrt(3.) * IdentityMatrix(3));
    KRATOS_EXPECT_MATRIX_NEAR(eye, IdentityMatrix(3), Defaults::absolute_tolerance)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateVonMisesStressHydrostatic)
{
    const auto stress_vector = UblasUtilities::CreateVector({-2.0, -2.0, -2.0, 0.0});
    EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateVonMisesStress(stress_vector));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateVonMisesStressPureShear)
{
    const auto stress_vector = UblasUtilities::CreateVector({0.0, 0.0, 0.0, 2.0});
    EXPECT_DOUBLE_EQ(2.0 * std::sqrt(3.0), StressStrainUtilities::CalculateVonMisesStress(stress_vector));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateLodeAngle)
{
    // Validate Triaxial Extension (TXE)
    auto stress_vector = UblasUtilities::CreateVector({0.5, -1.0, 0.5, 0.0});
    EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(30.0),
                            StressStrainUtilities::CalculateLodeAngle(stress_vector));
    // Validate Triaxial Compression (TXC)
    EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(-30.0),
                     StressStrainUtilities::CalculateLodeAngle(-1.0 * stress_vector));
    // Validate Shear (SHR)
    stress_vector = UblasUtilities::CreateVector({-1.0, 0.0, 1.0, 0.0});
    EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(0.0),
                            StressStrainUtilities::CalculateLodeAngle(stress_vector));

    // Regression tests with a small perturbation
    constexpr auto perturbation = 1.0e-8;
    stress_vector = UblasUtilities::CreateVector({0.5 + perturbation, -1.0, 0.5 + perturbation, 0.0});
    EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(30.0),
                            StressStrainUtilities::CalculateLodeAngle(stress_vector));
    EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(-30.0),
                            StressStrainUtilities::CalculateLodeAngle(-1.0 * stress_vector));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateMohrCoulombShearCapacityZeroStress)
{
    const Vector   stress_vector  = ZeroVector(4);
    constexpr auto cohesion       = 2.0;
    constexpr auto friction_angle = 0.0;
    EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateMohrCoulombShearCapacity(
                              stress_vector, cohesion, friction_angle));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateMohrCoulombShearCapacityZeroQMCResultsIn1)
{
    const Vector   stress_vector  = ZeroVector(4);
    constexpr auto cohesion       = 0.0;
    constexpr auto friction_angle = 0.0;
    EXPECT_DOUBLE_EQ(1.0, StressStrainUtilities::CalculateMohrCoulombShearCapacity(
                              stress_vector, cohesion, friction_angle));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateMohrCoulombShearCapacityHydrostatic)
{
    const auto     stress_vector = UblasUtilities::CreateVector({-2.0, -2.0, -2.0, 0.0});
    constexpr auto cohesion      = 0.0;
    constexpr auto friction_angle_in_radians = MathUtils<>::DegreesToRadians(90.0);
    EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateMohrCoulombShearCapacity(
                              stress_vector, cohesion, friction_angle_in_radians));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateMohrCoulombShearCapacityShearOnly)
{
    const auto     stress_vector  = UblasUtilities::CreateVector({-2.0, 0.0, 2.0, 0.0});
    constexpr auto cohesion       = 2.0;
    constexpr auto friction_angle = 0.0;
    EXPECT_DOUBLE_EQ(1.0, StressStrainUtilities::CalculateMohrCoulombShearCapacity(
                              stress_vector, cohesion, friction_angle));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateMohrCoulombShearCapacityThrowsWhenPhiIsOutOfBounds)
{
    const auto     stress_vector             = UblasUtilities::CreateVector({-2.0, 0.0, 2.0, 0.0});
    constexpr auto cohesion                  = 2.0;
    auto           friction_angle_in_radians = MathUtils<>::DegreesToRadians(-0.0001);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        StressStrainUtilities::CalculateMohrCoulombShearCapacity(stress_vector, cohesion, friction_angle_in_radians),
        "Friction angle must be in the range [0, 90] (degrees) : -0.0001");
    friction_angle_in_radians = MathUtils<>::DegreesToRadians(90.0001);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        StressStrainUtilities::CalculateMohrCoulombShearCapacity(stress_vector, cohesion, friction_angle_in_radians),
        "Friction angle must be in the range [0, 90] (degrees) : 90.0001");
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateMohrCoulombPressureCapacityZeroStress)
{
    Vector         stress_vector  = ZeroVector(4);
    constexpr auto cohesion       = 2.0;
    constexpr auto friction_angle = 0.0;
    EXPECT_DOUBLE_EQ(0.0, StressStrainUtilities::CalculateMohrCoulombPressureCapacity(
                              stress_vector, cohesion, friction_angle));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateMohrCoulombPressureCapacityShearOnly)
{
    const auto     stress_vector  = UblasUtilities::CreateVector({-2.0, 0.0, 2.0, 0.0});
    const auto     cohesion       = std::sqrt(2.0 * 2.0 * 3.0);
    constexpr auto friction_angle = MathUtils<>::DegreesToRadians(30.0);
    EXPECT_NEAR(3. * std::sin(friction_angle),
                StressStrainUtilities::CalculateMohrCoulombPressureCapacity(stress_vector, cohesion, friction_angle),
                Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateCauchyStrain)
{
    const auto B = UblasUtilities::CreateMatrix({{1.0, 2.0, 3.0, 4.0, 5.0},
                                                 {1.1, 1.2, 1.3, 1.4, 1.5},
                                                 {1.0, 0.9, 0.8, 0.7, 0.6},
                                                 {0.0, 1.0, 2.0, 3.0, 4.0},
                                                 {5.0, 4.0, 3.0, 2.0, 1.0}});

    const auto displacements = UblasUtilities::CreateVector({0.01, 0.02, 0.03, 0.04, 0.05});

    const auto strain = StressStrainUtilities::CalculateCauchyStrain(B, displacements);

    const auto expected_strain = UblasUtilities::CreateVector({0.55, 0.205, 0.11, 0.4, 0.35});

    KRATOS_EXPECT_VECTOR_NEAR(strain, expected_strain, Defaults::absolute_tolerance)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculateStrains)
{
    const auto B = UblasUtilities::CreateMatrix(
        {{1.0, 2.0, 3.0, 4.0}, {1.1, 1.2, 1.3, 1.4}, {1.0, 0.9, 0.8, 0.7}, {0.0, 1.0, 2.0, 3.0}});

    const auto displacements = UblasUtilities::CreateVector({0.01, 0.02, 0.03, 0.04});

    bool        use_hencky_strain = false;
    std::size_t voigt_size        = 4;

    std::vector<Matrix> Bs;
    Bs.push_back(B);
    Bs.push_back(B);

    std::vector<Matrix> deformation_gradients;
    Matrix              deformation_gradient = std::sqrt(3.0) * IdentityMatrix(3);
    deformation_gradients.push_back(deformation_gradient);
    deformation_gradients.push_back(deformation_gradient);

    auto strains = StressStrainUtilities::CalculateStrains(deformation_gradients, Bs, displacements,
                                                           use_hencky_strain, voigt_size);
    auto expected_strain = UblasUtilities::CreateVector({0.3, 0.13, 0.08, 0.2});
    std::vector<Vector> expected_strains;
    expected_strains.push_back(expected_strain);
    expected_strains.push_back(expected_strain);

    for (size_t i = 0; i < strains.size(); ++i)
        KRATOS_EXPECT_VECTOR_NEAR(strains[i], expected_strains[i], Defaults::absolute_tolerance)

    use_hencky_strain = true;
    strains = StressStrainUtilities::CalculateStrains(deformation_gradients, Bs, displacements,
                                                      use_hencky_strain, voigt_size);

    expected_strain = UblasUtilities::CreateVector({0.549306, 0.549306, 0.549306, 0});
    expected_strains.clear();
    expected_strains.push_back(expected_strain);
    expected_strains.push_back(expected_strain);

    for (size_t i = 0; i < strains.size(); ++i)
        KRATOS_EXPECT_VECTOR_NEAR(strains[i], expected_strains[i], 1.E-6)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckCalculatePrincipalStresses)
{
    // Arrange
    const auto cauchy_stresses = UblasUtilities::CreateVector({80.0, 50.0, 20.0, 40.0, 35.0, 45.0});

    // Act
    Vector actual_principal_stresses;
    Matrix actual_rotation_matrix;
    StressStrainUtilities::CalculatePrincipalStresses(cauchy_stresses, actual_principal_stresses,
                                                      actual_rotation_matrix);

    // Assert
    const auto expected_principal_stresses =
        UblasUtilities::CreateVector({135.736961146391, 22.5224297324582, -8.25939087884923});
    KRATOS_EXPECT_VECTOR_NEAR(actual_principal_stresses, expected_principal_stresses, Defaults::absolute_tolerance)

    const auto expected_rotation_matrix =
        UblasUtilities::CreateMatrix({{0.7304344631839778, -0.6094259674898185, -0.3083269127764114},
                                      {0.5210391553001751, 0.7890961797761399, -0.3253389274384210},
                                      {0.4415695796102967, 0.0769883706270019, 0.8939178357941993}});
    KRATOS_EXPECT_MATRIX_NEAR(actual_rotation_matrix, expected_rotation_matrix, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CheckCalculatePrincipalStressesAndRotationMatrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto cauchy_stresses = UblasUtilities::CreateVector({80.0, 50.0, 20.0, 40.0, 35.0, 45.0});

    // Act
    const auto& [actual_principal_stresses, actual_rotation_matrix] =
        StressStrainUtilities::CalculatePrincipalStressesAndRotationMatrix(cauchy_stresses);

    // Assert
    const auto expected_principal_stresses =
        Geo::PrincipalStresses{135.736961146391, 22.5224297324582, -8.25939087884923};
    KRATOS_EXPECT_VECTOR_NEAR(actual_principal_stresses.Values(),
                              expected_principal_stresses.Values(), Defaults::absolute_tolerance)

    const auto expected_rotation_matrix =
        UblasUtilities::CreateMatrix({{0.7304344631839778, -0.6094259674898185, -0.3083269127764114},
                                      {0.5210391553001751, 0.7890961797761399, -0.3253389274384210},
                                      {0.4415695796102967, 0.0769883706270019, 0.8939178357941993}});
    KRATOS_EXPECT_MATRIX_NEAR(actual_rotation_matrix, expected_rotation_matrix, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CheckTransformPrincipalStressesToPandQ, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto principal_stresses = Geo::PrincipalStresses{30.0, 20.0, 10.0};
    const auto p_q = StressStrainUtilities::TransformPrincipalStressesToPandQ(principal_stresses);
    const auto expected_solution_pq = Geo::PQ{20.0, std::sqrt(300.0)};
    KRATOS_EXPECT_VECTOR_NEAR(p_q.Values(), expected_solution_pq.Values(), Defaults::absolute_tolerance)
}

} // namespace Kratos::Testing
