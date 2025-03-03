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
    KRATOS_EXPECT_MATRIX_NEAR(eye, IdentityMatrix(3), 1e-12)
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
    // Validate Triaxial Extension (TXE)
    stress_vector <<= 0.5, -1.0, 0.5, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(30.),
                            StressStrainUtilities::CalculateLodeAngle(stress_vector));
    // Validate Triaxial Compression (TXC)
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(-30.),
                            StressStrainUtilities::CalculateLodeAngle(-1. * stress_vector));
    // Validate Shear (SHR)
    stress_vector <<= -1.0, 0.0, 1.0, 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(MathUtils<>::DegreesToRadians(0.),
                            StressStrainUtilities::CalculateLodeAngle(stress_vector));
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
                       1.E-10);
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

    KRATOS_EXPECT_VECTOR_NEAR(strain, expected_strain, 1.E-10)
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
        KRATOS_EXPECT_VECTOR_NEAR(strains[i], expected_strains[i], 1.E-10)

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

} // namespace Kratos::Testing
