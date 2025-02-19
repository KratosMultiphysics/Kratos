// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Marjan Fathian
//

#include "containers/model.h"
#include "custom_elements/three_dimensional_stress_state.h"
#include "custom_utilities/registration_utilities.h"
#include "includes/checks.h"
#include "includes/serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <sstream>
#include <string>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionalStressState_CalculateBMatrixReturnsCorrectResults,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto p_stress_state_policy = std::make_unique<ThreeDimensionalStressState>();

    Vector Np(4);
    Np <<= 1.0, 2.0, 3.0, 4.0;

    // clang-format off
    Matrix GradNpT(4, 3);
    GradNpT <<= 1.0,  2.0,  3.0,
                4.0,  5.0,  6.0,
                7.0,  8.0,  9.0,
                10.0, 11.0, 12.0;
    // clang-format on

    const auto calculated_matrix = p_stress_state_policy->CalculateBMatrix(
        GradNpT, Np, ModelSetupUtilities::Create3D4NTetrahedraGeometry());

    // clang-format off
    Matrix expected_matrix(6, 12);
    expected_matrix <<= // This row (INDEX_3D_XX) contains the first column of GradNpT as INDEX_X
                        1,  0,  0,  4,  0,  0,  7,  0,  0,  10, 0,  0,
                        // This row (INDEX_3D_YY) contains the second column of GradNpT as INDEX_Y
                        0,  2,  0,  0,  5,  0,  0,  8,  0,  0,  11, 0,
                        // This row (INDEX_3D_ZZ) contains the third column of GradNpT as INDEX_Z
                        0,  0,  3,  0,  0,  6,  0,  0,  9,  0,  0,  12,
                        // This row (INDEX_3D_XY) contains the first and second columns of GradNpT as INDEX_Y and INDEX_X respectively
                        2,  1,  0,  5,  4,  0,  8,  7,  0,  11, 10, 0,
                        // This row (INDEX_3D_YZ) contains the second and third column of GradNpT as INDEX_Z and INDEX_Y respectively
                        0,  3,  2,  0,  6,  5,  0,  9,  8,  0,  12, 11,
                        // This row (INDEX_3D_XZ) contains the first and third column of GradNpT as INDEX_Z and INDEX_X respectively
                        3,  0,  1,  6,  0,  4,  9,  0,  7,  12, 0,  10;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(calculated_matrix, expected_matrix, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionalStressState_ReturnsCorrectIntegrationCoefficient,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto p_stress_state_policy = std::make_unique<ThreeDimensionalStressState>();

    // The shape function values for this integration point are 0.2, 0.5 and 0.3 for nodes 1, 2 and 3 respectively
    Geometry<Node>::IntegrationPointType integration_point(0.5, 0.3, 0.0, 0.5);

    const auto detJ                   = 2.0;
    const auto calculated_coefficient = p_stress_state_policy->CalculateIntegrationCoefficient(
        integration_point, detJ, ModelSetupUtilities::Create3D4NTetrahedraGeometry());

    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) = 1.0
    KRATOS_EXPECT_NEAR(calculated_coefficient, 1.0, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionalStressState_GivesCorrectClone, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<ThreeDimensionalStressState>();
    const auto p_clone = p_stress_state_policy->Clone();

    KRATOS_EXPECT_NE(dynamic_cast<ThreeDimensionalStressState*>(p_clone.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionalStressState_CalculateGreenLagrangeStrainReturnsCorrectResults,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto p_stress_state_policy = std::make_unique<ThreeDimensionalStressState>();

    // clang-format off
    Matrix deformation_gradient = ZeroMatrix(3, 3);
    deformation_gradient <<= 1.0, 2.0, 3.0,
                             4.0, 5.0, 6.0,
                             7.0, 8.0, 9.0;
    // clang-format on

    // The expected strain is calculated as follows:
    // 0.5 * (F^T * F - I) and then converted to a vector
    Vector expected_vector = ZeroVector(6);
    expected_vector <<= 32.5, 46, 62.5, 78, 108, 90;

    const auto calculated_vector = p_stress_state_policy->CalculateGreenLagrangeStrain(deformation_gradient);

    KRATOS_CHECK_VECTOR_NEAR(expected_vector, calculated_vector, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionalStressState_GivesCorrectVoigtVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<ThreeDimensionalStressState>();
    Vector voigt_vector = p_stress_state_policy->GetVoigtVector();

    Vector expected_voigt_vector(6);
    expected_voigt_vector <<= 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(voigt_vector, expected_voigt_vector, 1.E-10)
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionalStressState_GivesCorrectVoigtSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<ThreeDimensionalStressState>();
    KRATOS_EXPECT_EQ(p_stress_state_policy->GetVoigtSize(), VOIGT_SIZE_3D);
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionalStressState_GivesCorrectStressTensorSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<ThreeDimensionalStressState>();
    KRATOS_EXPECT_EQ(p_stress_state_policy->GetStressTensorSize(), STRESS_TENSOR_SIZE_3D);
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionalStressState_CanBeSavedAndLoadedThroughInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    RegistrationUtilities::RegisterStressStatePolicies();
    const auto p_policy =
        std::unique_ptr<StressStatePolicy>{std::make_unique<ThreeDimensionalStressState>()};
    auto serializer = Serializer{new std::stringstream{}};

    // Act
    serializer.save("test_tag"s, p_policy);
    auto p_loaded_policy = std::unique_ptr<StressStatePolicy>{};
    serializer.load("test_tag"s, p_loaded_policy);

    // Assert
    ASSERT_NE(p_loaded_policy, nullptr);
    KRATOS_EXPECT_EQ(p_loaded_policy->GetVoigtSize(), VOIGT_SIZE_3D);
    auto expected_voigt_vector = Vector{6};
    expected_voigt_vector <<= 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(p_loaded_policy->GetVoigtVector(), expected_voigt_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
