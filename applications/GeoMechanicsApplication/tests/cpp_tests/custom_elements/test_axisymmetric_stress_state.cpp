// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Marjan Fathian
//

#include "containers/model.h"
#include "custom_elements/axisymmetric_stress_state.h"
#include "custom_elements/stress_state_policy.h"
#include "custom_utilities/registration_utilities.h"
#include "geometries/geometry.h"
#include "includes/checks.h"
#include "includes/stream_serializer.h"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <string>
#include <type_traits>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateBMatrixWithValidGeometryReturnsCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<AxisymmetricStressState>();

    Vector Np(3);
    Np <<= 1.0, 2.0, 3.0;

    // clang-format off
    Matrix GradNpT(3, 2);
    GradNpT <<= 1.0, 2.0,
                3.0, 4.0,
                5.0, 6.0;
    // clang-format on

    const Matrix calculated_matrix = p_stress_state_policy->CalculateBMatrix(
        GradNpT, Np, ModelSetupUtilities::Create2D3NTriangleGeometry());

    // clang-format off
    Matrix expected_matrix(4, 6);
    expected_matrix <<= 1   ,0  ,3   ,0 ,5   ,0, // This row contains the first column of GradNpT on columns 1, 3 and 5
                        0   ,2  ,0   ,4 ,0   ,6, // This row contains the second column of GradNpT on columns 2, 4 and 6
                        0.2 ,0  ,0.4 ,0 ,0.6 ,0, // This row contains Np/radius on columns 1, 3 and 5, where radius = 5
                        2   ,1  ,4   ,3 ,6   ,5; // This row contains the first and second columns of GradNpT, swapping x and y
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(calculated_matrix, expected_matrix, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(AxisymmetricStressState_CannotBeCopiedButItCanBeMoved, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    EXPECT_FALSE(std::is_copy_constructible_v<AxisymmetricStressState>);
    EXPECT_FALSE(std::is_copy_assignable_v<AxisymmetricStressState>);
    EXPECT_TRUE(std::is_move_constructible_v<AxisymmetricStressState>);
    EXPECT_TRUE(std::is_move_assignable_v<AxisymmetricStressState>);
}

KRATOS_TEST_CASE_IN_SUITE(TestCloneReturnsCorrectType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<AxisymmetricStressState>();
    const auto p_cloned_stress_state_policy = p_stress_state_policy->Clone();

    KRATOS_EXPECT_NE(dynamic_cast<AxisymmetricStressState*>(p_cloned_stress_state_policy.get()), nullptr);
    KRATOS_EXPECT_NE(p_cloned_stress_state_policy.get(), p_stress_state_policy.get());
}

KRATOS_TEST_CASE_IN_SUITE(TestCalculateGreenLagrangeStrainThrows, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<AxisymmetricStressState>();

    // Note: avoid a warning triggered by the `[[nodiscard]]` attribute of the `CalculateGreenLagrangeStrain()`
    // member function by assigning the return value to a dummy variable. In turn, the dummy
    // variable needs to be marked `[[maybe_unused]]` to avoid a warning about an unused variable.
    KRATOS_EXPECT_EXCEPTION_IS_THROWN([[maybe_unused]] const auto strain_vector =
                                          p_stress_state_policy->CalculateGreenLagrangeStrain(Matrix()),
                                      "The calculation of Green Lagrange strain is "
                                      "not implemented for axisymmetric configurations.")
}

KRATOS_TEST_CASE_IN_SUITE(AxisymmetricStressState_GivesCorrectVoigtVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<AxisymmetricStressState>();
    Vector voigt_vector = p_stress_state_policy->GetVoigtVector();

    Vector expected_voigt_vector(4);
    expected_voigt_vector <<= 1.0, 1.0, 1.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(voigt_vector, expected_voigt_vector, 1.E-10)
}

KRATOS_TEST_CASE_IN_SUITE(AxisymmetricStressState_GivesCorrectVoigtSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<AxisymmetricStressState>();
    KRATOS_EXPECT_EQ(p_stress_state_policy->GetVoigtSize(), VOIGT_SIZE_2D_PLANE_STRAIN);
}

KRATOS_TEST_CASE_IN_SUITE(AxisymmetricStressState_GivesCorrectStressTensorSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<AxisymmetricStressState>();
    KRATOS_EXPECT_EQ(p_stress_state_policy->GetStressTensorSize(), STRESS_TENSOR_SIZE_2D);
}

KRATOS_TEST_CASE_IN_SUITE(AxisymmetricStressState_CanBeSavedAndLoadedThroughInterface, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("AxisymmetricStressState"s, AxisymmetricStressState{})};
    const auto p_policy = std::unique_ptr<StressStatePolicy>{std::make_unique<AxisymmetricStressState>()};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_policy);
    auto p_loaded_policy = std::unique_ptr<StressStatePolicy>{};
    serializer.load("test_tag"s, p_loaded_policy);

    // Assert
    ASSERT_NE(p_loaded_policy, nullptr);
    KRATOS_EXPECT_EQ(p_loaded_policy->GetVoigtSize(), VOIGT_SIZE_2D_AXISYMMETRIC);
    auto expected_voigt_vector = Vector{4};
    expected_voigt_vector <<= 1.0, 1.0, 1.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(p_loaded_policy->GetVoigtVector(), expected_voigt_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
