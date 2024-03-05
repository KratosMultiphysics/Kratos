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
//

#include "containers/model.h"
#include "custom_elements/plane_strain_stress_state.h"
#include "includes/checks.h"
#include "testing/testing.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateBMatrixGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    auto p_stress_state_policy = std::make_unique<PlaneStrainStressState>();

    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    Vector Np(3);
    Np <<= 1.0, 2.0, 3.0;

    // clang-format off
    Matrix GradNpT(3, 2);
    GradNpT <<= 1.0, 2.0,
                3.0, 4.0,
                5.0, 6.0;
    // clang-format on

    const auto calculated_matrix =
        p_stress_state_policy->CalculateBMatrix(GradNpT, Np, r_model_part.GetElement(1).GetGeometry());

    // clang-format off
    Matrix expected_matrix(4, 6);
    expected_matrix <<= 1  ,0  ,3  ,0  ,5  ,0, // This row contains the first column of GradNpT on columns 1, 3 and 5
                        0  ,2  ,0  ,4  ,0  ,6, // This row contains the second column of GradNpT on columns 2, 4 and 6
                        0  ,0  ,0  ,0  ,0  ,0, // This row is not set, so remains 0.
                        2  ,1  ,4  ,3  ,6  ,5; // This row contains the first and second columns of GradNpT, swapping x and y
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(calculated_matrix, expected_matrix, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(PlaneStrainStressState_ReturnsCorrectIntegrationCoefficient, KratosGeoMechanicsFastSuite)
{
    const auto p_stress_state_policy = std::make_unique<PlaneStrainStressState>();

    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    // The shape function values for this integration point are 0.2, 0.5 and 0.3 for nodes 1, 2 and 3 respectively
    Geometry<Node>::IntegrationPointType integration_point(0.5, 0.3, 0.0, 0.5);

    const auto detJ                   = 2.0;
    const auto calculated_coefficient = p_stress_state_policy->CalculateIntegrationCoefficient(
        integration_point, detJ, r_model_part.GetElement(1).GetGeometry());

    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) = 1.0
    KRATOS_EXPECT_NEAR(calculated_coefficient, 1.0, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(PlaneStrainStressState_ReturnsCorrectGreenLagrangeStrain, KratosGeoMechanicsFastSuite)
{
    const auto p_stress_state_policy = std::make_unique<PlaneStrainStressState>();

    // clang-format off
    Matrix deformation_gradient = ZeroMatrix(2, 2);
    deformation_gradient <<= 1.0, 2.0,
                             3.0, 4.0;
    // clang-format on

    const auto calculated_strain = p_stress_state_policy->CalculateGreenLagrangeStrain(deformation_gradient);

    Vector expected_vector = ZeroVector(4);
    expected_vector <<= 4.5, 9.5, 0, 14;

    // The expected strain is calculated as follows:
    // 0.5 * (F^T * F - I) and then converted to a vector
    KRATOS_CHECK_VECTOR_NEAR(calculated_strain, expected_vector, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(PlaneStrainStressState_GivesCorrectClone, KratosGeoMechanicsFastSuite)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<PlaneStrainStressState>();
    const auto p_clone = p_stress_state_policy->Clone();

    KRATOS_EXPECT_NE(dynamic_cast<PlaneStrainStressState*>(p_clone.get()), nullptr);
}

} // namespace Kratos::Testing