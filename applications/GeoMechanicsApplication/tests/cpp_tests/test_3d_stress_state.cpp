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
#include "custom_elements/three_dimension_stress_state.h"
#include "includes/checks.h"
#include "testing/testing.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"
#include <boost/numeric/ublas/assignment.hpp>


using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionStressState_CalculateBMatrixReturnsCorrectResults, KratosGeoMechanicsFastSuite)
{
    auto stress_state = std::make_unique<ThreeDimensionStressState>();

    Model      model;
    ModelPart& model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model);

    Vector Np(4);
    Np <<= 1.0, 2.0, 3.0, 4.0;

    // clang-format off
    Matrix GradNpT(4, 3);
    GradNpT <<= 1.0,  2.0,  3.0,
                4.0,  5.0,  6.0,
                7.0,  8.0,  9.0,
                10.0, 11.0, 12.0;
    // clang-format on

    const Matrix calculated_matrix =
            stress_state->CalculateBMatrix(GradNpT, Np, model_part.GetElement(1).GetGeometry());

    // clang-format off
    Matrix expected_matrix(6, 12);
    expected_matrix <<= 1,  0,  0,  4,  0,  0,  7,  0,  0,  10, 0,  0,  // This row (INDEX_3D_XX) contains the first column of GradNpT as INDEX_X
                        0,  2,  0,  0,  5,  0,  0,  8,  0,  0,  11, 0,  // This row (INDEX_3D_YY) contains the second column of GradNpT as INDEX_Y
                        0,  0,  3,  0,  0,  6,  0,  0,  9,  0,  0,  12, // This row (INDEX_3D_ZZ) contains the third column of GradNpT as INDEX_Z
                        2,  1,  0,  5,  4,  0,  8,  7,  0,  11, 10, 0,  // This row (INDEX_3D_XY) contains the first and second columns of GradNpT as INDEX_Y and INDEX_X respectively
                        0,  3,  2,  0,  6,  5,  0,  9,  8,  0,  12, 11, // This row (INDEX_3D_YZ) contains the second and third column of GradNpT as INDEX_Z and INDEX_Y respectively
                        3,  0,  1,  6,  0,  4,  9,  0,  7,  12, 0,  10; // This row (INDEX_3D_XZ) contains the first and third column of GradNpT as INDEX_Z and INDEX_X respectively
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(calculated_matrix, expected_matrix, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionStressState_ReturnsCorrectIntegrationCoefficient, KratosGeoMechanicsFastSuite)
{
    const auto p_stress_state_policy = std::make_unique<ThreeDimensionStressState>();

    Model      model;
    ModelPart& model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    // The shape function values for this integration point are 0.2, 0.5 and 0.3 for nodes 1, 2 and 3 respectively
    Geometry<Node>::IntegrationPointType integration_point(0.5, 0.3, 0.0, 0.5);

    const double detJ                   = 2.0;
    const double calculated_coefficient = p_stress_state_policy->CalculateIntegrationCoefficient(
            integration_point, detJ, model_part.GetElement(1).GetGeometry());

    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) = 1.0
    KRATOS_EXPECT_NEAR(calculated_coefficient, 1.0, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionStressState_GivesCorrectClone, KratosGeoMechanicsFastSuite)
{
    const auto p_stress_state_policy = std::make_unique<ThreeDimensionStressState>();
    const auto p_clone = p_stress_state_policy->Clone();

    KRATOS_EXPECT_NE(dynamic_cast<ThreeDimensionStressState*>(p_clone.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(ThreeDimensionStressState_CalculateGreenLagrangeStrainGivesExpectedVector, KratosGeoMechanicsFastSuite)
{
    const auto p_stress_state_policy = std::make_unique<ThreeDimensionStressState>();

    const Vector expected_vector   = ZeroVector(0);
    const Vector calculated_vector = p_stress_state_policy->CalculateGreenLagrangeStrain(Matrix());

    KRATOS_CHECK_VECTOR_NEAR(expected_vector, calculated_vector, 1e-12)
}

}