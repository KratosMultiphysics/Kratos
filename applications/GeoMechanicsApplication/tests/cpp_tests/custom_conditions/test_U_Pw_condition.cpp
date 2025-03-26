// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

// Project includes
#include "custom_conditions/U_Pw_condition.hpp"

#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateLeftHandSideUPwCondition_ReturnsEmptyMatrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto p_dummy_geometry = std::make_shared<Geometry<Node>>();
    auto cond             = UPwCondition<2, 2>(1, p_dummy_geometry, nullptr);

    Matrix left_hand_side_matrix = ZeroMatrix(6, 6);
    cond.CalculateLeftHandSide(left_hand_side_matrix, ProcessInfo{});

    Matrix expected_matrix = ZeroMatrix(0, 0);
    KRATOS_CHECK_MATRIX_NEAR(left_hand_side_matrix, expected_matrix, 1.0e-6);
}

} // namespace Kratos::Testing
