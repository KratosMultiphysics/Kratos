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
#include "geometries/line_2d_2.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

using namespace Kratos;

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateLeftHandSideUPwCondition_ReturnsEmptyMatrix)
{
    auto p_line_geometry = std::make_shared<Line2D2<Node>>(make_intrusive<Node>(1, 0.0, 0.0, 0.0),
                                                           make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    auto cond            = UPwCondition<2, 2>(1, p_line_geometry, nullptr);

    Matrix left_hand_side_matrix = ZeroMatrix(6, 6);
    cond.CalculateLeftHandSide(left_hand_side_matrix, ProcessInfo{});

    Matrix expected_matrix = ZeroMatrix(0, 0);
    KRATOS_CHECK_MATRIX_NEAR(left_hand_side_matrix, expected_matrix, 1.0e-6);
}

} // namespace Kratos::Testing
