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

#include "geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

/// <summary>
/// Tests the calculation of the left hand side matrix of the default UPwCondition, which should be an empty matrix.
/// </summary>
/// <param name=""></param>
/// <param name=""></param>
KRATOS_TEST_CASE_IN_SUITE(CalculateLeftHandSideUPwCondition, KratosGeoMechanicsFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

    // create geometry as UPwCondition needs a geometry to be initialized
    auto                              node = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    std::vector<ModelPart::IndexType> element_nodes{1};
    auto p_geometry = r_model_part.CreateNewGeometry("Point2D", element_nodes);

    // create UPwCondition
    auto cond = UPwCondition<2, 2>(1, p_geometry, nullptr);

    // calculate left hand side matrix
    Matrix left_hand_side_matrix = ZeroMatrix(6, 6);
    cond.CalculateLeftHandSide(left_hand_side_matrix, ProcessInfo{});

    // set expected_results
    Matrix expected_matrix = ZeroMatrix(0, 0);

    // compare results
    KRATOS_CHECK_MATRIX_NEAR(left_hand_side_matrix, expected_matrix, 1.0e-6);
}

} // namespace Kratos::Testing
