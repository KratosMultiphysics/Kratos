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
namespace Kratos::Testing {


        KRATOS_TEST_CASE_IN_SUITE(CalculateLeftHandSideUPwCondition_ReturnsEmptyMatrix, KratosGeoMechanicsFastSuite)
        {

            Model current_model;
            auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);	

            // create geometry as UPwCondition needs a geometry to be initialized
            auto node = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            std::vector< ModelPart::IndexType> element_nodes{ 1};
            auto p_geometry = r_model_part.CreateNewGeometry("Point2D", element_nodes);
            
            auto cond = UPwCondition<2, 2>(1, p_geometry, nullptr);
            
            Matrix left_hand_side_matrix = ZeroMatrix(6, 6);
            cond.CalculateLeftHandSide(left_hand_side_matrix, ProcessInfo{});
            
            Matrix expected_matrix = ZeroMatrix(0, 0);
            KRATOS_CHECK_MATRIX_NEAR(left_hand_side_matrix, expected_matrix, 1.0e-6);
        }

}