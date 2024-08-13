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
#include "custom_elements/U_Pw_condition.hpp"

#include "geo_mechanics_fast_suite.h"

namespace Kratos::Testing {


        /// <summary>
		/// Tests the calculation of the left hand side matrix of the default UPwCondition, which should be an empty matrix.
        /// </summary>
        /// <param name=""></param>
        /// <param name=""></param>
        KRATOS_TEST_CASE_IN_SUITE(CalculateLeftHandSideUPwCondition, KratosGeoMechanicsFastSuite)
        {

            Model current_model;
            auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
			const auto& r_process_info = r_model_part.GetProcessInfo();
  
			auto p_cond = std::make_shared<UPwCondition<2, 2>>(1, nullptr);

			// calculate left hand side matrix
			Matrix left_hand_side_matrix = ZeroMatrix(6, 6);
			p_cond->CalculateLeftHandSide(left_hand_side_matrix, r_process_info);

			// set expected_results
			Matrix expected_matrix = ZeroMatrix(0, 0);

			// compare results
			KRATOS_CHECK_MATRIX_NEAR(left_hand_side_matrix, expected_matrix, 1.0e-6);
        }

}