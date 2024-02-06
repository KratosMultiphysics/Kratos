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

#include "custom_elements/axisymmetric_stress_state_policy.h"
#include "custom_elements/stress_state_policy.h"
#include "includes/checks.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateBMatrixGivesExpectedMatrix, KratosGeoMechanicsFastSuite)
{
    std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<AxisymmetricStressStatePolicy>();

    const Matrix expected_matrix   = ZeroMatrix(0, 0);
    const Matrix calculated_matrix = p_stress_state_policy->CalculateBMatrix(Matrix(), Vector());

    KRATOS_CHECK_MATRIX_NEAR(expected_matrix, calculated_matrix, 1e-12)
}

} // namespace Kratos::Testing
