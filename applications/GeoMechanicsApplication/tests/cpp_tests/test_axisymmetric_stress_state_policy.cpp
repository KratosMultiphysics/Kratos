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
#include "geo_mechanics_application_constants.h"
#include "includes/checks.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateBMatrixGivesCorrectMatrixSize, KratosGeoMechanicsFastSuite)
{
    std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<AxisymmetricStressStatePolicy>();

    const SizeType voigt_size = VOIGT_SIZE_2D_PLANE_STRAIN;
    const SizeType number_of_nodes = 3;
    const SizeType dimension       = 2;

    const Matrix expected_matrix   = ZeroMatrix(voigt_size, number_of_nodes * dimension);
    const Matrix calculated_matrix = p_stress_state_policy->CalculateBMatrix(Matrix(), Vector(), number_of_nodes);

    KRATOS_EXPECT_EQ(calculated_matrix.size1(), voigt_size);
    KRATOS_EXPECT_EQ(calculated_matrix.size2(), number_of_nodes * dimension);

}

} // namespace Kratos::Testing
