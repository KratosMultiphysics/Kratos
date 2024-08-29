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

#include "custom_elements/interface_stress_state.h"
#include "geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CanCreateInterfaceStressState, KratosGeoMechanicsFastSuite)
{
    auto p_stress_state_policy = std::make_unique<InterfaceStressState>();

    KRATOS_EXPECT_NE(p_stress_state_policy, nullptr);
}

} // namespace Kratos::Testing

// VoigtVector -> [1,0]
// B matrix -> (Voigt size * n u dof) size

// Bmatrix, op de eerste rij, normaal, komt node n1(4 - node 1), n2(node 5 - node 2), n3(node 6 - node 3)
