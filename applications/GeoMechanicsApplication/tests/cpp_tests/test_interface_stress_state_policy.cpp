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
#include "custom_geometries/line_interface_geometry.h"
#include "geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CanCreateInterfaceStressState, KratosGeoMechanicsFastSuite)
{
    auto p_stress_state_policy = std::make_unique<InterfaceStressState>();

    KRATOS_EXPECT_NE(p_stress_state_policy, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_CloneCreatesCorrectInstance, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    std::unique_ptr<StressStatePolicy> p_stress_state_policy = std::make_unique<InterfaceStressState>();

    KRATOS_EXPECT_NE(dynamic_cast<InterfaceStressState*>(p_stress_state_policy->Clone().get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_ReturnsEmptyBMatrixWhenInputtingEmptyShapeFunctionValues,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto stress_state_policy = InterfaceStressState{};

    const auto b_matrix = stress_state_policy.CalculateBMatrix({}, {}, {});

    KRATOS_EXPECT_EQ(b_matrix.size1(), 0);
    KRATOS_EXPECT_EQ(b_matrix.size2(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_ReturnsExpectedVoigtSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto stress_state_policy = InterfaceStressState{};

    KRATOS_EXPECT_EQ(stress_state_policy.GetVoigtSize(), 2);
}

} // namespace Kratos::Testing

// VoigtVector -> [1,0]
// B matrix -> (Voigt size * n u dof) size

// Bmatrix, op de eerste rij, normaal, komt node n1(4 - node 1), n2(node 5 - node 2), n3(node 6 - node 3)
