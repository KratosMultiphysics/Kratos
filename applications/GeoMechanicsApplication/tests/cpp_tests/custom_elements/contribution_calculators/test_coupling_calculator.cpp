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

#include "custom_elements/contribution_calculators/up_coupling_calculator.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

namespace Kratos::Testing
{
    TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestCanCreateUPCouplingCalculator)
    {
        constexpr auto number_of_u_dof = 4;
        constexpr auto number_of_pw_dof = 2;
        UPCouplingCalculator<number_of_u_dof, number_of_pw_dof> coupling_calculator;
    }
}
