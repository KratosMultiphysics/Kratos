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

    TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestCanCreateUpCouplingInputProvider)
    {
        constexpr auto number_of_u_dof = 4;
        constexpr auto number_of_pw_dof = 2;

        auto get_b_matrices = []()
        {
            return std::vector<Matrix>{};
        };
        auto get_integration_coefficients = []()
        {
            return std::vector<double>{};
        };
        const auto np_container = Matrix{};
        auto get_np_container = [&np_container]() -> const Matrix& { return np_container; };
        auto get_biot_coefficients = []()
        {
            return std::vector<double>{};
        };
        auto get_bishop_coefficients = []()
        {
            return std::vector<double>{};
        };
        auto get_voigt_vector = []()
        {
            return Vector{};
        };


        UPCouplingCalculator<number_of_u_dof, number_of_pw_dof>::InputProvider input_provider(
            get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients, get_biot_coefficients,
            get_bishop_coefficients);
    }
}
