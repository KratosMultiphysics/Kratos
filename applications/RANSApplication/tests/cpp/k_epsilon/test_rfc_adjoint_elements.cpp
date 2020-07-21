//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes
#include <functional>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"
#include "test_adjoint_utilities.h"

namespace Kratos
{
namespace Testing
{
// k-epsilon k element tests
KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_CalculateFirstDerivativesLHS, KratosRansFastSuite)
{
    const auto& calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, const ProcessInfo& rProcessInfo) {
            rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansKEpsilonKRFC2D3N", "RansKEpsilonKAdjointRFC2D3N",
        TURBULENT_KINETIC_ENERGY, calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_Calculate_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                          KratosRansFastSuite)
{
    const auto& calculate_sensitivity_matrix = [](Matrix& rOutput, ElementType& rElement,
                                                  const ProcessInfo& rProcessInfo) {
        rElement.Calculate(TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansKEpsilonKRFC2D3N", "RansKEpsilonKAdjointRFC2D3N",
        TURBULENT_ENERGY_DISSIPATION_RATE, calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_Calculate_VELOCITY_PARTIAL_DERIVATIVE, KratosRansFastSuite)
{
    const auto& calculate_sensitivity_matrix = [](Matrix& rOutput, ElementType& rElement,
                                                  const ProcessInfo& rProcessInfo) {
        rElement.Calculate(VELOCITY_PARTIAL_DERIVATIVE, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansKEpsilonKRFC2D3N", "RansKEpsilonKAdjointRFC2D3N", VELOCITY,
        calculate_sensitivity_matrix, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_CalculateSensitivityMatrix, KratosRansFastSuite)
{
    const auto& calculate_sensitivity_matrix = [](Matrix& rOutput, ElementType& rElement,
                                                  const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansKEpsilonKRFC2D3N", "RansKEpsilonKAdjointRFC2D3N",
        SHAPE_SENSITIVITY, calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_CalculateSecondDerivativesLHS, KratosRansFastSuite)
{
    const auto& calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, const ProcessInfo& rProcessInfo) {
            rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
            const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
            noalias(rOutput) = rOutput * (1.0 - bossak_alpha);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansKEpsilonKRFC2D3N", "RansKEpsilonKAdjointRFC2D3N",
        TURBULENT_KINETIC_ENERGY_RATE, calculate_sensitivity_matrix, 1e-7, 1e-5);
}

// k-epsilon epsilon element tests
KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_CalculateFirstDerivativesLHS, KratosRansFastSuite)
{
    const auto& calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, const ProcessInfo& rProcessInfo) {
            rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansKEpsilonEpsilonRFC2D3N", "RansKEpsilonEpsilonAdjointRFC2D3N",
        TURBULENT_ENERGY_DISSIPATION_RATE, calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_Calculate_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                          KratosRansFastSuite)
{
    const auto& calculate_sensitivity_matrix = [](Matrix& rOutput, ElementType& rElement,
                                                  const ProcessInfo& rProcessInfo) {
        rElement.Calculate(TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansKEpsilonEpsilonRFC2D3N", "RansKEpsilonEpsilonAdjointRFC2D3N",
        TURBULENT_KINETIC_ENERGY, calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_Calculate_VELOCITY_PARTIAL_DERIVATIVE,
                          KratosRansFastSuite)
{
    const auto& calculate_sensitivity_matrix = [](Matrix& rOutput, ElementType& rElement,
                                                  const ProcessInfo& rProcessInfo) {
        rElement.Calculate(VELOCITY_PARTIAL_DERIVATIVE, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansKEpsilonEpsilonRFC2D3N", "RansKEpsilonEpsilonAdjointRFC2D3N",
        VELOCITY, calculate_sensitivity_matrix, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_CalculateSensitivityMatrix, KratosRansFastSuite)
{
    const auto& calculate_sensitivity_matrix = [](Matrix& rOutput, ElementType& rElement,
                                                  const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansKEpsilonEpsilonRFC2D3N", "RansKEpsilonEpsilonAdjointRFC2D3N",
        SHAPE_SENSITIVITY, calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_CalculateSecondDerivativesLHS, KratosRansFastSuite)
{
    const auto& calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, const ProcessInfo& rProcessInfo) {
            rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
            const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
            noalias(rOutput) = rOutput * (1.0 - bossak_alpha);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansKEpsilonEpsilonRFC2D3N", "RansKEpsilonEpsilonAdjointRFC2D3N",
        TURBULENT_ENERGY_DISSIPATION_RATE_2, calculate_sensitivity_matrix, 1e-6, 1e-5);
}

} // namespace Testing
} // namespace Kratos