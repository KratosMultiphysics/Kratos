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
#include "test_k_epsilon_utilities.h"

// Element data includes
#include "custom_elements/element_data/evm_k_epsilon/evm_k_epsilon_epsilon_adjoint_element_data.h"
#include "custom_elements/element_data/evm_k_epsilon/evm_k_epsilon_epsilon_element_data.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateEffectiveKinematicViscosityDerivative_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, double, primal_element_data_type, adjoint_element_data_type>(
        TURBULENT_KINETIC_ENERGY,
        RansModellingApplicationTestUtilities::ElementDataMethods::CalculateEffectiveKinematicViscosity,
        1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateEffectiveKinematicViscosityDerivative_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, double, primal_element_data_type, adjoint_element_data_type>(
        TURBULENT_ENERGY_DISSIPATION_RATE,
        RansModellingApplicationTestUtilities::ElementDataMethods::CalculateEffectiveKinematicViscosity,
        1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateEffectiveKinematicViscosityDerivative_VELOCITY,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, array_1d<double, 3>, primal_element_data_type, adjoint_element_data_type>(
        VELOCITY, RansModellingApplicationTestUtilities::ElementDataMethods::CalculateEffectiveKinematicViscosity,
        1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateEffectiveKinematicViscosityDerivative_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, array_1d<double, 3>, primal_element_data_type, adjoint_element_data_type>(
        SHAPE_SENSITIVITY,
        RansModellingApplicationTestUtilities::ElementDataMethods::CalculateEffectiveKinematicViscosity,
        1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateReactionTermDerivative_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, double, primal_element_data_type, adjoint_element_data_type>(
        TURBULENT_KINETIC_ENERGY,
        RansModellingApplicationTestUtilities::ElementDataMethods::CalculateReactionTerm,
        1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateReactionTermDerivative_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, double, primal_element_data_type, adjoint_element_data_type>(
        TURBULENT_ENERGY_DISSIPATION_RATE,
        RansModellingApplicationTestUtilities::ElementDataMethods::CalculateReactionTerm,
        1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateReactionTermDerivative_VELOCITY,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, array_1d<double, 3>, primal_element_data_type, adjoint_element_data_type>(
        VELOCITY, RansModellingApplicationTestUtilities::ElementDataMethods::CalculateReactionTerm,
        1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateReactionTermDerivative_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, array_1d<double, 3>, primal_element_data_type, adjoint_element_data_type>(
        SHAPE_SENSITIVITY,
        RansModellingApplicationTestUtilities::ElementDataMethods::CalculateReactionTerm,
        1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateSourceTermDerivative_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, double, primal_element_data_type, adjoint_element_data_type>(
        TURBULENT_KINETIC_ENERGY,
        RansModellingApplicationTestUtilities::ElementDataMethods::CalculateSourceTerm,
        1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateSourceTermDerivative_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, double, primal_element_data_type, adjoint_element_data_type>(
        TURBULENT_ENERGY_DISSIPATION_RATE,
        RansModellingApplicationTestUtilities::ElementDataMethods::CalculateSourceTerm,
        1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateSourceTermDerivative_VELOCITY,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, array_1d<double, 3>, primal_element_data_type, adjoint_element_data_type>(
        VELOCITY, RansModellingApplicationTestUtilities::ElementDataMethods::CalculateSourceTerm,
        1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonElementData_CalculateSourceTermDerivative_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    using primal_element_data_type =
        EvmKEpsilonElementDataUtilities::EpsilonElementData<2>;
    using adjoint_element_data_type =
        EvmKEpsilonAdjointElementDataUtilities::EpsilonAdjointElementData<2, 3>;

    RansEvmKEpsilonModel::RunRansEvmKEpsilonElementDataTest<2, 3, array_1d<double, 3>, primal_element_data_type, adjoint_element_data_type>(
        SHAPE_SENSITIVITY,
        RansModellingApplicationTestUtilities::ElementDataMethods::CalculateSourceTerm,
        1e-8, 1e-6);
}

} // namespace Testing
} // namespace Kratos