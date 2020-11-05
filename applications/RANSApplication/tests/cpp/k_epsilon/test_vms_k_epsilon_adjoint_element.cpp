//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <functional>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "../stabilization_method_test_utilities.h"
#include "custom_elements/data_containers/k_epsilon/k_adjoint_element_data.h"
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_adjoint_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"
#include "includes/cfd_variables.h"
#include "test_utilities.h"
#include "custom_processes/rans_nut_k_epsilon_update_process.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "fluid_dynamics_application_variables.h"

#include "custom_elements/data_containers/k_epsilon/vms_rfc_k_epsilon_adjoint_element_data.h"
#include "custom_elements/two_equation_turbulence_model_adjoint_element.h"

#include "custom_utilities/adjoint_test_utilities.h"

namespace Kratos
{
namespace Testing
{

namespace
{
template<class TDataType>
void RunRansVMSKEpsilonAdjointElementTest(
    const std::string& rPrimalElementName,
    const Variable<TDataType>& rVariable,
    const std::function<void(Matrix&, ModelPart::ElementType&, const ProcessInfo&)>& rDerivativesRetrievalFunction,
    const IndexType EquationOffset,
    const IndexType DerivativesOffset,
    const double Delta,
    const double Tolerance)
{
    Model model;

    const auto& set_variable_values = [](ModelPart& rModelPart) {
        using namespace RansApplicationTestUtilities;

        RandomFillNodalHistoricalVariable(rModelPart, VELOCITY, -10.0, 10.0);
        RandomFillNodalHistoricalVariable(rModelPart, PRESSURE, -10.0, 10.0);
        RandomFillNodalHistoricalVariable(rModelPart, KINEMATIC_VISCOSITY, 1e-5, 1e-3);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 10.0, 100.0);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 1.0, 100.0);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 50.0);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE_2, 1.0, 50.0);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_2, 1.0, 10.0);

        auto& r_process_info = rModelPart.GetProcessInfo();
        r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 0.5);
        r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.5);
        r_process_info.SetValue(TURBULENCE_RANS_C1, 2.5);
        r_process_info.SetValue(TURBULENCE_RANS_C2, 3.5);
        r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);
        r_process_info.SetValue(DELTA_TIME, 1.5);
        r_process_info.SetValue(BOSSAK_ALPHA, 1.8);
        r_process_info.SetValue(DYNAMIC_TAU, 1.4);
        r_process_info.SetValue(RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, 2.0);
        r_process_info.SetValue(RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, 3.0);
    };

    // preparing primal model part
    const auto& add_primal_solution_step_variables = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(PRESSURE);
        rModelPart.AddNodalSolutionStepVariable(VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
    };
    const auto& add_primal_dofs = [](ModelPart::NodeType& rNode) {
        rNode.AddDof(TURBULENT_KINETIC_ENERGY);
        rNode.AddDof(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto& r_primal_model_part = RansApplicationTestUtilities::CreateTestModelPart(
        model, rPrimalElementName, "LineCondition2D2N",
        add_primal_solution_step_variables, add_primal_dofs);
    r_primal_model_part.GetProcessInfo()[DOMAIN_SIZE] = 2;
    set_variable_values(r_primal_model_part);

    // prepare adjoint model part
    const auto& add_adjoint_solution_step_variables = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(PRESSURE);
        rModelPart.AddNodalSolutionStepVariable(VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_1);
    };
    const auto& add_adjoint_dofs = [](ModelPart::NodeType& rNode) {
        rNode.AddDof(RANS_SCALAR_1_ADJOINT_1);
        rNode.AddDof(RANS_SCALAR_2_ADJOINT_1);
    };

    auto& r_adjoint_model_part = RansApplicationTestUtilities::CreateTestModelPart(
        model, "RansAdjointVMSKEpsilonRFC2D3N", "LineCondition2D2N",
        add_adjoint_solution_step_variables, add_adjoint_dofs);
    r_adjoint_model_part.GetProcessInfo()[DOMAIN_SIZE] = 2;
    set_variable_values(r_adjoint_model_part);

    const auto& update_function = [](ModelPart& rModelPart) {
        RansNutKEpsilonUpdateProcess nut_update(
            rModelPart.GetModel(), rModelPart.Name(), 2.1, 1e-12, 0);
        nut_update.ExecuteAfterCouplingSolveStep();
    };

    RansApplicationTestUtilities::RunAdjointElementTest<ModelPart::ElementsContainerType, TDataType>(
        r_primal_model_part, r_adjoint_model_part, update_function, rVariable,
        rDerivativesRetrievalFunction, EquationOffset, DerivativesOffset, Delta, Tolerance);
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateFirstDerivativesLHSKU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonKRFC2D3N", VELOCITY,
                                         derivatives_method, 3, 0, 2e-3, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateFirstDerivativesLHSKP, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonKRFC2D3N", PRESSURE,
                                         derivatives_method, 3, 2, 1e-5, 1e-7);
}

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateFirstDerivativesLHSKK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonKRFC2D3N", TURBULENT_KINETIC_ENERGY,
                                         derivatives_method, 3, 3, 2e-3, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateFirstDerivativesLHSKEpsilon, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonKRFC2D3N", TURBULENT_ENERGY_DISSIPATION_RATE,
                                         derivatives_method, 3, 4, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateFirstDerivativesLHSEpsilonU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonEpsilonRFC2D3N", VELOCITY,
                                         derivatives_method, 4, 0, 1e-4, 1e-2);
}

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateFirstDerivativesLHSEpsilonP, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonEpsilonRFC2D3N", PRESSURE,
                                         derivatives_method, 4, 2, 1e-5, 1e-7);
}

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateFirstDerivativesLHSEpsilonK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonEpsilonRFC2D3N", TURBULENT_KINETIC_ENERGY,
                                         derivatives_method, 4, 3, 1e-4, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateFirstDerivativesLHSEpsilonEpsilon, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonEpsilonRFC2D3N", TURBULENT_ENERGY_DISSIPATION_RATE,
                                         derivatives_method, 4, 4, 3e-4, 3e-3);
}

} // namespace Testing
} // namespace Kratos