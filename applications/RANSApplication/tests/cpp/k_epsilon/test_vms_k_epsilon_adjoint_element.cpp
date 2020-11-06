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

        RandomFillNodalHistoricalVariable(rModelPart, DISTANCE, 1e-2, 1.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 0.1, 1.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 5.0, 10.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 10.0, 20.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE_2, 5.0, 10.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, VELOCITY, 50.0, 100.0, 0);

        RandomFillNodalHistoricalVariable(rModelPart, PRESSURE, 5.0, 10.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, EXTERNAL_PRESSURE, 50.0, 100.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, ACCELERATION, 2.0, 3.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, KINEMATIC_VISCOSITY, 3e-2, 3e-2, 0);
        RandomFillNodalHistoricalVariable(rModelPart, DENSITY, 200.0, 200.0, 0);

        RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_1, 55.0, 120.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_3, 35.0, 150.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_SCALAR_1, 45.0, 100.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, AUX_ADJOINT_FLUID_VECTOR_1, 45.0, 100.0, 0);

        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_1_ADJOINT_1, 15.0, 20.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_2_ADJOINT_1, 15.0, 20.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_1_ADJOINT_2, 25.0, 30.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_2_ADJOINT_2, 25.0, 30.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_1_ADJOINT_3, 25.0, 30.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_2_ADJOINT_3, 25.0, 30.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_AUX_ADJOINT_SCALAR_1, 25.0, 30.0, 0);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_AUX_ADJOINT_SCALAR_2, 25.0, 30.0, 0);

        RandomFillNodalHistoricalVariable(rModelPart, DISTANCE, 1e-2, 1.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 0.1, 1.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 5.0, 10.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 10.0, 20.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE_2, 5.0, 10.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, VELOCITY, 5.0, 10.0, 1);

        RandomFillNodalHistoricalVariable(rModelPart, PRESSURE, 5.0, 10.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, EXTERNAL_PRESSURE, 50.0, 100.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, ACCELERATION, 2.0, 3.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, KINEMATIC_VISCOSITY, 3e-2, 3e-2, 1);
        RandomFillNodalHistoricalVariable(rModelPart, DENSITY, 200.0, 200.0, 1);

        RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_1, 55.0, 120.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_3, 35.0, 150.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_SCALAR_1, 45.0, 100.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, AUX_ADJOINT_FLUID_VECTOR_1, 45.0, 100.0, 1);

        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_1_ADJOINT_1, 15.0, 20.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_2_ADJOINT_1, 15.0, 20.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_1_ADJOINT_2, 25.0, 30.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_2_ADJOINT_2, 25.0, 30.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_1_ADJOINT_3, 25.0, 30.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_SCALAR_2_ADJOINT_3, 25.0, 30.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_AUX_ADJOINT_SCALAR_1, 25.0, 30.0, 1);
        RandomFillNodalHistoricalVariable(rModelPart, RANS_AUX_ADJOINT_SCALAR_2, 25.0, 30.0, 1);

        auto& r_process_info = rModelPart.GetProcessInfo();
        r_process_info.SetValue(DOMAIN_SIZE, 2);
        r_process_info.SetValue(DYNAMIC_TAU, 0.04);
        r_process_info.SetValue(DELTA_TIME, 0.01);
        r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);
        r_process_info.SetValue(TURBULENCE_RANS_C1, 1.44);
        r_process_info.SetValue(TURBULENCE_RANS_C2, 1.92);
        r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 1.03);
        r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.3);
        r_process_info.SetValue(BOSSAK_ALPHA, -0.03);
        r_process_info.SetValue(WALL_SMOOTHNESS_BETA, 5.2);
        r_process_info.SetValue(WALL_VON_KARMAN, 0.41);
        r_process_info.SetValue(OSS_SWITCH, 0);
        r_process_info.SetValue(RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, 1.5);
        r_process_info.SetValue(RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, 3.0);
    };

    // preparing primal model part
    const auto& add_solution_step_variables = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(DISTANCE);
        rModelPart.AddNodalSolutionStepVariable(PRESSURE);
        rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(DENSITY);
        rModelPart.AddNodalSolutionStepVariable(VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
        rModelPart.AddNodalSolutionStepVariable(RELAXED_ACCELERATION);
        rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
        rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(RANS_Y_PLUS);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_3);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_3);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_1);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_2);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_3);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_SCALAR_1);
        rModelPart.AddNodalSolutionStepVariable(AUX_ADJOINT_FLUID_VECTOR_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUX_ADJOINT_SCALAR_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUX_ADJOINT_SCALAR_2);
    };
    const auto& add_dofs = [](ModelPart::NodeType& rNode) {
        rNode.AddDof(TURBULENT_KINETIC_ENERGY);
        rNode.AddDof(TURBULENT_ENERGY_DISSIPATION_RATE);
        rNode.AddDof(RANS_SCALAR_1_ADJOINT_1);
        rNode.AddDof(RANS_SCALAR_2_ADJOINT_1);

        rNode.AddDof(PRESSURE);
        rNode.AddDof(VELOCITY_X);
        rNode.AddDof(VELOCITY_Y);
        rNode.AddDof(VELOCITY_Z);

        rNode.AddDof(ADJOINT_FLUID_SCALAR_1);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_X);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_Y);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_Z);
    };

    auto& r_primal_model_part = RansApplicationTestUtilities::CreateTestModelPart(
        model, rPrimalElementName, "LineCondition2D2N",
        add_solution_step_variables, add_dofs);
    set_variable_values(r_primal_model_part);

    // prepare adjoint model part
    auto& r_adjoint_model_part = RansApplicationTestUtilities::CreateTestModelPart(
        model, "RansAdjointVMSKEpsilonRFC2D3N", "LineCondition2D2N",
        add_solution_step_variables, add_dofs);
    set_variable_values(r_adjoint_model_part);

    const auto& update_function = [](ModelPart& rModelPart) {
        RansNutKEpsilonUpdateProcess nut_update(
            rModelPart.GetModel(), rModelPart.Name(), 2.1, 1e-12, 0);
        nut_update.ExecuteAfterCouplingSolveStep();

        const int number_of_nodes = rModelPart.NumberOfNodes();
        const double bossak_alpha = rModelPart.GetProcessInfo()[BOSSAK_ALPHA];

        for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
            auto& r_node = *(rModelPart.NodesBegin() + i_node);
            r_node.FastGetSolutionStepValue(VISCOSITY) = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) + r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) = RansApplicationTestUtilities::ComputeRelaxedVariableRate(bossak_alpha, TURBULENT_KINETIC_ENERGY_RATE, r_node);
            r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) = RansApplicationTestUtilities::ComputeRelaxedVariableRate(bossak_alpha, TURBULENT_ENERGY_DISSIPATION_RATE_2, r_node);
            r_node.FastGetSolutionStepValue(RELAXED_ACCELERATION) = RansApplicationTestUtilities::ComputeRelaxedVariableRate(bossak_alpha, ACCELERATION, r_node);
        }
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
                                         derivatives_method, 3, 0, 1e-5, 1e-5);
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
                                         derivatives_method, 3, 3, 1e-7, 1e-5);
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

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateSecondDerivativesLHSK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rMatrix) = rMatrix * (1.0 - bossak_alpha);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonKRFC2D3N", TURBULENT_KINETIC_ENERGY_RATE,
                                         derivatives_method, 3, 3, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateFirstDerivativesLHSEpsilonU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonEpsilonRFC2D3N", VELOCITY,
                                         derivatives_method, 4, 0, 1e-5, 1e-5);
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
                                         derivatives_method, 4, 3, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateFirstDerivativesLHSEpsilonEpsilon, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonEpsilonRFC2D3N", TURBULENT_ENERGY_DISSIPATION_RATE,
                                         derivatives_method, 4, 4, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RANSVMSKEpsilonAdjointCalculateSecondDerivativesLHSEpsilon, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rMatrix) = rMatrix * (1.0 - bossak_alpha);
    };
    RunRansVMSKEpsilonAdjointElementTest("RansKEpsilonEpsilonRFC2D3N", TURBULENT_ENERGY_DISSIPATION_RATE_2,
                                         derivatives_method, 4, 4, 1e-6, 1e-5);
}

} // namespace Testing
} // namespace Kratos