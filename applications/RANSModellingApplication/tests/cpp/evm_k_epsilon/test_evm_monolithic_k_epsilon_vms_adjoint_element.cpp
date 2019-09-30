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
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_adjoint_utilities.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_velocity_sensitivities_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_sensitivities_process.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "fluid_dynamics_application_variables.h"
#include "test_k_epsilon_utilities.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_CalculateFirstDerivativesLHS_VMS2D3N,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable_velocity = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_velocity[Dim];
    };

    auto perturb_variable_pressure = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(PRESSURE);
    };

    auto perturb_variable_tke = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto perturb_variable_epsilon = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    const double delta = 1e-5;
    const double tolerance = 1e-3;

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_velocity, delta, tolerance, 0, 0);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_pressure, delta, tolerance, 2, 0);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable_tke, delta, tolerance, 3, 0);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_epsilon, delta, tolerance, 4, 0);
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_CalculateFirstDerivativesLHS_EVMKElement2D3N,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RansEvmK2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable_velocity = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_velocity[Dim];
    };

    auto perturb_variable_pressure = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(PRESSURE);
    };

    auto perturb_variable_tke = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto perturb_variable_epsilon = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    const double delta = 1e-5;
    const double tolerance = 1e-3;

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_velocity, delta, tolerance, 0, 3);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_pressure, delta, tolerance, 2, 3);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable_tke, delta, tolerance, 3, 3);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_epsilon, delta, tolerance, 4, 3);
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_CalculateFirstDerivativesLHS_EVMEpsilonElement2D3N,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RansEvmEpsilon2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable_velocity = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_velocity[Dim];
    };

    auto perturb_variable_pressure = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(PRESSURE);
    };

    auto perturb_variable_tke = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto perturb_variable_epsilon = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    const double delta = 1e-5;
    const double tolerance = 1e-3;

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_velocity, delta, tolerance, 0, 4);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_pressure, delta, tolerance, 2, 4);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable_tke, delta, tolerance, 3, 4);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_epsilon, delta, tolerance, 4, 4);
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_CalculateSecondDerivativesLHS_VMS2D3N,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    const double one_minus_bossak_alpha =
        1.0 - r_primal_model_part.GetProcessInfo()[BOSSAK_ALPHA];

    auto perturb_variable_velocity = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(ACCELERATION);
        return r_velocity[Dim];
    };

    auto perturb_variable_tke = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE);
    };

    auto perturb_variable_epsilon = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2);
    };

    auto calculate_sensitivity_matrix = [one_minus_bossak_alpha](
                                            Matrix& rOutput, Element& rElement,
                                            ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    const double delta = 1e-5;
    const double tolerance = 1e-3;

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_velocity, delta, tolerance, 0, 0);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable_tke, delta, tolerance, 3, 0);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_epsilon, delta, tolerance, 4, 0);
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_CalculateSecondDerivativesLHS_EVMKElement2D3N,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RansEvmK2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    const double one_minus_bossak_alpha =
        1.0 - r_primal_model_part.GetProcessInfo()[BOSSAK_ALPHA];

    auto perturb_variable_velocity = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(ACCELERATION);
        return r_velocity[Dim];
    };

    auto perturb_variable_tke = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE);
    };

    auto perturb_variable_epsilon = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2);
    };

    auto calculate_sensitivity_matrix = [one_minus_bossak_alpha](
                                            Matrix& rOutput, Element& rElement,
                                            ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    const double delta = 1e-5;
    const double tolerance = 1e-3;

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_velocity, delta, tolerance, 0, 3);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable_tke, delta, tolerance, 3, 3);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_epsilon, delta, tolerance, 4, 3);
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_CalculateSecondDerivativesLHS_EVMEpsilonElement2D3N,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RansEvmEpsilon2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    const double one_minus_bossak_alpha =
        1.0 - r_primal_model_part.GetProcessInfo()[BOSSAK_ALPHA];

    auto perturb_variable_velocity = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(ACCELERATION);
        return r_velocity[Dim];
    };

    auto perturb_variable_tke = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE);
    };

    auto perturb_variable_epsilon = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2);
    };

    auto calculate_sensitivity_matrix = [one_minus_bossak_alpha](
                                            Matrix& rOutput, Element& rElement,
                                            ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    const double delta = 1e-5;
    const double tolerance = 1e-3;

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_velocity, delta, tolerance, 0, 4);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable_tke, delta, tolerance, 3, 4);

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        primal_nut_process, adjoint_y_plus_process, adjoint_nut_process,
        y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, calculate_sensitivity_matrix,
        perturb_variable_epsilon, delta, tolerance, 4, 4);
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_CalculateSensitivityMatrix_VMS2D3N,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_coordinates = rNode.Coordinates();
        return r_coordinates[Dim];
    };

    auto calculate_sensitivity_matrix_velocity_pressure =
        [](Matrix& rOutput, Element& rElement, ProcessInfo& rProcessInfo) {
            rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
        };

    const double delta = 1e-5;
    const double tolerance = 1e-3;

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix_velocity_pressure, perturb_variable, delta,
        tolerance, 0, 0);
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_CalculateSensitivityMatrix_EVMKElement2D3N,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RansEvmK2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_coordinates = rNode.Coordinates();
        return r_coordinates[Dim];
    };

    auto calculate_sensitivity_matrix_tke = [](Matrix& rOutput, Element& rElement,
                                               ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    const double delta = 1e-5;
    const double tolerance = 1e-3;

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix_tke, perturb_variable, delta, tolerance, 0, 3);
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_CalculateSensitivityMatrix_EVMEpsilonElement2D3N,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RansEvmEpsilon2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_coordinates = rNode.Coordinates();
        return r_coordinates[Dim];
    };

    auto calculate_sensitivity_matrix_epsilon = [](Matrix& rOutput, Element& rElement,
                                                   ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    const double delta = 1e-5;
    const double tolerance = 1e-3;

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix_epsilon, perturb_variable, delta, tolerance, 0, 4);
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_GetValuesVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector element_values;
        r_element.GetValuesVector(element_values);

        Vector values(15);
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const NodeType& r_node = r_geometry[i_node];
            const array_1d<double, 3>& r_vector =
                r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1);
            values[local_index++] = r_vector[0];
            values[local_index++] = r_vector[1];
            values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1);
            values[local_index++] = r_node.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_1);
            values[local_index++] = r_node.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_1);
        }

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_GetFirstDerivativesVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);

        Vector element_values;
        r_element.GetFirstDerivativesVector(element_values);

        Vector values = ZeroVector(15);

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_GetSecondDerivativesVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector element_values;
        r_element.GetSecondDerivativesVector(element_values);

        Vector values(15);
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const NodeType& r_node = r_geometry[i_node];
            const array_1d<double, 3>& r_vector =
                r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);
            values[local_index++] = r_vector[0];
            values[local_index++] = r_vector[1];
            values[local_index++] = 0.0;
            values[local_index++] = r_node.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_3);
            values[local_index++] = r_node.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_3);
        }

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_GetDofList,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    ProcessInfo& r_process_info = r_adjoint_model_part.GetProcessInfo();

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Element::DofsVectorType element_dofs;
        r_element.GetDofList(element_dofs, r_process_info);

        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = r_geometry[i_node];
            KRATOS_ERROR_IF(element_dofs[local_index++] != r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_X))
                << "Dofs mismatch";
            KRATOS_ERROR_IF(element_dofs[local_index++] != r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_Y))
                << "Dofs mismatch";
            KRATOS_ERROR_IF(element_dofs[local_index++] != r_node.pGetDof(ADJOINT_FLUID_SCALAR_1))
                << "Dofs mismatch";
            KRATOS_ERROR_IF(element_dofs[local_index++] != r_node.pGetDof(RANS_SCALAR_1_ADJOINT_1))
                << "Dofs mismatch";
            KRATOS_ERROR_IF(element_dofs[local_index++] != r_node.pGetDof(RANS_SCALAR_2_ADJOINT_1))
                << "Dofs mismatch";
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(EVMMonolithicKEpsilonVMSAdjointElement2D3N_EquationIdVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMMonolithicKEpsilonVMSAdjointElement2D3N");

    ProcessInfo& r_process_info = r_adjoint_model_part.GetProcessInfo();

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Element::EquationIdVectorType element_eq_ids;
        r_element.EquationIdVector(element_eq_ids, r_process_info);

        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = r_geometry[i_node];
            KRATOS_ERROR_IF(element_eq_ids[local_index++] !=
                            r_node.GetDof(ADJOINT_FLUID_VECTOR_1_X).EquationId())
                << "Equation id mismatch";
            KRATOS_ERROR_IF(element_eq_ids[local_index++] !=
                            r_node.GetDof(ADJOINT_FLUID_VECTOR_1_Y).EquationId())
                << "Equation id mismatch";
            KRATOS_ERROR_IF(element_eq_ids[local_index++] !=
                            r_node.GetDof(ADJOINT_FLUID_SCALAR_1).EquationId())
                << "Equation id mismatch";
            KRATOS_ERROR_IF(element_eq_ids[local_index++] !=
                            r_node.GetDof(RANS_SCALAR_1_ADJOINT_1).EquationId())
                << "Equation id mismatch";
            KRATOS_ERROR_IF(element_eq_ids[local_index++] !=
                            r_node.GetDof(RANS_SCALAR_2_ADJOINT_1).EquationId())
                << "Equation id mismatch";
        }
    }
}

} // namespace Testing
} // namespace Kratos