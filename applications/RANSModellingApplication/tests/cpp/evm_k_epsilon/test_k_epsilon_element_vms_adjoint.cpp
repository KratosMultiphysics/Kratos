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
#include "test_k_epsilon_utilities.h"
#include "custom_utilities/test_utilities.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjointElementVelocityDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKEpsilonVMSAdjoint2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_velocity[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjointElementShapeSensitivity,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKEpsilonVMSAdjoint2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
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

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjointElementTKEFirstDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKEpsilonVMSAdjoint2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-5, 6e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjointElementEpsilonFirstDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKEpsilonVMSAdjoint2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-5, 2e-4);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjointElementPressureFirstDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKEpsilonVMSAdjoint2D3N");

    const int domain_size = r_primal_model_part.GetProcessInfo()[DOMAIN_SIZE];

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(PRESSURE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-5, domain_size);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjointElementAccelerationDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKEpsilonVMSAdjoint2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "RansEvmKElementSensitivityMatrix"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(ACCELERATION);
        return r_velocity[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * (1.0 - bossak_alpha);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process, nut_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

} // namespace Testing
} // namespace Kratos