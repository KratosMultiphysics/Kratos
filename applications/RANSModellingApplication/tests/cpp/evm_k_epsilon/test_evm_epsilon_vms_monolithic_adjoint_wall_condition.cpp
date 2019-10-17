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
KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_EquationIdVector,
                          RANSEvModelsKEpsilonConditionMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_adjoint_model_part, "RansEvmVmsMonolithicAdjointWallCondition2D2N");

    for (auto& condition : r_adjoint_model_part.Conditions())
    {
        std::vector<std::size_t> equation_ids{};
        condition.EquationIdVector(equation_ids, r_adjoint_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(equation_ids.size(), condition.GetGeometry().PointsNumber() * 3);
        std::size_t local_index = 0;
        for (std::size_t i = 0; i < equation_ids.size() / 3; ++i)
        {
            auto& r_node = condition.GetGeometry()[i];
            KRATOS_ERROR_IF(equation_ids[local_index++] !=
                            r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_X)->EquationId())
                << "EquationId mismatch";
            KRATOS_ERROR_IF(equation_ids[local_index++] !=
                            r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_Y)->EquationId())
                << "EquationId mismatch";
            KRATOS_ERROR_IF(equation_ids[local_index++] !=
                            r_node.pGetDof(ADJOINT_FLUID_SCALAR_1)->EquationId())
                << "EquationId mismatch";
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_GetDofList,
                          RANSEvModelsKEpsilonConditionMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_adjoint_model_part, "RansEvmVmsMonolithicAdjointWallCondition2D2N");

    for (auto& condition : r_adjoint_model_part.Conditions())
    {
        auto dofs = Condition::DofsVectorType{};
        condition.GetDofList(dofs, r_adjoint_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(dofs.size(), condition.GetGeometry().PointsNumber() * 3);
        std::size_t local_index = 0;
        for (std::size_t i = 0; i < dofs.size() / 3; ++i)
        {
            auto& r_node = condition.GetGeometry()[i];
            KRATOS_ERROR_IF(dofs[local_index++] != r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_X))
                << "Dofs mismatch";
            KRATOS_ERROR_IF(dofs[local_index++] != r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_Y))
                << "Dofs mismatch";
            KRATOS_ERROR_IF(dofs[local_index++] != r_node.pGetDof(ADJOINT_FLUID_SCALAR_1))
                << "Dofs mismatch";
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_GetValuesVector,
                          RANSEvModelsKEpsilonConditionMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_adjoint_model_part, "RansEvmVmsMonolithicAdjointWallCondition2D2N");

    for (IndexType i_condition = 0;
         i_condition < r_adjoint_model_part.NumberOfConditions(); ++i_condition)
    {
        Condition& r_condition = *(r_adjoint_model_part.ConditionsBegin() + i_condition);
        GeometryType& r_geometry = r_condition.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector condition_values;
        r_condition.GetValuesVector(condition_values);

        Vector values(6);
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const NodeType& r_node = r_geometry[i_node];
            const array_1d<double, 3>& r_vector =
                r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1);
            values[local_index++] = r_vector[0];
            values[local_index++] = r_vector[1];
            values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1);
        }

        RansModellingApplicationTestUtilities::CheckNear(condition_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_GetFirstDerivativesVector,
                          RANSEvModelsKEpsilonConditionMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_adjoint_model_part, "RansEvmVmsMonolithicAdjointWallCondition2D2N");

    for (IndexType i_condition = 0;
         i_condition < r_adjoint_model_part.NumberOfConditions(); ++i_condition)
    {
        Condition& r_condition = *(r_adjoint_model_part.ConditionsBegin() + i_condition);

        Vector condition_values;
        r_condition.GetFirstDerivativesVector(condition_values);

        Vector values = ZeroVector(6);

        RansModellingApplicationTestUtilities::CheckNear(condition_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_GetSecondDerivativesVector,
                          RANSEvModelsKEpsilonConditionMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_adjoint_model_part, "RansEvmVmsMonolithicAdjointWallCondition2D2N");

    for (IndexType i_condition = 0;
         i_condition < r_adjoint_model_part.NumberOfConditions(); ++i_condition)
    {
        Condition& r_condition = *(r_adjoint_model_part.ConditionsBegin() + i_condition);
        GeometryType& r_geometry = r_condition.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector condition_values;
        r_condition.GetSecondDerivativesVector(condition_values);

        Vector values(6);
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const NodeType& r_node = r_geometry[i_node];
            const array_1d<double, 3>& r_vector =
                r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);
            values[local_index++] = r_vector[0];
            values[local_index++] = r_vector[1];
            values[local_index++] = 0.0;
        }

        RansModellingApplicationTestUtilities::CheckNear(condition_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_CalculateFirstDerivativesLHS,
                          RANSModellingApplicationConditionInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_primal_model_part, "RansEvmVmsMonolithicWallCondition2D2N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_adjoint_model_part, "RansEvmVmsMonolithicAdjointWallCondition2D2N");

    Process dummy_y_plus_process;

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
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

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Condition& rCondition,
                                           ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, dummy_y_plus_process,
        primal_nut_process, dummy_y_plus_process, adjoint_nut_process, dummy_y_plus_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_Calculate_RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                          RANSModellingApplicationConditionInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_primal_model_part, "RansEvmVmsMonolithicWallCondition2D2N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_adjoint_model_part, "RansEvmVmsMonolithicAdjointWallCondition2D2N");

    Process dummy_y_plus_process;

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Condition& rCondition,
                                           ProcessInfo& rProcessInfo) {
        rCondition.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                             rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, dummy_y_plus_process,
        primal_nut_process, dummy_y_plus_process, adjoint_nut_process, dummy_y_plus_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_Calculate_RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                          RANSModellingApplicationConditionInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_primal_model_part, "RansEvmVmsMonolithicWallCondition2D2N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_adjoint_model_part, "RansEvmVmsMonolithicAdjointWallCondition2D2N");

    Process dummy_y_plus_process;

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Condition& rCondition,
                                           ProcessInfo& rProcessInfo) {
        rCondition.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                             rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, dummy_y_plus_process,
        primal_nut_process, dummy_y_plus_process, adjoint_nut_process, dummy_y_plus_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_CalculateSecondDerivativesLHS,
                          RANSModellingApplicationConditionInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_primal_model_part, "RansEvmVmsMonolithicWallCondition2D2N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_adjoint_model_part, "RansEvmVmsMonolithicAdjointWallCondition2D2N");

    Process dummy_y_plus_process;

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION);
        return r_acceleration[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Condition& rCondition,
                                           ProcessInfo& rProcessInfo) {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, dummy_y_plus_process,
        primal_nut_process, dummy_y_plus_process, adjoint_nut_process, dummy_y_plus_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_CalculateSensitivityMatrix,
                          RANSModellingApplicationConditionInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_primal_model_part, "RansEvmVmsMonolithicWallCondition2D2N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonConditionTestModelPart(
        r_adjoint_model_part, "RansEvmVmsMonolithicAdjointWallCondition2D2N");

    Process dummy_y_plus_process;

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

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Condition& rCondition,
                                           ProcessInfo& rProcessInfo) {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, dummy_y_plus_process, primal_nut_process,
        dummy_y_plus_process, adjoint_nut_process, dummy_y_plus_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-9, 1e-5);
}

} // namespace Testing
} // namespace Kratos