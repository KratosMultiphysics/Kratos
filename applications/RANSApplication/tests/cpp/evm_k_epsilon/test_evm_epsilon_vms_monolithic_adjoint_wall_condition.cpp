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
#include "fluid_dynamics_application_variables.h"
#include "rans_application_variables.h"
#include "test_k_epsilon_utilities.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_EquationIdVector,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
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

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_GetDofList, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
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

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_GetValuesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
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
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
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
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
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
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmVmsMonolithicAdjointWallCondition2D2N", VELOCITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_Calculate_RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                             rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmVmsMonolithicAdjointWallCondition2D2N",
        TURBULENT_KINETIC_ENERGY, calculate_sensitivity_matrix, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_Calculate_RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                             rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmVmsMonolithicAdjointWallCondition2D2N",
        TURBULENT_ENERGY_DISSIPATION_RATE, calculate_sensitivity_matrix, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_CalculateSecondDerivativesLHS,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmVmsMonolithicAdjointWallCondition2D2N", ACCELERATION,
        calculate_sensitivity_matrix, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmVmsMonolithicAdjointWallCondition2D2N_CalculateSensitivityMatrix,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmVmsMonolithicAdjointWallCondition2D2N", SHAPE_SENSITIVITY,
        calculate_sensitivity_matrix, 1e-9, 1e-5);
}
} // namespace Testing
} // namespace Kratos