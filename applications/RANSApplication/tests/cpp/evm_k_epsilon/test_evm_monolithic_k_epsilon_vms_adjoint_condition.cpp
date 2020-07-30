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
KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_GetValuesVector,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
        r_adjoint_model_part,
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfConditions(); ++i_element)
    {
        Condition& r_element = *(r_adjoint_model_part.ConditionsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector element_values;
        r_element.GetValuesVector(element_values);

        Vector values(10);
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

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_GetFirstDerivativesVector,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
        r_adjoint_model_part,
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfConditions(); ++i_element)
    {
        Condition& r_element = *(r_adjoint_model_part.ConditionsBegin() + i_element);

        Vector element_values;
        r_element.GetFirstDerivativesVector(element_values);

        Vector values = ZeroVector(10);

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_GetSecondDerivativesVector,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
        r_adjoint_model_part,
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfConditions(); ++i_element)
    {
        Condition& r_element = *(r_adjoint_model_part.ConditionsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector element_values;
        r_element.GetSecondDerivativesVector(element_values);

        Vector values(10);
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

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_GetDofList,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
        r_adjoint_model_part,
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N");

    ProcessInfo& r_process_info = r_adjoint_model_part.GetProcessInfo();

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfConditions(); ++i_element)
    {
        Condition& r_element = *(r_adjoint_model_part.ConditionsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Condition::DofsVectorType element_dofs;
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

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_EquationIdVector,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
        r_adjoint_model_part,
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N");

    ProcessInfo& r_process_info = r_adjoint_model_part.GetProcessInfo();

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfConditions(); ++i_element)
    {
        Condition& r_element = *(r_adjoint_model_part.ConditionsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Condition::EquationIdVectorType element_eq_ids;
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

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonVmsMonolithicWall2D2N_VELOCITY,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", VELOCITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonVmsMonolithicWall2D2N_PRESSURE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 2;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", PRESSURE,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonVmsMonolithicWall2D2N_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", TURBULENT_KINETIC_ENERGY,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonVmsMonolithicWall2D2N_TURBULENT_ENERGY_DISSIPATON_RATE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", TURBULENT_ENERGY_DISSIPATION_RATE,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonEpsilonWall2D2N_VELOCITY,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 4;
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonEpsilonWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", VELOCITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonEpsilonWall2D2N_PRESSURE,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 4;
    constexpr int derivative_offset = 2;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonEpsilonWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", PRESSURE,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonEpsilonWall2D2N_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 4;
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonEpsilonWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", TURBULENT_KINETIC_ENERGY,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonEpsilonWall2D2N_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 4;
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonEpsilonWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        TURBULENT_ENERGY_DISSIPATION_RATE, calculate_sensitivity_matrix, 1e-7,
        1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonKCondition2D2N_VELOCITY,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 3;
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "LineCondition2D2N", "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", VELOCITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonKCondition2D2N_PRESSURE,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 3;
    constexpr int derivative_offset = 2;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "LineCondition2D2N", "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", PRESSURE,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonKCondition2D2N_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 3;
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "LineCondition2D2N", "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        TURBULENT_KINETIC_ENERGY, calculate_sensitivity_matrix, 1e-7, 1e-5,
        derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateFirstDerivativesLHS_RansEvmKEpsilonKCondition2D2N_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 3;
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "LineCondition2D2N", "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        TURBULENT_ENERGY_DISSIPATION_RATE, calculate_sensitivity_matrix, 1e-7,
        1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSecondDerivativesLHS_RansEvmKEpsilonVmsMonolithicWall2D2N_ACCELERATION,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 0;
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", ACCELERATION,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSecondDerivativesLHS_RansEvmKEpsilonVmsMonolithicWall2D2N_TURBULENT_KINETIC_ENERGY_RATE,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 0;
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        TURBULENT_KINETIC_ENERGY_RATE, calculate_sensitivity_matrix, 1e-7, 1e-5,
        derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSecondDerivativesLHS_RansEvmKEpsilonVmsMonolithicWall2D2N_TURBULENT_ENERGY_DISSIPATION_RATE_2,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 0;
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        TURBULENT_ENERGY_DISSIPATION_RATE_2, calculate_sensitivity_matrix, 1e-7,
        1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSecondDerivativesLHS_RansEvmKEpsilonKCondition2D2N_ACCELERATION,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 3;
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "LineCondition2D2N", "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        ACCELERATION, calculate_sensitivity_matrix, 1e-7, 1e-5,
        derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSecondDerivativesLHS_RansEvmKEpsilonKCondition2D2N_TURBULENT_KINETIC_ENERGY_RATE,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 3;
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "LineCondition2D2N", "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        TURBULENT_KINETIC_ENERGY_RATE, calculate_sensitivity_matrix, 1e-7, 1e-5,
        derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSecondDerivativesLHS_RansEvmKEpsilonKCondition2D2N_TURBULENT_ENERGY_DISSIPATION_RATE_2,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 3;
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "LineCondition2D2N", "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        TURBULENT_ENERGY_DISSIPATION_RATE_2, calculate_sensitivity_matrix, 1e-7,
        1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSecondDerivativesLHS_RansEvmKEpsilonEpsilonWall2D2N_ACCELERATION,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 4;
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonEpsilonWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", ACCELERATION,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSecondDerivativesLHS_RansEvmKEpsilonEpsilonWall2D2N_TURBULENT_KINETIC_ENERGY_RATE,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 4;
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonEpsilonWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        TURBULENT_KINETIC_ENERGY_RATE, calculate_sensitivity_matrix, 1e-7, 1e-5,
        derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSecondDerivativesLHS_RansEvmKEpsilonEpsilonWall2D2N_TURBULENT_ENERGY_DISSIPATION_RATE_2,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 4;
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonEpsilonWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        TURBULENT_ENERGY_DISSIPATION_RATE_2, calculate_sensitivity_matrix, 1e-7,
        1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSensitivityMatrix_RansEvmKEpsilonVmsMonolithicWall2D2N,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 0;
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonVmsMonolithicWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", SHAPE_SENSITIVITY,
        calculate_sensitivity_matrix, 1e-8, 1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSensitivityMatrix_RansEvmKEpsilonEpsilonWall2D2N,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 4;
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "RansEvmKEpsilonEpsilonWall2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", SHAPE_SENSITIVITY,
        calculate_sensitivity_matrix, 1e-8, 1e-5, derivative_offset, equation_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N_CalculateSensitivityMatrix_RansEvmKEpsilonKCondition2D2N,
                          KratosRansFastSuite)
{
    constexpr int equation_offset = 4;
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ConditionType& rCondition, ProcessInfo& rProcessInfo) -> void {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
        "LineCondition2D2N",
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N", SHAPE_SENSITIVITY,
        calculate_sensitivity_matrix, 1e-8, 1e-5, derivative_offset, equation_offset);
}

} // namespace Testing
} // namespace Kratos
