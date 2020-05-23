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
KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_GetValuesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmMonolithicKEpsilonVMSAdjoint2D");

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

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_GetFirstDerivativesVector,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmMonolithicKEpsilonVMSAdjoint2D");

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

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_GetSecondDerivativesVector,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmMonolithicKEpsilonVMSAdjoint2D");

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

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_GetDofList, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmMonolithicKEpsilonVMSAdjoint2D");

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

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_EquationIdVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmMonolithicKEpsilonVMSAdjoint2D");

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

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_VMS2D3N_VELOCITY,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", VELOCITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_VMS2D3N_PRESSURE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 2;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", PRESSURE,
        calculate_sensitivity_matrix, 1e-6, 1e-4, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_VMS2D3N_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", TURBULENT_KINETIC_ENERGY,
        calculate_sensitivity_matrix, 1e-7, 5e-4, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_VMS2D3N_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", TURBULENT_ENERGY_DISSIPATION_RATE,
        calculate_sensitivity_matrix, 1e-5, 1e-3, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_RansEvmKEpsilonK2D3N_VELOCITY,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        VELOCITY, calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, 3);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_RansEvmKEpsilonK2D3N_PRESSURE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 2;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        PRESSURE, calculate_sensitivity_matrix, 1e-6, 1e-4, derivative_offset, 3);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_RansEvmKEpsilonK2D3N_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", TURBULENT_KINETIC_ENERGY,
        calculate_sensitivity_matrix, 1e-7, 5e-4, derivative_offset, 3);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_RansEvmKEpsilonK2D3N_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        TURBULENT_ENERGY_DISSIPATION_RATE, calculate_sensitivity_matrix, 1e-5,
        1e-3, derivative_offset, 3);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_RansEvmKEpsilonEpsilon2D3N_VELOCITY,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonEpsilon2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        VELOCITY, calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, 4);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_RansEvmKEpsilonEpsilon2D3N_PRESSURE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 2;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonEpsilon2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        PRESSURE, calculate_sensitivity_matrix, 1e-6, 1e-4, derivative_offset, 4);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_RansEvmKEpsilonEpsilon2D3N_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonEpsilon2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        TURBULENT_KINETIC_ENERGY, calculate_sensitivity_matrix, 1e-7, 5e-4,
        derivative_offset, 4);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateFirstDerivativesLHS_RansEvmKEpsilonEpsilon2D3N_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonEpsilon2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        TURBULENT_ENERGY_DISSIPATION_RATE, calculate_sensitivity_matrix, 1e-5,
        1e-3, derivative_offset, 4);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSecondDerivativesLHS_VMS2D3N_ACCELERATION,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", ACCELERATION,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSecondDerivativesLHS_VMS2D3N_TURBULENT_KINETIC_ENERGY_RATE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", TURBULENT_KINETIC_ENERGY_RATE,
        calculate_sensitivity_matrix, 1e-7, 5e-4, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSecondDerivativesLHS_VMS2D3N_TURBULENT_ENERGY_DISSIPATION_RATE_2,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", TURBULENT_ENERGY_DISSIPATION_RATE_2,
        calculate_sensitivity_matrix, 1e-5, 1e-3, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSecondDerivativesLHS_RansEvmKEpsilonK2D3N_ACCELERATION,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", ACCELERATION,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, 3);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSecondDerivativesLHS_RansEvmKEpsilonK2D3N_TURBULENT_KINETIC_ENERGY_RATE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        TURBULENT_KINETIC_ENERGY_RATE, calculate_sensitivity_matrix, 1e-7, 5e-4,
        derivative_offset, 3);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSecondDerivativesLHS_RansEvmKEpsilonK2D3N_TURBULENT_ENERGY_DISSIPATION_RATE_2,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        TURBULENT_ENERGY_DISSIPATION_RATE_2, calculate_sensitivity_matrix, 1e-5,
        1e-3, derivative_offset, 3);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSecondDerivativesLHS_RansEvmKEpsilonEpsilon2D3N_ACCELERATION,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonEpsilon2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        ACCELERATION, calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, 4);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSecondDerivativesLHS_RansEvmKEpsilonEpsilon2D3N_TURBULENT_KINETIC_ENERGY_RATE,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 3;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonEpsilon2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        TURBULENT_KINETIC_ENERGY_RATE, calculate_sensitivity_matrix, 1e-7, 5e-4,
        derivative_offset, 4);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSecondDerivativesLHS_RansEvmKEpsilonEpsilon2D3N_TURBULENT_ENERGY_DISSIPATION_RATE_2,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 4;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double one_minus_bossak_alpha = 1.0 - rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * one_minus_bossak_alpha;
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonEpsilon2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D",
        TURBULENT_ENERGY_DISSIPATION_RATE_2, calculate_sensitivity_matrix, 1e-5,
        1e-3, derivative_offset, 4);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSensitivityMatrix_VMS2D3N,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", SHAPE_SENSITIVITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSensitivityMatrix_RansEvmKEpsilonK2D3N,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", SHAPE_SENSITIVITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, 3);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmMonolithicKEpsilonVMSAdjoint2D_CalculateSensitivityMatrix_RansEvmKEpsilonEpsilon2D3N,
                          KratosRansFastSuite)
{
    constexpr int derivative_offset = 0;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonEpsilon2D3N", "RansEvmMonolithicKEpsilonVMSAdjoint2D", SHAPE_SENSITIVITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5, derivative_offset, 4);
}

} // namespace Testing
} // namespace Kratos