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

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_EquationIdVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmKAdjoint2D3N");

    for (auto& element : r_adjoint_model_part.Elements())
    {
        std::vector<std::size_t> equation_ids{};
        element.EquationIdVector(equation_ids, r_adjoint_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(equation_ids.size(), element.GetGeometry().PointsNumber());

        for (std::size_t i = 0; i < equation_ids.size(); ++i)
        {
            KRATOS_ERROR_IF(
                equation_ids[i] !=
                element.GetGeometry()[i].GetDof(RANS_SCALAR_1_ADJOINT_1).EquationId())
                << "Equation id mismatch.";
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_GetDofList, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmKAdjoint2D3N");

    for (auto& element : r_adjoint_model_part.Elements())
    {
        auto dofs = Element::DofsVectorType{};
        element.GetDofList(dofs, r_adjoint_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(dofs.size(), element.GetGeometry().PointsNumber());
        for (std::size_t i = 0; i < dofs.size(); ++i)
        {
            KRATOS_ERROR_IF(dofs[i] != element.GetGeometry()[i].pGetDof(RANS_SCALAR_1_ADJOINT_1))
                << "Dofs mismatch.";
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_GetValuesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmKAdjoint2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector element_values;
        r_element.GetValuesVector(element_values);

        Vector values(3);
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const NodeType& r_node = r_geometry[i_node];
            values[local_index++] = r_node.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_1);
        }

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_GetFirstDerivativesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmKAdjoint2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);

        Vector element_values;
        r_element.GetFirstDerivativesVector(element_values);

        Vector values = ZeroVector(3);

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_GetSecondDerivativesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmKAdjoint2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector element_values;
        r_element.GetSecondDerivativesVector(element_values);

        Vector values(3);
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const NodeType& r_node = r_geometry[i_node];
            values[local_index++] = r_node.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_3);
        }

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateFirstDerivativesLHS, KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) {
            rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmKAdjoint2D3N", TURBULENT_KINETIC_ENERGY,
        calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_Calculate_RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmKAdjoint2D3N",
        TURBULENT_ENERGY_DISSIPATION_RATE, calculate_sensitivity_matrix, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateSecondDerivativesLHS, KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * (1.0 - bossak_alpha);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmKAdjoint2D3N",
        TURBULENT_KINETIC_ENERGY_RATE, calculate_sensitivity_matrix, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateSensitivityMatrix, KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmKAdjoint2D3N", SHAPE_SENSITIVITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_Calculate_RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.Calculate(RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "RansEvmKEpsilonK2D3N", "RansEvmKAdjoint2D3N", VELOCITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5);
}

} // namespace Testing
} // namespace Kratos