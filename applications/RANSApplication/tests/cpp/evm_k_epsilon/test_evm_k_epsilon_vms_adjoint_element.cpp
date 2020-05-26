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
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_EquationIdVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmKEpsilonVMSAdjoint2D3N");

    for (auto& element : r_adjoint_model_part.Elements())
    {
        std::vector<std::size_t> equation_ids{};
        element.EquationIdVector(equation_ids, r_adjoint_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(equation_ids.size(), element.GetGeometry().PointsNumber() * 3);
        std::size_t local_index = 0;
        for (std::size_t i = 0; i < equation_ids.size() / 3; ++i)
        {
            auto& r_node = element.GetGeometry()[i];
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

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_GetDofList, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmKEpsilonVMSAdjoint2D3N");

    for (auto& element : r_adjoint_model_part.Elements())
    {
        auto dofs = Element::DofsVectorType{};
        element.GetDofList(dofs, r_adjoint_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(dofs.size(), element.GetGeometry().PointsNumber() * 3);
        std::size_t local_index = 0;
        for (std::size_t i = 0; i < dofs.size() / 3; ++i)
        {
            auto& r_node = element.GetGeometry()[i];
            KRATOS_ERROR_IF(dofs[local_index++] != r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_X))
                << "Dofs mismatch";
            KRATOS_ERROR_IF(dofs[local_index++] != r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_Y))
                << "Dofs mismatch";
            KRATOS_ERROR_IF(dofs[local_index++] != r_node.pGetDof(ADJOINT_FLUID_SCALAR_1))
                << "Dofs mismatch";
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_GetValuesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmKEpsilonVMSAdjoint2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector element_values;
        r_element.GetValuesVector(element_values);

        Vector values(9);
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

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_GetFirstDerivativesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmKEpsilonVMSAdjoint2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);

        Vector element_values;
        r_element.GetFirstDerivativesVector(element_values);

        Vector values = ZeroVector(9);

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_GetSecondDerivativesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
        r_adjoint_model_part, "RansEvmKEpsilonVMSAdjoint2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector element_values;
        r_element.GetSecondDerivativesVector(element_values);

        Vector values(9);
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

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_CalculateFirstDerivativesLHS_VELOCITY,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) {
            rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmKEpsilonVMSAdjoint2D3N", VELOCITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_CalculateSensitivityMatrix, KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) {
            rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmKEpsilonVMSAdjoint2D3N", SHAPE_SENSITIVITY,
        calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_Calculate_RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) {
            rElement.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                               rOutput, rProcessInfo);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmKEpsilonVMSAdjoint2D3N", TURBULENT_KINETIC_ENERGY,
        calculate_sensitivity_matrix, 1e-5, 6e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_Calculate_RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) {
            rElement.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                               rOutput, rProcessInfo);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmKEpsilonVMSAdjoint2D3N", TURBULENT_ENERGY_DISSIPATION_RATE,
        calculate_sensitivity_matrix, 1e-5, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_CalculateFirstDerivativesLHS_PRESSURE,
                          KratosRansFastSuite)
{
    constexpr int domain_size = 2;
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) {
            rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmKEpsilonVMSAdjoint2D3N", PRESSURE,
        calculate_sensitivity_matrix, 1e-6, 5e-5, domain_size);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_CalculateSecondDerivativesLHS, KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) {
            rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
            const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
            noalias(rOutput) = rOutput * (1.0 - bossak_alpha);
        };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmKEpsilonVMSAdjoint2D3N", ACCELERATION,
        calculate_sensitivity_matrix, 1e-7, 1e-5);
}

} // namespace Testing
} // namespace Kratos