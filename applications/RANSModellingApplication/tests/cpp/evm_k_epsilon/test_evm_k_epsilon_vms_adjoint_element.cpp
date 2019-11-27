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
<<<<<<< HEAD
#include "custom_utilities/test_utilities.h"
#include "fluid_dynamics_application_variables.h"
#include "rans_modelling_application_variables.h"
#include "test_k_epsilon_utilities.h"
=======
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_adjoint_utilities.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_velocity_sensitivities_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_sensitivities_process.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "test_k_epsilon_utilities.h"
#include "fluid_dynamics_application_variables.h"
>>>>>>> origin/rans-modelling-application-k-epsilon-adjoints

namespace Kratos
{
namespace Testing
{
<<<<<<< HEAD
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_EquationIdVector, KratosRansFastSuite)
=======
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_EquationIdVector,
                          KratosRansFastSuite)
>>>>>>> origin/rans-modelling-application-k-epsilon-adjoints
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
<<<<<<< HEAD
            KRATOS_ERROR_IF(equation_ids[local_index++] !=
                            r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_X)->EquationId())
                << "EquationId mismatch";
            KRATOS_ERROR_IF(equation_ids[local_index++] !=
                            r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_Y)->EquationId())
                << "EquationId mismatch";
            KRATOS_ERROR_IF(equation_ids[local_index++] !=
                            r_node.pGetDof(ADJOINT_FLUID_SCALAR_1)->EquationId())
=======
            KRATOS_ERROR_IF(equation_ids[local_index++] != r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_X)->EquationId())
                << "EquationId mismatch";
            KRATOS_ERROR_IF(equation_ids[local_index++] != r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_Y)->EquationId())
                << "EquationId mismatch";
            KRATOS_ERROR_IF(equation_ids[local_index++] != r_node.pGetDof(ADJOINT_FLUID_SCALAR_1)->EquationId())
>>>>>>> origin/rans-modelling-application-k-epsilon-adjoints
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

<<<<<<< HEAD
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_GetValuesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
=======
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_GetValuesVector,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("test");
>>>>>>> origin/rans-modelling-application-k-epsilon-adjoints
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

<<<<<<< HEAD
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_GetFirstDerivativesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
=======
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_GetFirstDerivativesVector,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("test");
>>>>>>> origin/rans-modelling-application-k-epsilon-adjoints
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

<<<<<<< HEAD
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_GetSecondDerivativesVector, KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
=======
KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_GetSecondDerivativesVector,
                          KratosRansFastSuite)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("test");
>>>>>>> origin/rans-modelling-application-k-epsilon-adjoints
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
<<<<<<< HEAD
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
=======
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmKEpsilonVMSAdjoint2D3N",
        VELOCITY, calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_CalculateFirstDerivativesLHS_PRESSURE,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmKEpsilonVMSAdjoint2D3N",
        PRESSURE, calculate_sensitivity_matrix, 1e-5, 1e-5, 2);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_CalculateSensitivityMatrix,
                          KratosRansFastSuite)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmKEpsilonVMSAdjoint2D3N",
        SHAPE_SENSITIVITY, calculate_sensitivity_matrix, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_Calculate_RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                          KratosRansFastSuite1)
{
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> calculate_sensitivity_matrix =
        [](Matrix& rOutput, ElementType& rElement, ProcessInfo& rProcessInfo) -> void {
        rElement.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
        KRATOS_WATCH(rOutput);
    };

    RansEvmKEpsilonModel::RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
        "VMS2D3N", "RansEvmKEpsilonVMSAdjoint2D3N",
        TURBULENT_KINETIC_ENERGY, calculate_sensitivity_matrix, 1e-7, 1e-5);
}



// KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_Calculate_RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
//                           KratosRansFastSuite)
// {
//     Model primal_model;
//     ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
//     RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
//         r_primal_model_part, "VMS2D3N");

//     Model adjoint_model;
//     ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
//     RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
//         r_adjoint_model_part, "RansEvmKEpsilonVMSAdjoint2D3N");

//     Parameters empty_y_plus_parameters = Parameters(R"({
//         "model_part_name" : "test"
//     })");
//     RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
//         adjoint_model, empty_y_plus_parameters);
//     RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
//         adjoint_model, empty_y_plus_parameters);
//     RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
//         primal_model, empty_y_plus_parameters);

//     Parameters empty_nut_parameters = Parameters(R"({
//         "model_part_name" : "test"
//     })");
//     RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
//         adjoint_model, empty_nut_parameters);
//     RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
//         adjoint_model, empty_nut_parameters);
//     RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

//     auto perturb_variable = [](NodeType& rNode) -> double& {
//         return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
//     };

//     auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
//                                            ProcessInfo& rProcessInfo) {
//         rElement.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
//                            rOutput, rProcessInfo);
//     };

//     RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
//         r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
//         adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
//         nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
//         calculate_sensitivity_matrix, perturb_variable, 1e-5, 1e-3);
// }

// KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_CalculateFirstDerivativesLHS_PRESSURE,
//                           KratosRansFastSuite)
// {
//     Model primal_model;
//     ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
//     RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
//         r_primal_model_part, "VMS2D3N");

//     Model adjoint_model;
//     ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
//     RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
//         r_adjoint_model_part, "RansEvmKEpsilonVMSAdjoint2D3N");

//     const int domain_size = r_primal_model_part.GetProcessInfo()[DOMAIN_SIZE];

//     Parameters empty_y_plus_parameters = Parameters(R"({
//         "model_part_name" : "test"
//     })");
//     RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
//         adjoint_model, empty_y_plus_parameters);
//     RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
//         adjoint_model, empty_y_plus_parameters);
//     RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
//         primal_model, empty_y_plus_parameters);

//     Parameters empty_nut_parameters = Parameters(R"({
//         "model_part_name" : "test"
//     })");
//     RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
//         adjoint_model, empty_nut_parameters);
//     RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
//         adjoint_model, empty_nut_parameters);
//     RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

//     auto perturb_variable = [](NodeType& rNode) -> double& {
//         return rNode.FastGetSolutionStepValue(PRESSURE);
//     };

//     auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
//                                            ProcessInfo& rProcessInfo) {
//         rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
//     };

//     RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
//         r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
//         adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
//         nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
//         calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-3, domain_size);
// }

// KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjoint2D3N_CalculateSecondDerivativesLHS,
//                           KratosRansFastSuite)
// {
//     Model primal_model;
//     ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
//     RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
//         r_primal_model_part, "VMS2D3N");

//     Model adjoint_model;
//     ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
//     RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
//         r_adjoint_model_part, "RansEvmKEpsilonVMSAdjoint2D3N");

//     Parameters empty_y_plus_parameters = Parameters(R"({
//         "model_part_name" : "test"
//     })");
//     RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
//         adjoint_model, empty_y_plus_parameters);
//     RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
//         adjoint_model, empty_y_plus_parameters);
//     RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
//         primal_model, empty_y_plus_parameters);

//     Parameters empty_nut_parameters = Parameters(R"({
//         "model_part_name" : "test"
//     })");
//     RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
//         adjoint_model, empty_nut_parameters);
//     RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
//         adjoint_model, empty_nut_parameters);
//     RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

//     auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
//         array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(ACCELERATION);
//         return r_velocity[Dim];
//     };

//     auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
//                                            ProcessInfo& rProcessInfo) {
//         rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
//         const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
//         noalias(rOutput) = rOutput * (1.0 - bossak_alpha);
//     };

//     RansModellingApplicationTestUtilities::RunResidualVectorSensitivityTest(
//         r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
//         adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
//         nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
//         calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
// }
>>>>>>> origin/rans-modelling-application-k-epsilon-adjoints

} // namespace Testing
} // namespace Kratos