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
#include "test_k_epsilon_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_EquationIdVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

    for (auto& element : r_adjoint_model_part.Elements())
    {
        std::vector<std::size_t> equation_ids{};
        element.EquationIdVector(equation_ids, r_adjoint_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(equation_ids.size(), element.GetGeometry().PointsNumber() * 3);
        std::size_t local_index = 0;
        for (std::size_t i = 0; i < equation_ids.size() / 3; ++i)
        {
            auto& r_node = element.GetGeometry()[i];
            KRATOS_ERROR_IF(equation_ids[local_index++] != r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_X)->EquationId())
                << "EquationId mismatch";
            KRATOS_ERROR_IF(equation_ids[local_index++] != r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_Y)->EquationId())
                << "EquationId mismatch";
            KRATOS_ERROR_IF(equation_ids[local_index++] != r_node.pGetDof(ADJOINT_FLUID_SCALAR_1)->EquationId())
                << "EquationId mismatch";
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_GetDofList, RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

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

KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_GetValuesVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

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

KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_GetFirstDerivativesVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

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

KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_GetSecondDerivativesVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

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

KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_CalculateFirstDerivativesLHS_VELOCITY,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

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
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_velocity[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_CalculateSensitivityMatrix,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

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

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_Calculate_RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

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
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-5, 6e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_Calculate_RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

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
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-5, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_CalculateFirstDerivativesLHS_PRESSURE,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

    const int domain_size = r_primal_model_part.GetProcessInfo()[DOMAIN_SIZE];

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

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(PRESSURE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-3, domain_size);
}

KRATOS_TEST_CASE_IN_SUITE(RansRansEvmKEpsilonVMSAdjoint2D3N_CalculateSecondDerivativesLHS,
                          RANSModellingApplicationInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "EVMKEpsilonVMSAdjointElement2D3N");

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
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

} // namespace Testing
} // namespace Kratos