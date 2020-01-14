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

#include "test_k_epsilon_utilities.h"

// System includes
#include <random>

// External includes
#include "fluid_dynamics_application_variables.h"

// Project includes
#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"

#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_sensitivities_process.h"

namespace Kratos
{
namespace Testing
{
typedef ModelPart::NodeType NodeType;

typedef ModelPart::ElementType ElementType;

typedef Geometry<NodeType> GeometryType;

typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

/**
 * Auxiliar function to generate a triangular element to be tested.
 */
namespace RansEvmKEpsilonModel
{
void AddVariablesToModelPart(ModelPart& rModelPart)
{
    // Set buffer size
    rModelPart.SetBufferSize(2);

    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(RELAXED_ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
    rModelPart.AddNodalSolutionStepVariable(RANS_Y_PLUS);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_2);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_2);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_3);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_3);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_1);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_2);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_3);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_SCALAR_1);
    rModelPart.AddNodalSolutionStepVariable(AUX_ADJOINT_FLUID_VECTOR_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUX_ADJOINT_SCALAR_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUX_ADJOINT_SCALAR_2);
}

void InitializeProcessInfo(ModelPart& rModelPart)
{
    // Process info creation
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rModelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.04);
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, 0.01);
    rModelPart.GetProcessInfo().SetValue(TURBULENCE_RANS_C_MU, 0.09);
    rModelPart.GetProcessInfo().SetValue(TURBULENCE_RANS_C1, 1.44);
    rModelPart.GetProcessInfo().SetValue(TURBULENCE_RANS_C2, 1.92);
    rModelPart.GetProcessInfo().SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 1.03);
    rModelPart.GetProcessInfo().SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.3);
    rModelPart.GetProcessInfo().SetValue(BOSSAK_ALPHA, -0.03);
    rModelPart.GetProcessInfo().SetValue(WALL_SMOOTHNESS_BETA, 5.2);
    rModelPart.GetProcessInfo().SetValue(WALL_VON_KARMAN, 0.41);
    rModelPart.GetProcessInfo().SetValue(OSS_SWITCH, 0);
    rModelPart.GetProcessInfo().SetValue(IS_CO_SOLVING_PROCESS_ACTIVE, true);

    // Set the element properties
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    p_elem_prop->SetValue(KINEMATIC_VISCOSITY, 3.0e-02);
}

void CreateModelPartNodes(ModelPart& rModelPart)
{
    // Element creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);

    VariableUtils().AddDof<Variable<double>>(TURBULENT_KINETIC_ENERGY, rModelPart);
    VariableUtils().AddDof<Variable<double>>(TURBULENT_ENERGY_DISSIPATION_RATE, rModelPart);
    VariableUtils().AddDof<Variable<double>>(RANS_SCALAR_1_ADJOINT_1, rModelPart);
    VariableUtils().AddDof<Variable<double>>(RANS_SCALAR_2_ADJOINT_1, rModelPart);

    VariableUtils().AddDof<Variable<double>>(PRESSURE, rModelPart);
    VariableUtils().AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>(
        VELOCITY_X, rModelPart);
    VariableUtils().AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>(
        VELOCITY_Y, rModelPart);
    VariableUtils().AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>(
        VELOCITY_Z, rModelPart);

    VariableUtils().AddDof<Variable<double>>(ADJOINT_FLUID_SCALAR_1, rModelPart);
    VariableUtils().AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>(
        ADJOINT_FLUID_VECTOR_1_X, rModelPart);
    VariableUtils().AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>(
        ADJOINT_FLUID_VECTOR_1_Y, rModelPart);
    VariableUtils().AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>(
        ADJOINT_FLUID_VECTOR_1_Z, rModelPart);
}

void CreateModelPartElements(ModelPart& rModelPart, const std::string& ElementName)
{
    Properties::Pointer p_elem_prop = rModelPart.pGetProperties(0);

    std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
    rModelPart.CreateNewElement(ElementName,
                                rModelPart.GetRootModelPart().NumberOfElements() + 1,
                                elem_nodes, p_elem_prop);
}

void CreateModelPartConditions(ModelPart& rModelPart, const std::string& ConditionName)
{
    Properties::Pointer p_cond_prop = rModelPart.pGetProperties(0);

    std::vector<ModelPart::IndexType> cond_nodes_1{1, 2};
    rModelPart
        .CreateNewCondition(ConditionName,
                            rModelPart.GetRootModelPart().NumberOfConditions() + 1,
                            cond_nodes_1, p_cond_prop)
        ->Set(SLIP, true);

    std::vector<ModelPart::IndexType> cond_nodes_2{2, 3};
    rModelPart
        .CreateNewCondition(ConditionName,
                            rModelPart.GetRootModelPart().NumberOfConditions() + 1,
                            cond_nodes_2, p_cond_prop)
        ->Set(SLIP, true);

    std::vector<ModelPart::IndexType> cond_nodes_3{3, 1};
    rModelPart
        .CreateNewCondition(ConditionName,
                            rModelPart.GetRootModelPart().NumberOfConditions() + 1,
                            cond_nodes_3, p_cond_prop)
        ->Set(SLIP, false);
}

void InitializeNodalVariables(ModelPart& rModelPart)
{
    using namespace RansModellingApplicationTestUtilities;

    InitializeVariableWithRandomValues(rModelPart, DISTANCE, 1e-2, 1.0, 2);
    InitializeVariableWithRandomValues(rModelPart, TURBULENT_KINETIC_ENERGY, 0.1, 1.0, 2);
    InitializeVariableWithRandomValues(
        rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 5.0, 10.0, 2);
    InitializeVariableWithRandomValues(
        rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 10.0, 20.0, 2);
    InitializeVariableWithRandomValues(
        rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE_2, 5.0, 10.0, 2);
    InitializeVariableWithRandomValues(rModelPart, VELOCITY, 5.0, 10.0, 2);

    InitializeVariableWithRandomValues(rModelPart, PRESSURE, 5.0, 10.0, 2);
    InitializeVariableWithRandomValues(rModelPart, EXTERNAL_PRESSURE, 50.0, 100.0, 2);
    InitializeVariableWithRandomValues(rModelPart, ACCELERATION, 2.0, 3.0, 2);
    InitializeYPlus(rModelPart);

    InitializeVariableWithValues(rModelPart, KINEMATIC_VISCOSITY, 3e-2);
    InitializeVariableWithValues(rModelPart, DENSITY, 200.0);

    InitializeVariableWithRandomValues(rModelPart, ADJOINT_FLUID_VECTOR_1, 55.0, 120.0, 2);
    InitializeVariableWithRandomValues(rModelPart, ADJOINT_FLUID_VECTOR_3, 35.0, 150.0, 2);
    InitializeVariableWithRandomValues(rModelPart, ADJOINT_FLUID_SCALAR_1, 45.0, 100.0, 2);
    InitializeVariableWithRandomValues(rModelPart, AUX_ADJOINT_FLUID_VECTOR_1,
                                       45.0, 100.0, 2);

    InitializeVariableWithRandomValues(rModelPart, RANS_SCALAR_1_ADJOINT_1, 15.0, 20.0, 2);
    InitializeVariableWithRandomValues(rModelPart, RANS_SCALAR_2_ADJOINT_1, 15.0, 20.0, 2);
    InitializeVariableWithRandomValues(rModelPart, RANS_SCALAR_1_ADJOINT_2, 25.0, 30.0, 2);
    InitializeVariableWithRandomValues(rModelPart, RANS_SCALAR_2_ADJOINT_2, 25.0, 30.0, 2);
    InitializeVariableWithRandomValues(rModelPart, RANS_SCALAR_1_ADJOINT_3, 25.0, 30.0, 2);
    InitializeVariableWithRandomValues(rModelPart, RANS_SCALAR_2_ADJOINT_3, 25.0, 30.0, 2);
    InitializeVariableWithRandomValues(rModelPart, RANS_AUX_ADJOINT_SCALAR_1,
                                       25.0, 30.0, 2);
    InitializeVariableWithRandomValues(rModelPart, RANS_AUX_ADJOINT_SCALAR_2,
                                       25.0, 30.0, 2);
}

void InitializeYPlus(ModelPart& rModelPart)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    const double c_mu_25 = std::pow(r_process_info[TURBULENCE_RANS_C_MU], 0.25);
    const double von_karman = r_process_info[WALL_VON_KARMAN];
    const double beta = r_process_info[WALL_SMOOTHNESS_BETA];

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        const array_1d<double, 3>& r_velocity = r_node.FastGetSolutionStepValue(VELOCITY);
        const double velocity_magnitude = norm_2(r_velocity);
        const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);

        double& y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);

        y_plus = std::exp(von_karman *
                          (velocity_magnitude / (c_mu_25 * std::sqrt(tke)) - beta)) *
                 std::pow(-1, i_node) * 1.1;
    }
}

void CreateEquationIds(ModelPart& rModelPart)
{
    std::size_t global_equation_id = 0;
    for (auto& r_node : rModelPart.Nodes())
    {
        for (auto& r_dof : r_node.GetDofs())
        {
            r_dof->SetEquationId(global_equation_id++);
        }
    }
}

template <>
void GenerateRansEvmKEpsilonTestModelPart<ModelPart::ElementsContainerType>(
    ModelPart& rModelPart, const std::string& ElementName)
{
    AddVariablesToModelPart(rModelPart);
    InitializeProcessInfo(rModelPart);
    CreateModelPartNodes(rModelPart);
    CreateModelPartElements(rModelPart, ElementName);
    CreateEquationIds(rModelPart);
    InitializeNodalVariables(rModelPart);
}

template <>
void GenerateRansEvmKEpsilonTestModelPart<ModelPart::ConditionsContainerType>(
    ModelPart& rModelPart, const std::string& ConditionName)
{
    AddVariablesToModelPart(rModelPart);
    InitializeProcessInfo(rModelPart);
    CreateModelPartNodes(rModelPart);
    CreateModelPartConditions(rModelPart, ConditionName);
    CreateEquationIds(rModelPart);
    InitializeNodalVariables(rModelPart);
}

void UpdateVariablesInModelPart(ModelPart& rModelPart)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    const double bossak_alpha = rModelPart.GetProcessInfo()[BOSSAK_ALPHA];

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i_node);

        const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const double& nu_t = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        r_node.FastGetSolutionStepValue(VISCOSITY) = nu + nu_t;

        const double tke_rate =
            r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE);
        const double old_tke_rate =
            r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE, 1);
        const double epsilon_rate =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2);
        const double old_epsilon_rate =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2, 1);
        r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) =
            (1 - bossak_alpha) * tke_rate + bossak_alpha * old_tke_rate;
        r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) =
            (1 - bossak_alpha) * epsilon_rate + bossak_alpha * old_epsilon_rate;

        const array_1d<double, 3>& acceleration =
            r_node.FastGetSolutionStepValue(ACCELERATION);
        const array_1d<double, 3>& old_acceleration =
            r_node.FastGetSolutionStepValue(ACCELERATION, 1);

        r_node.FastGetSolutionStepValue(RELAXED_ACCELERATION) =
            acceleration * (1 - bossak_alpha) + bossak_alpha * old_acceleration;
    }
}

template <typename TDataType, typename TContainer>
void RunRansEvmKEpsilonTest(const std::string& PrimalName,
                            const std::string& AdjointName,
                            const Variable<TDataType>& rPerturbVariable,
                            std::function<void(Matrix&, typename TContainer::data_type&, ProcessInfo&)> CalculateElementResidualScalarSensitivity,
                            const double Delta,
                            const double Tolerance,
                            const int DerivativesOffset,
                            const int EquationOffset)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<TContainer>(
        r_primal_model_part, PrimalName);

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart<TContainer>(
        r_adjoint_model_part, AdjointName);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    std::vector<Process*> primal_processes;
    std::vector<Process*> adjoint_processes;

    primal_processes.push_back(&primal_nut_process);
    adjoint_processes.push_back(&adjoint_nut_process);
    adjoint_processes.push_back(&nut_sensitivities_process);

    auto perturbation_variable =
        RansModellingApplicationTestUtilities::GetPerturbationMethod(rPerturbVariable);

    RansModellingApplicationTestUtilities::RunResidualSensitivityTest<TContainer>(
        r_primal_model_part, r_adjoint_model_part, primal_processes,
        adjoint_processes, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        CalculateElementResidualScalarSensitivity, perturbation_variable, Delta,
        Tolerance, DerivativesOffset, EquationOffset);
}

// templated method instantiation

template void RunRansEvmKEpsilonTest<double, ModelPart::ElementsContainerType>(
    const std::string&,
    const std::string&,
    const Variable<double>&,
    std::function<void(Matrix&, ElementType&, ProcessInfo&)>,
    const double,
    const double,
    const int,
    const int);
template void RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ElementsContainerType>(
    const std::string&,
    const std::string&,
    const Variable<array_1d<double, 3>>&,
    std::function<void(Matrix&, ElementType&, ProcessInfo&)>,
    const double,
    const double,
    const int,
    const int);
template void RunRansEvmKEpsilonTest<double, ModelPart::ConditionsContainerType>(
    const std::string&,
    const std::string&,
    const Variable<double>&,
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)>,
    const double,
    const double,
    const int,
    const int);
template void RunRansEvmKEpsilonTest<array_1d<double, 3>, ModelPart::ConditionsContainerType>(
    const std::string&,
    const std::string&,
    const Variable<array_1d<double, 3>>&,
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)>,
    const double,
    const double,
    const int,
    const int);

} // namespace RansEvmKEpsilonModel
} // namespace Testing
} // namespace Kratos
