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
#include "includes/cfd_variables.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "custom_utilities/test_utilities.h"
#include "rans_modelling_application_variables.h"

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

void CreateModelPartElements(ModelPart& rModelPart, std::string ElementName)
{
    Properties::Pointer p_elem_prop = rModelPart.pGetProperties(0);

    std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
    rModelPart.CreateNewElement(ElementName,
                                rModelPart.GetRootModelPart().NumberOfElements() + 1,
                                elem_nodes, p_elem_prop);
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
    InitializeVariableWithRandomValues(rModelPart, ACCELERATION, 2.0, 3.0, 2);

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

void GenerateRansEvmKEpsilonTestModelPart(ModelPart& rModelPart, std::string ElementName)
{
    AddVariablesToModelPart(rModelPart);
    InitializeProcessInfo(rModelPart);
    CreateModelPartNodes(rModelPart);
    CreateModelPartElements(rModelPart, ElementName);
    InitializeNodalVariables(rModelPart);
}

void UpdateVariablesInModelPart(ModelPart& rModelPart)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    const double c_mu = rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU];
    const double bossak_alpha = rModelPart.GetProcessInfo()[BOSSAK_ALPHA];

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        const double epsilon =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);

        double& nu_t = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            c_mu, tke, epsilon, 1.0);
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

void UpdateVariablesInModelPartLowRe(ModelPart& rModelPart)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    const double c_mu = rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU];
    const double bossak_alpha = rModelPart.GetProcessInfo()[BOSSAK_ALPHA];

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
        const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        const double epsilon =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
        const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);

        double& nu_t = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            c_mu, tke, epsilon, f_mu);
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

void CalculatePrimalQuantities(std::vector<double>& rValues,
                               const ElementType& rElement,
                               const Vector& rGaussShapeFunctions,
                               const Matrix& rGaussShapeFunctionDerivatives,
                               const ProcessInfo& rCurrentProcessInfo)
{
    RansCalculationUtilities rans_calculation_utilities;

    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

    const GeometryType& r_geometry = rElement.GetGeometry();

    const double y_plus = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, RANS_Y_PLUS, rGaussShapeFunctions);
    const double tke = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, TURBULENT_KINETIC_ENERGY, rGaussShapeFunctions);
    const double epsilon = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, rGaussShapeFunctions);
    const double nu_t = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, TURBULENT_VISCOSITY, rGaussShapeFunctions);
    const double nu = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, KINEMATIC_VISCOSITY, rGaussShapeFunctions);

    const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
    const double Re_t = std::pow(tke, 2) / (nu * epsilon);
    const double theta = EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);
    const double f2 = EvmKepsilonModelUtilities::CalculateF2(tke, nu, epsilon);

    BoundedMatrix<double, 2, 2> velocity_gradient_matrix;
    rans_calculation_utilities.CalculateGradient<2>(
        velocity_gradient_matrix, r_geometry, VELOCITY, rGaussShapeFunctionDerivatives);
    const double P_k = EvmKepsilonModelUtilities::CalculateSourceTerm<2>(
        velocity_gradient_matrix, nu_t, tke);

    rValues.clear();
    rValues.push_back(nu_t);
    rValues.push_back(P_k);
    rValues.push_back(theta);
    rValues.push_back(Re_t);
    rValues.push_back(f2);
    rValues.push_back(f_mu);
    rValues.push_back(nu);
    rValues.push_back(epsilon);
    rValues.push_back(tke);
    rValues.push_back(y_plus);
}

void ReadNodalDataFromElement(Vector& rYPlus,
                              Vector& rTKE,
                              Vector& rEpsilon,
                              Vector& rNut,
                              Vector& rFmu,
                              const Element& rElement)
{
    const auto& r_geometry = rElement.GetGeometry();
    std::size_t number_of_nodes = r_geometry.PointsNumber();

    RansVariableUtils rans_variable_utils;

    rans_variable_utils.GetNodalArray(rTKE, rElement, TURBULENT_KINETIC_ENERGY);
    rans_variable_utils.GetNodalArray(rEpsilon, rElement, TURBULENT_ENERGY_DISSIPATION_RATE);
    rans_variable_utils.GetNodalArray(rYPlus, rElement, RANS_Y_PLUS);
    rans_variable_utils.GetNodalArray(rNut, rElement, TURBULENT_VISCOSITY);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
        rFmu[i_node] = EvmKepsilonModelUtilities::CalculateFmu(rYPlus[i_node]);
}
} // namespace RansEvmKEpsilonModel
} // namespace Testing
} // namespace Kratos
