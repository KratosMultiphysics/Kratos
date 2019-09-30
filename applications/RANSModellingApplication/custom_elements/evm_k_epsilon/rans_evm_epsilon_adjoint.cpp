//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    , KratosAppGenerator
//

// System includes

// External includes

// Include Base h
#include "rans_evm_epsilon_adjoint.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_adjoint_utilities.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "rans_modelling_application_variables.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * Constructor.
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmEpsilonAdjoint<TDim, TNumNodes>::RansEvmEpsilonAdjoint(IndexType NewId)
    : BaseType(NewId)
{
}

/**
 * Constructor using an array of nodes
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmEpsilonAdjoint<TDim, TNumNodes>::RansEvmEpsilonAdjoint(IndexType NewId,
                                                                    const NodesArrayType& ThisNodes)
    : BaseType(NewId, ThisNodes)
{
}

/**
 * Constructor using Geometry
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmEpsilonAdjoint<TDim, TNumNodes>::RansEvmEpsilonAdjoint(IndexType NewId,
                                                                    GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmEpsilonAdjoint<TDim, TNumNodes>::RansEvmEpsilonAdjoint(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

/**
 * Copy Constructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmEpsilonAdjoint<TDim, TNumNodes>::RansEvmEpsilonAdjoint(
    RansEvmEpsilonAdjoint<TDim, TNumNodes> const& rOther)
    : BaseType(rOther)
{
}

/**
 * Destructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmEpsilonAdjoint<TDim, TNumNodes>::~RansEvmEpsilonAdjoint()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmEpsilonAdjoint<TDim, TNumNodes>& RansEvmEpsilonAdjoint<TDim, TNumNodes>::operator=(
    RansEvmEpsilonAdjoint<TDim, TNumNodes> const& rOther)
{
    BaseType::operator=(rOther);
    Flags::operator=(rOther);
    // mpProperties = rOther.mpProperties;
    return *this;
}

///@}
///@name Operations
///@{

/**
 * ELEMENTS inherited from this class have to implement next
 * Create and Clone methods: MANDATORY
 */

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer RansEvmEpsilonAdjoint<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmEpsilonAdjoint>(
        NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer RansEvmEpsilonAdjoint<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmEpsilonAdjoint>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer RansEvmEpsilonAdjoint<TDim, TNumNodes>::Clone(IndexType NewId,
                                                                  NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmEpsilonAdjoint>(
        NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod RansEvmEpsilonAdjoint<TDim, TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

/**
 * This method provides the place to perform checks on the completeness of the
 * input and the compatibility with the problem options as well as the
 * contitutive laws selected It is designed to be called only once (or anyway,
 * not often) typically at the beginning of the calculations, so to verify that
 * nothing is missing from the input or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
template <unsigned int TDim, unsigned int TNumNodes>
int RansEvmEpsilonAdjoint<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Check(rCurrentProcessInfo);

    KRATOS_CHECK_VARIABLE_KEY(TURBULENCE_RANS_C_MU);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENCE_RANS_C1);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENCE_RANS_C2);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY_SIGMA);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE);
    KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
    KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS);
    KRATOS_CHECK_VARIABLE_KEY(RANS_AUXILIARY_VARIABLE_2);
    KRATOS_CHECK_VARIABLE_KEY(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE);
    KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS_VELOCITY_DERIVATIVES);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(RANS_SCALAR_2_ADJOINT_1);
    KRATOS_CHECK_VARIABLE_KEY(RANS_SCALAR_2_ADJOINT_3);

    for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
    {
        NodeType& r_node = this->GetGeometry()[iNode];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_3, r_node);

        KRATOS_CHECK_DOF_IN_NODE(RANS_SCALAR_2_ADJOINT_1, r_node);
    }
    return 0;

    KRATOS_CATCH("");
}

///@}
///@name Access
///@{

///@}
///@name Inquiry
///@{

///@}
///@name Input and output
///@{

/// Turn back information as a string.

template <unsigned int TDim, unsigned int TNumNodes>
std::string RansEvmEpsilonAdjoint<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "RansEvmEpsilonAdjoint #" << Element::Id();
    return buffer.str();
}

/// Print information about this object.

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "RansEvmEpsilonAdjoint #" << Element::Id();
}

/// Print object's data.

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
    Element::pGetGeometry()->PrintData(rOStream);
}

///@}
///@name Friends
///@{

///@}

///@name Protected static Member Variables
///@{

///@}
///@name Protected member Variables
///@{

///@}
///@name Protected Operators
///@{

///@}
///@name Protected Operations
///@{

///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{

///@}
///@name Protected LifeCycle
///@{

///@}

///@name Static Member Variables
///@{

///@}
///@name Member Variables
///@{

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{
template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& RansEvmEpsilonAdjoint<TDim, TNumNodes>::GetPrimalVariable() const
{
    return TURBULENT_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& RansEvmEpsilonAdjoint<TDim, TNumNodes>::GetPrimalRelaxedRateVariable() const
{
    return RANS_AUXILIARY_VARIABLE_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& RansEvmEpsilonAdjoint<TDim, TNumNodes>::GetAdjointVariable() const
{
    return RANS_SCALAR_2_ADJOINT_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& RansEvmEpsilonAdjoint<TDim, TNumNodes>::GetAdjointSecondVariable() const
{
    return RANS_SCALAR_2_ADJOINT_3;
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateElementData(
    RansEvmEpsilonAdjointData& rData,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rData.ShapeFunctionDerivatives = rShapeFunctionDerivatives;
    rData.ShapeFunctions = rShapeFunctions;

    const double& c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    const double& epsilon_sigma =
        rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    const double& c1 = rCurrentProcessInfo[TURBULENCE_RANS_C1];
    const double& c2 = rCurrentProcessInfo[TURBULENCE_RANS_C2];

    const double& nu_t = this->EvaluateInPoint(TURBULENT_VISCOSITY, rShapeFunctions);
    const double& tke = this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, rShapeFunctions);
    const double& epsilon =
        this->EvaluateInPoint(TURBULENT_ENERGY_DISSIPATION_RATE, rShapeFunctions);
    const double& nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, rShapeFunctions);
    const double& wall_distance = this->EvaluateInPoint(DISTANCE, rShapeFunctions);
    const double& y_plus = this->EvaluateInPoint(RANS_Y_PLUS, rShapeFunctions);
    const double& gamma = EvmKepsilonModelUtilities::CalculateGamma(c_mu, 1.0, tke, nu_t);

    rData.TurbulentKinematicViscositySensitivitiesK.resize(TNumNodes);
    rData.TurbulentKinematicViscositySensitivitiesEpsilon.resize(TNumNodes);

    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const NodeType& r_node = this->GetGeometry()[i_node];
        const Vector& turbulent_kinematic_viscosity_sensitivities =
            r_node.GetValue(RANS_NUT_PARTIAL_DERIVATIVES);
        rData.TurbulentKinematicViscositySensitivitiesK[i_node] =
            turbulent_kinematic_viscosity_sensitivities[0];
        rData.TurbulentKinematicViscositySensitivitiesEpsilon[i_node] =
            turbulent_kinematic_viscosity_sensitivities[1];
    }

    rData.TurbulentKinematicViscosity = nu_t;
    rData.TurbulentKineticEnergy = tke;
    rData.TurbulentEnergyDissipationRate = epsilon;
    rData.KinematicViscosity = nu;
    rData.WallDistance = wall_distance;
    rData.Gamma = gamma;
    rData.EffectiveKinematicViscosity = nu + nu_t / epsilon_sigma;
    rData.y_plus = y_plus;
    rData.C1 = c1;
    rData.C2 = c2;
    rData.VelocityDivergence =
        this->GetDivergenceOperator(VELOCITY, rShapeFunctionDerivatives);

    RansVariableUtils rans_variable_utils;

    rans_variable_utils.GetNodalArray(rData.NodalTurbulentKineticEnergy, *this,
                                      TURBULENT_KINETIC_ENERGY);
    rans_variable_utils.GetNodalArray(rData.NodalTurbulentEnergyDissipationRate,
                                      *this, TURBULENT_ENERGY_DISSIPATION_RATE);
    rans_variable_utils.GetNodalArray(rData.NodalYPlus, *this, RANS_Y_PLUS);

    std::size_t number_of_nodes = rData.NodalYPlus.size();

    if (rData.NodalFmu.size() != number_of_nodes)
        rData.NodalFmu.resize(rData.NodalYPlus.size());

    if (rData.NodalVelocity.size1() != TNumNodes || rData.NodalVelocity.size2() != TDim)
        rData.NodalVelocity.resize(TNumNodes, TDim);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        rData.NodalFmu[i_node] = 1.0;
        const array_1d<double, 3>& rVelocity =
            this->GetGeometry()[i_node].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
        {
            rData.NodalVelocity(i_node, i_dim) = rVelocity[i_dim];
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateEffectiveKinematicViscosity(
    const RansEvmEpsilonAdjointData& rCurrentData, const ProcessInfo& rCurrentProcessInfo) const
{
    return rCurrentData.EffectiveKinematicViscosity;
}

template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateReactionTerm(
    const RansEvmEpsilonAdjointData& rData, const ProcessInfo& rCurrentProcessInfo) const
{
    return rData.C2 * rData.Gamma + rData.C1 * 2.0 * rData.VelocityDivergence / 3.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateSourceTerm(
    const RansEvmEpsilonAdjointData& rData, const ProcessInfo& rCurrentProcessInfo) const
{
    BoundedMatrix<double, TDim, TDim> velocity_gradient_matrix;
    this->CalculateGradient(velocity_gradient_matrix, VELOCITY, rData.ShapeFunctionDerivatives);
    double production = EvmKepsilonModelUtilities::CalculateSourceTerm<TDim>(
        velocity_gradient_matrix, rData.TurbulentKinematicViscosity, rData.TurbulentKineticEnergy);

    production *= (rData.C1 * rData.Gamma);

    return production;
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityScalarDerivatives(
    Vector& rOutput,
    const Variable<double>& rDerivativeVariable,
    const RansEvmEpsilonAdjointData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        const double epsilon_sigma =
            rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
        // EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityTKESensitivities(
        //     rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
        //     rCurrentData.NodalTurbulentEnergyDissipationRate, rCurrentData.NodalFmu);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rCurrentData.TurbulentKinematicViscositySensitivitiesK,
            rCurrentData.ShapeFunctions);

        noalias(rOutput) = rOutput / epsilon_sigma;
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        const double epsilon_sigma =
            rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
        // const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        // EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityEpsilonSensitivities(
        //     rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
        //     rCurrentData.NodalTurbulentEnergyDissipationRate, rCurrentData.NodalFmu);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rCurrentData.TurbulentKinematicViscositySensitivitiesEpsilon,
            rCurrentData.ShapeFunctions);

        noalias(rOutput) = rOutput / epsilon_sigma;
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name() << " used in RansEvmEpsilonAdjoint::CalculateEffectiveKinematicViscosityDerivatives method.";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateReactionTermScalarDerivatives(
    Vector& rOutput,
    const Variable<double>& rDerivativeVariable,
    const RansEvmEpsilonAdjointData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
        const double c2 = rCurrentData.C2;

        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rCurrentData.TurbulentKinematicViscositySensitivitiesK,
            rCurrentData.ShapeFunctions);

        Vector theta_sensitivities(rOutput.size());
        EvmKepsilonModelAdjointUtilities::CalculateThetaTKESensitivity(
            theta_sensitivities, c_mu, 1.0, rCurrentData.TurbulentKineticEnergy,
            rCurrentData.TurbulentKinematicViscosity, rOutput, rCurrentData.ShapeFunctions);

        noalias(rOutput) = theta_sensitivities * c2;
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
        const double c2 = rCurrentData.C2;

        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rCurrentData.TurbulentKinematicViscositySensitivitiesEpsilon,
            rCurrentData.ShapeFunctions);

        Vector theta_sensitivities(rOutput.size());
        EvmKepsilonModelAdjointUtilities::CalculateThetaEpsilonSensitivity(
            theta_sensitivities, c_mu, 1.0, rCurrentData.TurbulentKineticEnergy,
            rCurrentData.TurbulentKinematicViscosity, rOutput);

        noalias(rOutput) = theta_sensitivities * c2;
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name() << " used in RansEvmEpsilonAdjoint::CalculateReactionTermScalarDerivatives method.";
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateSourceTermScalarDerivatives(
    Vector& rOutput,
    const Variable<double>& rDerivativeVariable,
    const RansEvmEpsilonAdjointData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
        const double c1 = rCurrentData.C1;

        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rCurrentData.TurbulentKinematicViscositySensitivitiesK,
            rCurrentData.ShapeFunctions);

        BoundedMatrix<double, TDim, TDim> velocity_gradient;
        this->CalculateGradient(velocity_gradient, VELOCITY,
                                rCurrentData.ShapeFunctionDerivatives);

        const double p_k = EvmKepsilonModelUtilities::CalculateSourceTerm<TDim>(
            velocity_gradient, rCurrentData.TurbulentKinematicViscosity,
            rCurrentData.TurbulentKineticEnergy);

        Vector p_k_sensitivities(rOutput.size());
        EvmKepsilonModelAdjointUtilities::CalculateProductionScalarSensitivities<TDim>(
            p_k_sensitivities, rOutput, velocity_gradient);

        Vector theta_sensitivities(rOutput.size());
        EvmKepsilonModelAdjointUtilities::CalculateThetaTKESensitivity(
            theta_sensitivities, c_mu, 1.0, rCurrentData.TurbulentKineticEnergy,
            rCurrentData.TurbulentKinematicViscosity, rOutput, rCurrentData.ShapeFunctions);

        noalias(rOutput) = theta_sensitivities * (p_k * c1) +
                           p_k_sensitivities * (rCurrentData.Gamma * c1);
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
        const double c1 = rCurrentData.C1;

        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rCurrentData.TurbulentKinematicViscositySensitivitiesEpsilon,
            rCurrentData.ShapeFunctions);

        BoundedMatrix<double, TDim, TDim> velocity_gradient;
        this->CalculateGradient(velocity_gradient, VELOCITY,
                                rCurrentData.ShapeFunctionDerivatives);

        const double p_k = EvmKepsilonModelUtilities::CalculateSourceTerm<TDim>(
            velocity_gradient, rCurrentData.TurbulentKinematicViscosity,
            rCurrentData.TurbulentKineticEnergy);

        Vector p_k_sensitivities(rOutput.size());
        EvmKepsilonModelAdjointUtilities::CalculateProductionScalarSensitivities<TDim>(
            p_k_sensitivities, rOutput, velocity_gradient);

        Vector theta_sensitivities(rOutput.size());
        EvmKepsilonModelAdjointUtilities::CalculateThetaEpsilonSensitivity(
            theta_sensitivities, c_mu, 1.0, rCurrentData.TurbulentKineticEnergy,
            rCurrentData.TurbulentKinematicViscosity, rOutput);

        noalias(rOutput) = theta_sensitivities * (p_k * c1) +
                           p_k_sensitivities * (rCurrentData.Gamma * c1);
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name() << " used in RansEvmEpsilonAdjoint::CalculateSourceTermScalarDerivatives method.";
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::Calculate(const Variable<Matrix>& rVariable,
                                                          Matrix& rOutput,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE)
    {
        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculateResidualScalarDerivatives(
            TURBULENT_KINETIC_ENERGY, local_matrix, rCurrentProcessInfo);
        if (rOutput.size1() != local_matrix.size1() ||
            rOutput.size2() != local_matrix.size2())
            rOutput.resize(local_matrix.size1(), local_matrix.size2(), false);

        noalias(rOutput) = local_matrix;
    }
    else
    {
        BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityVelocityDerivatives(
    Matrix& rOutput,
    const RansEvmEpsilonAdjointData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateReactionTermVelocityDerivatives(
    Matrix& rOutput,
    const RansEvmEpsilonAdjointData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    noalias(rOutput) =
        rCurrentData.ShapeFunctionDerivatives * ((2.0 / 3.0) * rCurrentData.C1);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateSourceTermVelocityDerivatives(
    Matrix& rOutput,
    const RansEvmEpsilonAdjointData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const double c1 = rCurrentData.C1;

    BoundedMatrix<double, TDim, TDim> velocity_gradient;
    this->CalculateGradient(velocity_gradient, VELOCITY, rCurrentData.ShapeFunctionDerivatives);

    EvmKepsilonModelAdjointUtilities::CalculateProductionVelocitySensitivities<TDim>(
        rOutput, rCurrentData.TurbulentKinematicViscosity,
        ZeroMatrix(rOutput.size1(), rOutput.size2()), velocity_gradient,
        rCurrentData.ShapeFunctionDerivatives);

    noalias(rOutput) = rOutput * (c1 * rCurrentData.Gamma);
}

template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityShapeSensitivity(
    const RansEvmEpsilonAdjointData& rCurrentData,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
    const ProcessInfo& rCurrentProcessInfo) const
{
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateReactionTermShapeSensitivity(
    const RansEvmEpsilonAdjointData& rCurrentData,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
    const ProcessInfo& rCurrentProcessInfo) const
{
    return (2.0 * rCurrentData.C1 / 3.0) * this->GetDivergenceOperator(VELOCITY, rDN_Dx_deriv);
}

template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmEpsilonAdjoint<TDim, TNumNodes>::CalculateSourceTermShapeSensitivity(
    const RansEvmEpsilonAdjointData& rCurrentData,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
    const ProcessInfo& rCurrentProcessInfo) const
{
    double value;

    const double c1 = rCurrentData.C1;

    BoundedMatrix<double, TDim, TDim> velocity_gradient;
    this->CalculateGradient(velocity_gradient, VELOCITY, rCurrentData.ShapeFunctionDerivatives);

    EvmKepsilonModelAdjointUtilities::CalculateProductionShapeSensitivities<TDim>(
        value, rCurrentData.TurbulentKinematicViscosity, 0.0, rCurrentData.NodalVelocity,
        velocity_gradient, rCurrentData.ShapeFunctionDerivatives, rDN_Dx_deriv);

    return value * c1 * rCurrentData.Gamma;
}

///@}
///@name Serialization
///@{

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmEpsilonAdjoint<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}

///@}
///@name Private  Access
///@{

///@}
///@name Private Inquiry
///@{

///@}
///@name Un accessible methods
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function

template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                RansEvmEpsilonAdjoint<TDim, TNumNodes>& rThis);

/// output stream function

template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEvmEpsilonAdjoint<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

// Class template instantiation
// Epsilon only
template class RansEvmEpsilonAdjoint<2, 3>;
template class RansEvmEpsilonAdjoint<3, 4>;
template class RansEvmEpsilonAdjoint<2, 4>;
template class RansEvmEpsilonAdjoint<3, 8>;

} // namespace Kratos.
