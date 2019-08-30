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
#include "evm_k_epsilon_vms_adjoint_element.h"

#include "custom_elements/evm_k_epsilon/evm_k_epsilon_adjoint_utilities.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "includes/cfd_variables.h"
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
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::EvmKEpsilonVMSAdjointElement(IndexType NewId)
    : BaseType(NewId)
{
}

/**
 * Constructor using Geometry
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::EvmKEpsilonVMSAdjointElement(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::EvmKEpsilonVMSAdjointElement(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

/**
 * Destructor
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::~EvmKEpsilonVMSAdjointElement()
{
}

///@}
///@name Operators
///@{

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
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
Element::Pointer EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EvmKEpsilonVMSAdjointElement>(
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
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
Element::Pointer EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EvmKEpsilonVMSAdjointElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
Element::Pointer EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EvmKEpsilonVMSAdjointElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
    KRATOS_CATCH("");
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
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
int EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::Check(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Check(rCurrentProcessInfo);

    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE);
    KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENCE_RANS_C_MU);
    KRATOS_CHECK_VARIABLE_KEY(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE);
    KRATOS_CHECK_VARIABLE_KEY(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE);
    KRATOS_CHECK_VARIABLE_KEY(RELAXED_ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS_VELOCITY_DERIVATIVES);

    for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
    {
        NodeType& r_node = this->GetGeometry()[iNode];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RELAXED_ACCELERATION, r_node);
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

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
std::string EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::Info() const
{
    std::stringstream buffer;
    buffer << "EvmKEpsilonVMSAdjointElement #" << Element::Id();
    return buffer.str();
}

/// Print information about this object.

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
void EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::PrintInfo(
    std::ostream& rOStream) const
{
    rOStream << "EvmKEpsilonVMSAdjointElement #" << Element::Id();
}

/// Print object's data.

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
void EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::PrintData(
    std::ostream& rOStream) const
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

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
void EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::CalculateElementData(
    RANSEvmVMSAdjointElementData& rData,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rData.ShapeFunctionDerivatives = rShapeFunctionDerivatives;
    rData.ShapeFunctions = rShapeFunctions;

    RansVariableUtils rans_variable_utils;

    rans_variable_utils.GetNodalArray(rData.NodalTurbulentKineticEnergy, *this,
                                      TURBULENT_KINETIC_ENERGY);
    rans_variable_utils.GetNodalArray(rData.NodalTurbulentEnergyDissipationRate,
                                      *this, TURBULENT_ENERGY_DISSIPATION_RATE);
    rans_variable_utils.GetNodalArray(rData.NodalYPlus, *this, RANS_Y_PLUS);

    const std::size_t number_of_nodes = rData.NodalYPlus.size();

    if (rData.NodalFmu.size() != number_of_nodes)
        rData.NodalFmu.resize(rData.NodalYPlus.size());

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        rData.NodalFmu[i_node] =
            1.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
void EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::CalculateTurbulentKinematicViscosityScalarDerivatives(
    Vector& rOutput,
    const Variable<double>& rDerivativeVariable,
    const RANSEvmVMSAdjointElementData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityTKESensitivities(
            rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
            rCurrentData.NodalTurbulentEnergyDissipationRate, rCurrentData.NodalFmu);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rOutput, rCurrentData.ShapeFunctions);
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityEpsilonSensitivities(
            rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
            rCurrentData.NodalTurbulentEnergyDissipationRate, rCurrentData.NodalFmu);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rOutput, rCurrentData.ShapeFunctions);
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name() << " used in EvmKAdjointElement::CalculateEffectiveKinematicViscosityDerivatives method.";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
void EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::Calculate(
    const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE)
    {
        this->CalculateElementTotalSteadyResidualScalarDerivatives(
            Output, TURBULENT_ENERGY_DISSIPATION_RATE, rCurrentProcessInfo);
        this->AddElementTotalMassResidualScalarDerivatives(
            Output, RELAXED_ACCELERATION, TURBULENT_ENERGY_DISSIPATION_RATE,
            -1.0, rCurrentProcessInfo);
    }
    else if (rVariable == RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE)
    {
        this->CalculateElementTotalSteadyResidualScalarDerivatives(
            Output, TURBULENT_KINETIC_ENERGY, rCurrentProcessInfo);
        this->AddElementTotalMassResidualScalarDerivatives(
            Output, RELAXED_ACCELERATION, TURBULENT_KINETIC_ENERGY, -1.0, rCurrentProcessInfo);
    }
    else
    {
        KRATOS_ERROR << "Unsupported variable "
                     << rVariable.Name() << " requested at StabilizedConvectionDiffusionReactionAdjoint::Calculate.";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
void EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::CalculateTurbulentKinematicViscosityVelocityDerivatives(
    Matrix& rOutput,
    const RANSEvmVMSAdjointElementData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rOutput.clear();
}

///@}
///@name Serialization
///@{

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
void EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
void EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>::load(Serializer& rSerializer)
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

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
inline std::istream& operator>>(
    std::istream& rIStream,
    EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>& rThis);

/// output stream function

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TMonolithicAssemblyNodalDofSize>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const EvmKEpsilonVMSAdjointElement<TDim, TNumNodes, TMonolithicAssemblyNodalDofSize>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

// Class template instantiation

template class EvmKEpsilonVMSAdjointElement<2>;
template class EvmKEpsilonVMSAdjointElement<3>;

template class EvmKEpsilonVMSAdjointElement<2, 3, 5>;
template class EvmKEpsilonVMSAdjointElement<3, 4, 6>;

} // namespace Kratos.
