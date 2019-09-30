//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes

// Include Base h
#include "rans_evm_k_element.h"
#include "custom_utilities/rans_calculation_utilities.h"

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
RansEvmKElement<TDim, TNumNodes>::RansEvmKElement(IndexType NewId)
    : StabilizedConvectionDiffusionReactionElement<TDim, TNumNodes, RansEvmKElementData>(NewId)
{
}

/**
 * Constructor using an array of nodes
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmKElement<TDim, TNumNodes>::RansEvmKElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : StabilizedConvectionDiffusionReactionElement<TDim, TNumNodes, RansEvmKElementData>(
          NewId, ThisNodes)
{
}

/**
 * Constructor using Geometry
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmKElement<TDim, TNumNodes>::RansEvmKElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : StabilizedConvectionDiffusionReactionElement<TDim, TNumNodes, RansEvmKElementData>(
          NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmKElement<TDim, TNumNodes>::RansEvmKElement(IndexType NewId,
                                                  GeometryType::Pointer pGeometry,
                                                  PropertiesType::Pointer pProperties)
    : StabilizedConvectionDiffusionReactionElement<TDim, TNumNodes, RansEvmKElementData>(
          NewId, pGeometry, pProperties)
{
}

/**
 * Copy Constructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmKElement<TDim, TNumNodes>::RansEvmKElement(RansEvmKElement<TDim, TNumNodes> const& rOther)
    : StabilizedConvectionDiffusionReactionElement<TDim, TNumNodes, RansEvmKElementData>(rOther)
{
}

/**
 * Destructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmKElement<TDim, TNumNodes>::~RansEvmKElement()
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
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer RansEvmKElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                          NodesArrayType const& ThisNodes,
                                                          PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmKElement>(
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
Element::Pointer RansEvmKElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                          GeometryType::Pointer pGeom,
                                                          PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmKElement>(NewId, pGeom, pProperties);
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
Element::Pointer RansEvmKElement<TDim, TNumNodes>::Clone(IndexType NewId,
                                                         NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmKElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                        ProcessInfo& CurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (unsigned int i = 0; i < TNumNodes; i++)
        rResult[i] = Element::GetGeometry()[i].GetDof(TURBULENT_KINETIC_ENERGY).EquationId();
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                  ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != TNumNodes)
        rElementalDofList.resize(TNumNodes);

    for (unsigned int i = 0; i < TNumNodes; i++)
        rElementalDofList[i] = Element::GetGeometry()[i].pGetDof(TURBULENT_KINETIC_ENERGY);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKElement<TDim, TNumNodes>::GetValuesVector(VectorType& rValues, int Step)
{
    this->GetFirstDerivativesVector(rValues, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKElement<TDim, TNumNodes>::GetFirstDerivativesVector(VectorType& rValues, int Step)
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] =
            rGeom[iNode].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKElement<TDim, TNumNodes>::GetSecondDerivativesVector(VectorType& rValues, int Step)
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] =
            rGeom[iNode].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod RansEvmKElement<TDim, TNumNodes>::GetIntegrationMethod() const
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
int RansEvmKElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Check(rCurrentProcessInfo);

    for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
    {
        NodeType& r_node = this->GetGeometry()[iNode];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_1, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_KINETIC_ENERGY, r_node);
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
std::string RansEvmKElement<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "RansEvmKElement #" << Element::Id();
    return buffer.str();
}

/// Print information about this object.

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKElement<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "RansEvmKElement #" << Element::Id();
}

/// Print object's data.

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKElement<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
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
void RansEvmKElement<TDim, TNumNodes>::CalculateConvectionDiffusionReactionData(
    RansEvmKElementData& rData,
    double& rEffectiveKinematicViscosity,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo,
    const int Step) const
{
    KRATOS_TRY

    rData.ShapeFunctionDerivatives = rShapeFunctionDerivatives;

    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    const double tke_sigma = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];

    const double nu_t = this->EvaluateInPoint(TURBULENT_VISCOSITY, rShapeFunctions);
    const double tke = this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, rShapeFunctions);
    const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, rShapeFunctions);
    const double gamma = EvmKepsilonModelUtilities::CalculateGamma(c_mu, 1.0, tke, nu_t);

    rData.TurbulentKinematicViscosity = nu_t;
    rData.TurbulentKineticEnergy = tke;
    rData.Gamma = gamma;
    rData.VelocityDivergence =
        this->GetDivergenceOperator(VELOCITY, rShapeFunctionDerivatives);
    rEffectiveKinematicViscosity = nu + nu_t / tke_sigma;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKElement<TDim, TNumNodes>::CalculateConvectionDiffusionReactionData(
    RansEvmKElementData& rData,
    double& rEffectiveKinematicViscosity,
    double& rVariableGradientNorm,
    double& rVariableRelaxedAcceleration,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo,
    const int Step) const
{
    this->CalculateConvectionDiffusionReactionData(
        rData, rEffectiveKinematicViscosity, rShapeFunctions,
        rShapeFunctionDerivatives, rCurrentProcessInfo, Step);

    array_1d<double, 3> tke_gradient;
    this->CalculateGradient(tke_gradient, TURBULENT_KINETIC_ENERGY, rShapeFunctionDerivatives);
    rVariableGradientNorm = norm_2(tke_gradient);

    rVariableRelaxedAcceleration =
        this->EvaluateInPoint(RANS_AUXILIARY_VARIABLE_1, rShapeFunctions);
}

template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmKElement<TDim, TNumNodes>::CalculateReactionTerm(const RansEvmKElementData& rData,
                                                               const ProcessInfo& rCurrentProcessInfo,
                                                               const int Step) const
{
    return std::max(rData.Gamma + (2.0 / 3.0) * rData.VelocityDivergence, 0.0);
}

template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmKElement<TDim, TNumNodes>::CalculateSourceTerm(const RansEvmKElementData& rData,
                                                             const ProcessInfo& rCurrentProcessInfo,
                                                             const int Step) const
{
    double production = 0.0;

    BoundedMatrix<double, TDim, TDim> velocity_gradient_matrix;
    this->CalculateGradient(velocity_gradient_matrix, VELOCITY, rData.ShapeFunctionDerivatives);

    production = EvmKepsilonModelUtilities::CalculateSourceTerm<TDim>(
        velocity_gradient_matrix, rData.TurbulentKinematicViscosity, rData.TurbulentKineticEnergy);

    return production;
}

///@}
///@name Serialization
///@{

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKElement<TDim, TNumNodes>::load(Serializer& rSerializer)
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
                                RansEvmKElement<TDim, TNumNodes>& rThis);

/// output stream function

template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEvmKElement<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

// Class template instantiation

template class RansEvmKElement<2, 3>;
template class RansEvmKElement<3, 4>;
template class RansEvmKElement<2, 4>;
template class RansEvmKElement<3, 8>;

} // namespace Kratos.
