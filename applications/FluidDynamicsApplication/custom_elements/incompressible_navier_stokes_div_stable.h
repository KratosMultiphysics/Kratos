//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#pragma once

// System includes


// External indludes


// Project includes
#include "geometries/geometry.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"

// Application includes
#include "custom_elements/fluid_element.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

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

template< unsigned int TDim >
class IncompressibleNavierStokesDivStable : public Element
{
public:
    ///@name Type Definitions
    ///@{

    static constexpr std::size_t StrainSize = 3*(TDim-1);

    static constexpr std::size_t VelocityNumNodes = TDim == 2 ? 6 : 10;

    static constexpr std::size_t PressureNumNodes = TDim == 2 ? 3 : 4;

    static constexpr std::size_t LocalSize = VelocityNumNodes*TDim + PressureNumNodes;

    static constexpr IntegrationMethod IntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(IncompressibleNavierStokesDivStable);

    using BaseType = Element;

    using NodeType = BaseType::NodeType;

    using GeometryType = BaseType::GeometryType;

    using NodesArrayType = BaseType::NodesArrayType;

    using VectorType = BaseType::VectorType;

    using MatrixType = BaseType::MatrixType;

    using IndexType = BaseType::IndexType;

    using SizeType = BaseType::SizeType;

    using EquationIdVectorType = BaseType::EquationIdVectorType;

    using DofsVectorType = BaseType::DofsVectorType;

    using DofsArrayType = BaseType::DofsArrayType;

    struct ElementDataContainer
    {
        // Gauss point kinematics
        array_1d<double, VelocityNumNodes> N_v;
        BoundedMatrix<double, VelocityNumNodes, TDim> DN_v;
        GeometryType::ShapeFunctionsSecondDerivativesType DDN_v;

        array_1d<double, PressureNumNodes> N_p;
        BoundedMatrix<double, PressureNumNodes, TDim> DN_p;

        double N_e;
        BoundedMatrix<double, 1, TDim> DN_e;

        // Nodal values
        array_1d<double, PressureNumNodes> Pressure;
        BoundedMatrix<double, VelocityNumNodes, TDim> Velocity;
        BoundedMatrix<double, VelocityNumNodes, TDim> VelocityOld1;
        BoundedMatrix<double, VelocityNumNodes, TDim> VelocityOld2;
        BoundedMatrix<double, VelocityNumNodes, TDim> MeshVelocity;
        BoundedMatrix<double, VelocityNumNodes, TDim> BodyForce;

        // Material response variables
        Vector StrainRate = ZeroVector(StrainSize);
        Vector ShearStress = ZeroVector(StrainSize);
        Matrix ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);

        // Auxiliary values
        double BDF0;
        double BDF1;
        double BDF2;
        double Weight;
        double DeltaTime;

        // Stabilization values
        double StabC1;
        double StabC2;
        double DynamicTau;
        double ElementSize;

        // Material parameters
        double Density;
        double EffectiveViscosity;
    };

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    IncompressibleNavierStokesDivStable(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    IncompressibleNavierStokesDivStable(
        IndexType NewId,
        const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    IncompressibleNavierStokesDivStable(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    IncompressibleNavierStokesDivStable(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Destructor.
    virtual ~IncompressibleNavierStokesDivStable();

    /// Copy constructor.
    IncompressibleNavierStokesDivStable(IncompressibleNavierStokesDivStable const &rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressibleNavierStokesDivStable &operator=(IncompressibleNavierStokesDivStable const &rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    void Initialize(const ProcessInfo &rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType &rResult,
        const ProcessInfo &rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType &rElementalDofList,
        const ProcessInfo &rCurrentProcessInfo) const override;

    void CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    const Parameters GetSpecifications() const override;

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    ///@}
protected:
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
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    /// Pointer to the viscous constitutive model
    ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void SetElementData(
        const ProcessInfo& rProcessInfo,
        ElementDataContainer& rElementData);

    void CalculateKinematics(
        Vector& rGaussWeights,
        Matrix& rVelocityN,
        Matrix& rPressureN,
        GeometryType::ShapeFunctionsGradientsType& rVelocityDNDX,
        GeometryType::ShapeFunctionsGradientsType& rPressureDNDX,
        Vector& rVelocityBubble,
        std::vector<BoundedMatrix<double, 1, TDim>>& rVelocityBubbleGrad,
        DenseVector<GeometryType::ShapeFunctionsSecondDerivativesType>& rVelocityDDNDDX);

    void CalculateStrainRate(ElementDataContainer& rData);

    void ComputeGaussPointLHSContribution(
        const ElementDataContainer& rData,
        MatrixType& rLHS);

    void ComputeGaussPointRHSContribution(
        const ElementDataContainer& rData,
        VectorType& rRHS);

    void ComputeGaussPointEnrichmentContribution(
        const ElementDataContainer& rData,
        array_1d<double, TDim>& rRHSee,
        BoundedMatrix<double, LocalSize, TDim>& rKue,
        BoundedMatrix<double, TDim, LocalSize>& rKeu,
        BoundedMatrix<double, TDim, TDim>& rKee);

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
}; // Class IncompressibleNavierStokesDivStable

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template< unsigned int TDim >
inline std::istream& operator >>(
    std::istream& rIStream,
    IncompressibleNavierStokesDivStable<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const IncompressibleNavierStokesDivStable<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} // Fluid Dynamics Application group

} // namespace Kratos.
