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
class IncompressibleNavierStokesP2P1Continuous : public Element
{
public:
    ///@name Type Definitions
    ///@{

    static constexpr std::size_t StrainSize = 3*(TDim-1);

    static constexpr std::size_t VelocityNumNodes = TDim == 2 ? 6 : 10;

    static constexpr std::size_t PressureNumNodes = TDim == 2 ? 3 : 4;

    static constexpr std::size_t LocalSize = VelocityNumNodes*TDim + PressureNumNodes;

    static constexpr IntegrationMethod IntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(IncompressibleNavierStokesP2P1Continuous);

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

    /// Default constructor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    IncompressibleNavierStokesP2P1Continuous(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    IncompressibleNavierStokesP2P1Continuous(
        IndexType NewId,
        const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    IncompressibleNavierStokesP2P1Continuous(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    /// Constructor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    IncompressibleNavierStokesP2P1Continuous(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Destructor.
    virtual ~IncompressibleNavierStokesP2P1Continuous();

    /// Copy constructor.
    IncompressibleNavierStokesP2P1Continuous(IncompressibleNavierStokesP2P1Continuous const &rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressibleNavierStokesP2P1Continuous &operator=(IncompressibleNavierStokesP2P1Continuous const &rOther) = delete;

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

    /**
     * @brief Set the element data container
     * This method accesses the database and fills the provided data container
     * @param rProcessInfo Reference to the current ProcessInfo container
     * @param rElementData Reference to the data container to be filled
     */
    void SetElementData(
        const ProcessInfo& rProcessInfo,
        ElementDataContainer& rElementData);

    /**
     * @brief Calculate the kinematics
     * This method calculates the current element kinematics
     * @param rGaussWeights Integration points weights
     * @param rVelocityN Velocity shape functions at the integration points
     * @param rPressureN Pressure shape functions at the integrationo points
     * @param rVelocityDNDX Velocity shape functions derivatives at the integration points
     * @param rPressureDNDX Pressure shape functions derivatives at the integration points
     * @param rVelocityDDNDDX Velocity shape functions second derivatives at the integration points
     */
    void CalculateKinematics(
        Vector& rGaussWeights,
        Matrix& rVelocityN,
        Matrix& rPressureN,
        GeometryType::ShapeFunctionsGradientsType& rVelocityDNDX,
        GeometryType::ShapeFunctionsGradientsType& rPressureDNDX,
        DenseVector<GeometryType::ShapeFunctionsSecondDerivativesType>& rVelocityDDNDDX);

    /**
     * @brief Calculates the strain rate
     * This method calculates the current element strain rate and stores it in the provided data container
     * @param rData Reference to the data container
     */
    void CalculateStrainRate(ElementDataContainer& rData);

    /**
     * @brief Adds the current Gauss point contribution to LHS
     * This function adds the current Gauss point contribution to the provided Left Hand Side (LHS) matrix
     * @param rData Reference to the data container from which the contribution is computed
     * @param rLHS Reference to the LHS matrix
     */
    void AddGaussPointLeftHandSideContribution(
        const ElementDataContainer& rData,
        MatrixType& rLHS);

    /**
     * @brief Adds the current Gauss point contribution to RHS
     * This function adds the current Gauss point contribution to the provided Right Hand Side (RHS) vector
     * @param rData Reference to the data container from which the contribution is computed
     * @param rRHS Reference to the RHS vector
     */
    void AddGaussPointRightHandSideContribution(
        const ElementDataContainer& rData,
        VectorType& rRHS);

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
}; // Class IncompressibleNavierStokesP2P1Continuous

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
    IncompressibleNavierStokesP2P1Continuous<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const IncompressibleNavierStokesP2P1Continuous<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} // Fluid Dynamics Application group

} // namespace Kratos.
