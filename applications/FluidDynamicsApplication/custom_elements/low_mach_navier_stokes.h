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

template< class TElementData >
class LowMachNavierStokes : public FluidElement<TElementData>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LowMachNavierStokes
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LowMachNavierStokes);

    /// Base type definition
    using BaseType = FluidElement<TElementData>;

    /// Node type (default is: Node)
    using NodeType = typename BaseType::NodeType;

    /// Geometry type (using with given NodeType)
    using GeometryType = typename BaseType::GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = typename BaseType::NodesArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = typename BaseType::VectorType;

    /// Matrix type for local contributions to the linear system
    using MatrixType = typename BaseType::MatrixType;

    using IndexType = typename BaseType::IndexType;

    using SizeType = typename BaseType::SizeType;

    using EquationIdVectorType = typename BaseType::EquationIdVectorType;

    using DofsVectorType = typename BaseType::DofsVectorType;

    using DofsArrayType = typename BaseType::DofsArrayType;

    /// Type for shape function values container
    using ShapeFunctionsType = typename BaseType::ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    using ShapeFunctionDerivativesType = typename BaseType::ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    using ShapeFunctionDerivativesArrayType = typename BaseType::ShapeFunctionDerivativesType;

    constexpr static SizeType Dim = BaseType::Dim;
    constexpr static SizeType NumNodes = BaseType::NumNodes;
    constexpr static SizeType BlockSize = BaseType::BlockSize;
    constexpr static SizeType LocalSize = BaseType::LocalSize;
    constexpr static SizeType StrainSize = (Dim-1)*3;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    LowMachNavierStokes(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    LowMachNavierStokes(
        IndexType NewId,
        const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    LowMachNavierStokes(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    LowMachNavierStokes(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Destructor.
    virtual ~LowMachNavierStokes();

    /// Copy constructor.
    LowMachNavierStokes(LowMachNavierStokes const &rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    LowMachNavierStokes &operator=(LowMachNavierStokes const &rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties) const override;

    void Initialize(const ProcessInfo &rCurrentProcessInfo) override;

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

    void CalculateMaterialResponse(TElementData &rData) const override;

    void AddTimeIntegratedSystem(
        TElementData& rData,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    void AddTimeIntegratedLHS(
        TElementData& rData,
        MatrixType& rLHS) override;

    void AddTimeIntegratedRHS(
        TElementData& rData,
        VectorType& rRHS) override;

    void AddBoundaryTraction(
        TElementData& rData,
        const Vector& rUnitNormal,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    void ComputeGaussPointLHSContribution(
        TElementData& rData,
        MatrixType& rLHS);

    void ComputeGaussPointRHSContribution(
        TElementData& rData,
        VectorType& rRHS);

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
}; // Class LowMachNavierStokes

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template< class TElementData >
inline std::istream& operator >>(
    std::istream& rIStream,
    LowMachNavierStokes<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const LowMachNavierStokes<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} // Fluid Dynamics Application group

} // namespace Kratos.
