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

#ifndef KRATOS_SYMBOLIC_STOKES_H
#define KRATOS_SYMBOLIC_STOKES_H

// External includes

// System includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"

// Application includes
#include "includes/cfd_variables.h"
#include "custom_elements/fluid_element.h"

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
class SymbolicStokes : public FluidElement<TElementData>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SymbolicStokes
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SymbolicStokes);

    /// Node type (default is: Node)
    typedef Node NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    /// Type for shape function values container
    typedef typename FluidElement<TElementData>::ShapeFunctionsType ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef typename FluidElement<TElementData>::ShapeFunctionDerivativesType ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef typename FluidElement<TElementData>::ShapeFunctionDerivativesArrayType ShapeFunctionDerivativesArrayType;

    constexpr static unsigned int Dim = FluidElement<TElementData>::Dim;
    constexpr static unsigned int NumNodes = FluidElement<TElementData>::NumNodes;
    constexpr static unsigned int BlockSize = FluidElement<TElementData>::BlockSize;
    constexpr static unsigned int LocalSize = FluidElement<TElementData>::LocalSize;

    constexpr static unsigned int StrainSize = (Dim-1)*3;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constructor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    SymbolicStokes(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    SymbolicStokes(
        IndexType NewId,
        const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    SymbolicStokes(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    /// Constructor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    SymbolicStokes(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Destructor.
    virtual ~SymbolicStokes();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new SymbolicStokes element, created using given input.
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    /// Create a new element of this type using given geometry
    /**
     * Returns a pointer to a new FluidElement element, created using given input.
     * @param NewId the ID of the new element
     * @param pGeom a pointer to the geometry to be used to create the element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

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

    void Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

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

    /// Assignment operator.
    SymbolicStokes& operator=(SymbolicStokes const& rOther);

    /// Copy constructor.
    SymbolicStokes(SymbolicStokes const& rOther);

    ///@}
}; // Class SymbolicStokes

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
    SymbolicStokes<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const SymbolicStokes<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group
} // namespace Kratos.

#endif // KRATOS_SYMBOLIC_STOKES_H
