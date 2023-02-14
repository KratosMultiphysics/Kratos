//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//                   Pablo Becker
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "geometries/geometry_data.h"

namespace Kratos
{

/*The "DistanceCalculationFluxBasedElement" element formulation is meant to compute the pseudo filltime or flowlength for a filling problem
Currently only supports very viscous problems (polymers or RTM)
It requires two steps to be performed:
1. Solves a transient diffusion problem to get a field equivalent to a potential.
2. from the field, compute gradients to obtain velocity and use it to compute flowlenght/filltime.
*/

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

template< unsigned int TDim,  unsigned int TNumNodes >
class KRATOS_API(KRATOS_CORE) DistanceCalculationFluxBasedElement : public Element
{
public:

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(DistanceCalculationFluxBasedElement);

    ///@name Type Definitions
    ///@{

    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    static constexpr unsigned int BlockSize = 1;

    static constexpr unsigned int LocalSize = TNumNodes * BlockSize;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constuctor.
    /**
    * @param NewId Index number of the new element (optional)
    */
    DistanceCalculationFluxBasedElement(IndexType NewId = 0) :
        Element(NewId)
    {}

    /// Constructor using an array of nodes.
    /**
    * @param NewId Index of the new element
    * @param ThisNodes An array containing the nodes of the new element
    */
    DistanceCalculationFluxBasedElement(
        IndexType NewId, 
        const NodesArrayType& ThisNodes):
        Element(NewId, ThisNodes)
    {}

    /// Constructor using a geometry object.
    /**
    * @param NewId Index of the new element
    * @param pGeometry Pointer to a geometry object
    */
    DistanceCalculationFluxBasedElement(
        IndexType NewId,
         GeometryType::Pointer pGeometry):
        Element(NewId, pGeometry)
    {}

    /// Constuctor using geometry and properties.
    /**
    * @param NewId Index of the new element
    * @param pGeometry Pointer to a geometry object
    * @param pProperties Pointer to the element's properties
    */
    DistanceCalculationFluxBasedElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~DistanceCalculationFluxBasedElement() override = default;


    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
    * Returns a pointer to a new DistanceCalculationFluxBasedElement element, created using given input.
    * @param NewId the ID of the new element
    * @param ThisNodes the nodes of the new element
    * @param pProperties the properties assigned to the new element
    * @return a Pointer to the new element
    */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties)const override
    {
        return Kratos::make_intrusive<DistanceCalculationFluxBasedElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new element of this type using given geometry
    /**
    * Returns a pointer to a new DistanceCalculationFluxBasedElement element, created using given input.
    * @param NewId the ID of the new element
    * @param pGeom a pointer to the geomerty to be used to create the element
    * @param pProperties the properties assigned to the new element
    * @return a Pointer to the new element
    */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override 
    {
        return Kratos::make_intrusive<DistanceCalculationFluxBasedElement>(NewId, pGeom, pProperties);
    }

    /// Computes the elemental LHS and RHS elemental contributions
    /**
     * Computes the Left Hand Side (LHS)
     * and Right Hand Side elemental contributions for the element.
     * @param rLeftHandSideMatrix elemental stiffness matrix
     * @param rRightHandSideVector elemental residual vector
     * @param rCurrentProcessInfo reference to the current process info
     */
    void CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;

    /// Computes the elemental contribution to assemble the velocity
    /**
     * @param rCurrentProcessInfo reference to the current process info
     */
    void AddExplicitContribution(
        const ProcessInfo &rCurrentProcessInfo) override;


    int Check(
        const ProcessInfo &rCurrentProcessInfo) const override;

    /**
     * @brief EquationIdVector Returns the global system rows corresponding to each local row.
     * @param rResult rResult[i] is the global index of local row i (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief GetDofList Returns a list of the element's Dofs.
     * @param rElementalDofList List of DOFs. (output)
     * @param rCurrentProcessInfo Current ProcessInfo instance. (input)
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;


    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;
    
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
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


    /**
     * @brief Assembles the local system to get a gradient in the direction of the flow
     * @return A transient diffusion system whose solution can be used to construct a potential field
     * this potential will then be used to compute the distance in a second step
     */
    void CalculatePotentialFlowSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief Assembles the local system to get a the distance field
     * @return A pure convection + sorce term system that solves the distance.
     */
    void CalculateDistanceSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief Computes gauss points
     * @return gauss weights and positions
     */
    void CalculateGaussPointsData(
        const GeometryType& rGeometry,
        BoundedVector<double,TNumNodes> &rGaussWeights,
        BoundedMatrix<double,TNumNodes,TNumNodes> &rNContainer,
        array_1d<BoundedMatrix<double, TNumNodes, TDim>, TNumNodes>& rDN_DXContainer);
    
    /**
     * @brief Computes gauss points for simplex geometries
     * @return gauss points coords
     */
    void GetSimplexShapeFunctionsOnGauss(BoundedMatrix<double,TNumNodes, TNumNodes>& Ncontainer);


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
    ///@name Serialization
    ///@{

private:

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    DistanceCalculationFluxBasedElement & operator=(DistanceCalculationFluxBasedElement const& rOther) = delete;

    /// Copy constructor.
    DistanceCalculationFluxBasedElement(DistanceCalculationFluxBasedElement const& rOther) = delete;

    ///@}

}; // Class DistanceCalculationFluxBasedElement
///@name Type Definitions
///@{



} // namespace Kratos.
