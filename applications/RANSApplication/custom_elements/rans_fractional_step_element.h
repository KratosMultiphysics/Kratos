//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//  Extended by :    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_FRACTIONAL_STEP_ELEMENT_H_INCLUDED)
#define KRATOS_RANS_FRACTIONAL_STEP_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/fractional_step.h"

// Application includes

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

/// A stabilized element for the incompressible Navier-Stokes equations.
/**
 */
template <unsigned int TDim>
class RansFractionalStepElement : public FractionalStep<TDim>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansFractionalStepElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(RansFractionalStepElement);

    /// Node type (default is: Node<3>)
    using NodeType = Node<3>;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    using IndexType = std::size_t;

    using SizeType = std::size_t;

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    RansFractionalStepElement(
        IndexType NewId = 0)
    : FractionalStep<TDim>(NewId)
    {
    }

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    RansFractionalStepElement(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
    : FractionalStep<TDim>(NewId, ThisNodes)
    {
    }

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    RansFractionalStepElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : FractionalStep<TDim>(NewId, pGeometry)
    {
    }

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    RansFractionalStepElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        Element::PropertiesType::Pointer pProperties)
    : FractionalStep<TDim>(NewId, pGeometry, pProperties)
    {
    }

    /// Destructor.
    ~RansFractionalStepElement() override = default;

    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new RansFractionalStepElement element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Element::PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<RansFractionalStepElement<TDim>>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    /**
     * Returns a pointer to a new RansFractionalStepElement element, created using given input
     * @param NewId the ID of the new element
     * @param pGeom a pointer to the geometry
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */

    Element::Pointer Create(
        IndexType NewId,
        Element::GeometryType::Pointer pGeom,
        Element::PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<RansFractionalStepElement<TDim>>(
            NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */

    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes) const override
    {
        Element::Pointer p_new_element = Create(
            NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

        p_new_element->SetData(this->GetData());
        p_new_element->SetFlags(this->GetFlags());

        return p_new_element;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "RansFractionalStepElement #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RansFractionalStepElement" << TDim << "D";
    }

    ///@}

private:
    ///@name Prvate Access
    ///@{

    void CalculateLocalFractionalVelocitySystem(
        Matrix& rLeftHandSideMatrix,
        Vector& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalPressureSystem(
        Matrix& rLeftHandSideMatrix,
        Vector& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RansFractionalStepElement& operator=(RansFractionalStepElement const& rOther);

    /// Copy constructor.
    RansFractionalStepElement(RansFractionalStepElement const& rOther);

    ///@}

}; // Class RansFractionalStepElement

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim>
inline std::istream& operator>>(
    std::istream& rIStream,
    RansFractionalStepElement<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansFractionalStepElement<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_FRACTIONAL_STEP_ELEMENT_H_INCLUDED  defined
