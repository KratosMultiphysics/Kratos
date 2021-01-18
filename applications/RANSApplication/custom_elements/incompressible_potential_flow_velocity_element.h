//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_VELOCITY_ELEMENT_H_INCLUDED)
#define KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_VELOCITY_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_elements/laplace_element.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class IncompressiblePotentialFlowVelocityElement
: public LaplaceElement<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = LaplaceElement<TDim, TNumNodes>;

    using NodeType = Node<3>;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    using PropertiesType = Properties;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = GeometryType::PointsArrayType;

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using IndexType = std::size_t;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of IncompressiblePotentialFlowVelocityElement
    KRATOS_CLASS_POINTER_DEFINITION(IncompressiblePotentialFlowVelocityElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit IncompressiblePotentialFlowVelocityElement(
        IndexType NewId = 0)
    : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    IncompressiblePotentialFlowVelocityElement(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
    : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    IncompressiblePotentialFlowVelocityElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    IncompressiblePotentialFlowVelocityElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    IncompressiblePotentialFlowVelocityElement(
        IncompressiblePotentialFlowVelocityElement const& rOther)
    : BaseType(rOther)
    {
    }

    /**
     * Destructor
     */
    ~IncompressiblePotentialFlowVelocityElement() override = default;

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
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<IncompressiblePotentialFlowVelocityElement>(
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
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<IncompressiblePotentialFlowVelocityElement>(
            NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<IncompressiblePotentialFlowVelocityElement>(
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
        KRATOS_CATCH("");
    }

    const Variable<double>& GetVariable() const override;

    void GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "IncompressiblePotentialFlowVelocityElement #" << this->Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IncompressiblePotentialFlowVelocityElement #" << this->Id();
    }

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }

    ///@}
}; // Class IncompressiblePotentialFlowVelocityElement

///@}
///@name Input and output
///@{

template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(
    std::istream& rIStream,
    IncompressiblePotentialFlowVelocityElement<TDim, TNumNodes>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const IncompressiblePotentialFlowVelocityElement<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_VELOCITY_ELEMENT_H_INCLUDED defined
