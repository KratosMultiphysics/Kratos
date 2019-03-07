//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Riccardo Rossi
//

#if !defined(KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_TRAILING_EDGE_ELEMENT_H)
#define KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_TRAILING_EDGE_ELEMENT_H

// Project includes
#include "custom_elements/incompressible_potential_flow_wake_element.h"
#include "utilities/enrichment_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <int Dim, int NumNodes>
class IncompressiblePotentialFlowTrailingEdgeElement : public IncompressiblePotentialFlowWakeElement<2,3>
{
public:
    ///@name Pointer Definitions
    /// Pointer definition of IncompressiblePotentialFlowTrailingEdgeElement
    KRATOS_CLASS_POINTER_DEFINITION(IncompressiblePotentialFlowTrailingEdgeElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    explicit IncompressiblePotentialFlowTrailingEdgeElement(IndexType NewId = 0){}

    /**
     * Constructor using an array of nodes
     */
    IncompressiblePotentialFlowTrailingEdgeElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : IncompressiblePotentialFlowWakeElement(NewId, ThisNodes){}

    /**
     * Constructor using Geometry
     */
    IncompressiblePotentialFlowTrailingEdgeElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : IncompressiblePotentialFlowWakeElement(NewId, pGeometry){}

    /**
     * Constructor using Properties
     */
    IncompressiblePotentialFlowTrailingEdgeElement(IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties)
        : IncompressiblePotentialFlowWakeElement(NewId, pGeometry, pProperties){}

    /**
     * Copy Constructor
     */
    IncompressiblePotentialFlowTrailingEdgeElement(IncompressiblePotentialFlowTrailingEdgeElement const& rOther) = delete;

    /**
     * Move Constructor
     */
    IncompressiblePotentialFlowTrailingEdgeElement(IncompressiblePotentialFlowTrailingEdgeElement&& rOther) = delete;

    /**
     * Destructor
     */
    ~IncompressiblePotentialFlowTrailingEdgeElement() override{}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressiblePotentialFlowTrailingEdgeElement& operator=(IncompressiblePotentialFlowTrailingEdgeElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator=(rOther);
        return *this;
    }

    /// Move operator.
    IncompressiblePotentialFlowTrailingEdgeElement& operator=(IncompressiblePotentialFlowTrailingEdgeElement&& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{
    void GetValueOnIntegrationPoints(const Variable<bool>& rVariable,
                                     std::vector<bool>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;
    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Private Operators
    ///@{
    void CalculateLocalSystemSubdividedElement(Matrix& lhs_positive, Matrix& lhs_negative);

    void AssignLocalSystemSubdividedElement(MatrixType& rLeftHandSideMatrix,
                                            Matrix& lhs_positive,
                                            Matrix& lhs_negative,
                                            Matrix& lhs_total,
                                            const ElementalData<NumNodes, Dim>& data) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
}; // Class IncompressiblePotentialFlowTrailingEdgeElement

///@}

} // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_TRAILING_EDGE_ELEMENT_H  defined
