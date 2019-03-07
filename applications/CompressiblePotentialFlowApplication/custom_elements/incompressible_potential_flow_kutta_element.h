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

#if !defined(KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_KUTTA_ELEMENT_H)
#define KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_KUTTA_ELEMENT_H

// Project includes
#include "custom_elements/incompressible_potential_flow_element.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <int Dim, int NumNodes>
class IncompressiblePotentialFlowKuttaElement : public IncompressiblePotentialFlowElement<2,3>
{
public:
    ///@name Pointer Definitions
    /// Pointer definition of IncompressiblePotentialFlowKuttaElement
    KRATOS_CLASS_POINTER_DEFINITION(IncompressiblePotentialFlowKuttaElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    explicit IncompressiblePotentialFlowKuttaElement(IndexType NewId = 0){}

    /**
     * Constructor using an array of nodes
     */
    IncompressiblePotentialFlowKuttaElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : IncompressiblePotentialFlowElement(NewId, ThisNodes){}

    /**
     * Constructor using Geometry
     */
    IncompressiblePotentialFlowKuttaElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : IncompressiblePotentialFlowElement(NewId, pGeometry){}

    /**
     * Constructor using Properties
     */
    IncompressiblePotentialFlowKuttaElement(IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties)
        : IncompressiblePotentialFlowElement(NewId, pGeometry, pProperties){}

    /**
     * Copy Constructor
     */
    IncompressiblePotentialFlowKuttaElement(IncompressiblePotentialFlowKuttaElement const& rOther) = delete;

    /**
     * Move Constructor
     */
    IncompressiblePotentialFlowKuttaElement(IncompressiblePotentialFlowKuttaElement&& rOther) = delete;

    /**
     * Destructor
     */
    ~IncompressiblePotentialFlowKuttaElement() override{}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressiblePotentialFlowKuttaElement& operator=(IncompressiblePotentialFlowKuttaElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator=(rOther);
        return *this;
    }

    /// Move operator.
    IncompressiblePotentialFlowKuttaElement& operator=(IncompressiblePotentialFlowKuttaElement&& rOther)
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

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

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

protected:
    ///@name Protected Operators
    ///@{

    void GetPotential(array_1d<double, NumNodes>& phis) const override;

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

}; // Class IncompressiblePotentialFlowKuttaElement

///@}

} // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_KUTTA_ELEMENT_H  defined
