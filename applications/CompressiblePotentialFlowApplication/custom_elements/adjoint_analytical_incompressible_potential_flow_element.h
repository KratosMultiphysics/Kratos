//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//
//  Main authors:    Inigo Lopez and Marc Nunez, based on A. Geiser, M. Fusseder and R. Rossi work
//

#if !defined(KRATOS_ADJOINT_ANALYTICAL_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED )
#define KRATOS_ADJOINT_ANALYTICAL_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED


// Project includes
#include "includes/element.h"
#include "adjoint_base_potential_flow_element.h"

namespace Kratos
{

template <class TPrimalElement>
class AdjointAnalyticalIncompressiblePotentialFlowElement : public AdjointBasePotentialFlowElement<TPrimalElement>
{
public:

    ///@name Type Definitions
    ///@{

    typedef AdjointBasePotentialFlowElement<TPrimalElement> BaseType;
    using BaseType::NumNodes;
    using BaseType::Dim;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::NodesArrayType NodesArrayType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of AdjointAnalyticalIncompressiblePotentialFlowElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointAnalyticalIncompressiblePotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit AdjointAnalyticalIncompressiblePotentialFlowElement(IndexType NewId = 0)
     : AdjointBasePotentialFlowElement<TPrimalElement>(NewId){};

    AdjointAnalyticalIncompressiblePotentialFlowElement(IndexType NewId,
                        typename GeometryType::Pointer pGeometry)
     : AdjointBasePotentialFlowElement<TPrimalElement>(NewId, pGeometry){};

    AdjointAnalyticalIncompressiblePotentialFlowElement(IndexType NewId,
                        typename GeometryType::Pointer pGeometry,
                        typename PropertiesType::Pointer pProperties)
     : AdjointBasePotentialFlowElement<TPrimalElement>(NewId, pGeometry, pProperties){};
    /**
     * Copy Constructor
     */
    AdjointAnalyticalIncompressiblePotentialFlowElement(AdjointAnalyticalIncompressiblePotentialFlowElement const& rOther) {};

    /**
     * Destructor
     */
    ~AdjointAnalyticalIncompressiblePotentialFlowElement() override {};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AdjointAnalyticalIncompressiblePotentialFlowElement & operator=(AdjointAnalyticalIncompressiblePotentialFlowElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator =(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const override;

    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;


protected:

private:

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class AdjointAnalyticalIncompressiblePotentialFlowElement


} // namespace Kratos.

#endif
