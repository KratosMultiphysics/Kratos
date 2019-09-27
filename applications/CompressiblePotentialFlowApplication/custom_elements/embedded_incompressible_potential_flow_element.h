//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez, based on Iñigo Lopez and Riccardo Rossi work
//

#if !defined(KRATOS_EMBEDDED_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H)
#define KRATOS_EMBEDDED_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H

// System includes

// External includes

// Project includes
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "incompressible_potential_flow_element.h"

namespace Kratos
{
template <int Dim, int NumNodes>
class EmbeddedIncompressiblePotentialFlowElement : public IncompressiblePotentialFlowElement<Dim,NumNodes>
{
public:
   ///@name Type Definitions
    ///@{

    typedef IncompressiblePotentialFlowElement<Dim,NumNodes> BaseType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;

    ///@name Pointer Definitions
    /// Pointer definition of EmbeddedIncompressiblePotentialFlowElement
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedIncompressiblePotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    explicit EmbeddedIncompressiblePotentialFlowElement(IndexType NewId = 0){}

    /**
     * Constructor using an array of nodes
     */
    EmbeddedIncompressiblePotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : IncompressiblePotentialFlowElement<Dim,NumNodes>(NewId, ThisNodes){}

    /**
     * Constructor using Geometry
     */
    EmbeddedIncompressiblePotentialFlowElement(IndexType NewId, typename GeometryType::Pointer pGeometry)
        : IncompressiblePotentialFlowElement<Dim,NumNodes>(NewId, pGeometry){}

    /**
     * Constructor using Properties
     */
    EmbeddedIncompressiblePotentialFlowElement(IndexType NewId,
                                       typename GeometryType::Pointer pGeometry,
                                       typename PropertiesType::Pointer pProperties)
        : IncompressiblePotentialFlowElement<Dim,NumNodes>(NewId, pGeometry, pProperties){}

    /**
     * Copy Constructor
     */
    EmbeddedIncompressiblePotentialFlowElement(EmbeddedIncompressiblePotentialFlowElement const& rOther) = delete;

    /**
     * Move Constructor
     */
    EmbeddedIncompressiblePotentialFlowElement(EmbeddedIncompressiblePotentialFlowElement&& rOther) = delete;

    /**
     * Destructor
     */
    ~EmbeddedIncompressiblePotentialFlowElement() override{}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    EmbeddedIncompressiblePotentialFlowElement& operator=(EmbeddedIncompressiblePotentialFlowElement const& rOther) = delete;

    /// Move operator.
    EmbeddedIncompressiblePotentialFlowElement& operator=(EmbeddedIncompressiblePotentialFlowElement&& rOther) = delete;

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            typename PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            typename GeometryType::Pointer pGeom,
                            typename PropertiesType::Pointer pProperties) const override;

    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;
    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;


protected:

    ModifiedShapeFunctions::Pointer pGetModifiedShapeFunctions(Vector& rDistances);

private:

    void CalculateEmbeddedLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo);
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class EmbeddedIncompressiblePotentialFlowElement

} // namespace Kratos.

#endif // KRATOS_EMBEDDED_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H  defined
