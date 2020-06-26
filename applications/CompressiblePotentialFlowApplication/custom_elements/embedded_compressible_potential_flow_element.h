//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez, based on Inigo Lopez and Riccardo Rossi work
//

#if !defined(KRATOS_EMBEDDED_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H)
#define KRATOS_EMBEDDED_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H

// System includes

// External includes

// Project includes
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "compressible_potential_flow_element.h"

namespace Kratos
{
template <int Dim, int NumNodes>
class EmbeddedCompressiblePotentialFlowElement : public CompressiblePotentialFlowElement<Dim,NumNodes>
{
public:
   ///@name Type Definitions
    ///@{

    typedef CompressiblePotentialFlowElement<Dim,NumNodes> BaseType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;

    ///@name Pointer Definitions
    /// Pointer definition of EmbeddedCompressiblePotentialFlowElement
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedCompressiblePotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    explicit EmbeddedCompressiblePotentialFlowElement(IndexType NewId = 0){}

    /**
     * Constructor using an array of nodes
     */
    EmbeddedCompressiblePotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : CompressiblePotentialFlowElement<Dim,NumNodes>(NewId, ThisNodes){}

    /**
     * Constructor using Geometry
     */
    EmbeddedCompressiblePotentialFlowElement(IndexType NewId, typename GeometryType::Pointer pGeometry)
        : CompressiblePotentialFlowElement<Dim,NumNodes>(NewId, pGeometry){}

    /**
     * Constructor using Properties
     */
    EmbeddedCompressiblePotentialFlowElement(IndexType NewId,
                                       typename GeometryType::Pointer pGeometry,
                                       typename PropertiesType::Pointer pProperties)
        : CompressiblePotentialFlowElement<Dim,NumNodes>(NewId, pGeometry, pProperties){}

    /**
     * Copy Constructor
     */
    EmbeddedCompressiblePotentialFlowElement(EmbeddedCompressiblePotentialFlowElement const& rOther) = delete;

    /**
     * Move Constructor
     */
    EmbeddedCompressiblePotentialFlowElement(EmbeddedCompressiblePotentialFlowElement&& rOther) = delete;

    /**
     * Destructor
     */
    ~EmbeddedCompressiblePotentialFlowElement() override{}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    EmbeddedCompressiblePotentialFlowElement& operator=(EmbeddedCompressiblePotentialFlowElement const& rOther) = delete;

    /// Move operator.
    EmbeddedCompressiblePotentialFlowElement& operator=(EmbeddedCompressiblePotentialFlowElement&& rOther) = delete;

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

}; // Class EmbeddedCompressiblePotentialFlowElement

} // namespace Kratos.

#endif // KRATOS_EMBEDDED_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H  defined
