//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

#if !defined(KRATOS_SURFACE_SMOOTHING_ELEMENT_H_INCLUDED )
#define KRATOS_SURFACE_SMOOTHING_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/element.h"

// Include Base headers
#include "includes/checks.h"

namespace Kratos
{

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

class SurfaceSmoothingElement : public Element
{
public:

    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of CellNodeProjectionElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SurfaceSmoothingElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    SurfaceSmoothingElement(IndexType NewId = 0);

    /**
     * Constructor using an array of nodes
     */
    SurfaceSmoothingElement(IndexType NewId, const NodesArrayType& ThisNodes);

    /**
     * Constructor using Geometry
     */
    SurfaceSmoothingElement(IndexType NewId, GeometryType::Pointer pGeometry);

    /**
     * Constructor using Properties
     */
    SurfaceSmoothingElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
     * Copy Constructor
     */
    SurfaceSmoothingElement(SurfaceSmoothingElement const& rOther);

    /**
     * Destructor
     */
    ~SurfaceSmoothingElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SurfaceSmoothingElement & operator=(SurfaceSmoothingElement const& rOther);

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
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
     * they can be managed internally with a private method to do the same calculations
     * only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

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

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

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

}; // Class SurfaceSmoothingElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_SURFACE_SMOOTHING_ELEMENT_H_INCLUDED