//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#if !defined(KRATOS_MASS_ELEMENT_H_INCLUDED)
#define KRATOS_MASS_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class MassElement : public Element
{
public:

    ///@name Type Definitions
    ///@{

    typedef Variable< array_1d< double, 3> > ArrayVariableType;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of MassElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MassElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    MassElement(IndexType NewId = 0)
        : Element(NewId) {}

    /**
     * Constructor using an array of nodes
     */
    MassElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes) {}

    /**
     * Constructor using Geometry
     */
    MassElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry) {}

    /**
     * Constructor using Properties
     */
    MassElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties) {}

    /**
     * Copy Constructor
     */
    MassElement(MassElement const& rOther)
        : Element(rOther) {}

    /**
     * Destructor
     */
    ~MassElement() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

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
     * this determines the elemental equation ID vector for all elemental DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;


    /**
     * Getting method to obtain the variable which defines the degrees of freedom
     */
    void GetValuesVector(Vector& values, int Step = 0) const override;

    /**
     * Getting method to obtain the time derivative of variable which defines the degrees of freedom
     */
    void GetFirstDerivativesVector(Vector& values, int Step = 0) const override;

    /**
     * Getting method to obtain the second time derivative of variable which defines the degrees of freedom
     */
    void GetSecondDerivativesVector(Vector& values, int Step = 0) const override;

    /**
     * is called to initialize the element
     * if the element needs to perform any operation before any calculation is done
     * the elemental variables will be initialized and set using this method
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

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
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental damping matrix
     * @param rDampingMatrix: the elemental damping matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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

    ///@name Member Variables
    ///@{

    double mMass = 0.0;

    ///@}

    ///@name Private Operations
    ///@{

    void GenericGetValuesVector(Vector& rValues, int Step, const ArrayVariableType& rVariable) const;

    double GetElementMass() const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    ///@}

}; // Class MassElement

///@}

} // namespace Kratos.

#endif // KRATOS_MASS_ELEMENT_H_INCLUDED defined
