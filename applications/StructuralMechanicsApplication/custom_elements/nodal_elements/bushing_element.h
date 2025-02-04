// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) BushingElement: public Element
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BushingElement);

    ///@}

public:

    ///@name Life Cycle
    ///@{

    /// Default constructors
    BushingElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    BushingElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    ///Copy constructor
    BushingElement(BushingElement const& rOther);

    /// Destructor.
    ~BushingElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BushingElement& operator=(BushingElement const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& ThisNodes) const override;

    //************* GETTING METHODS

    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(
        Vector& rValues,
        int Step = 0) const override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;


    //************* COMPUTING  METHODS

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
     * this calculates just the RHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this calculates just the LHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    static constexpr unsigned int msNumNodes = 2;

    static constexpr unsigned int msLocalSize = 6;

    static constexpr unsigned int msElementSize = msLocalSize * msNumNodes;

    ///@}

    ///@name Protected Operators
    ///@{
    BushingElement() : Element()
    {
    }

    ///@}

private:
    ///@name Private classes
    ///@{

    class Stiffness
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(Stiffness);

        virtual ~Stiffness() = default;

        virtual double GetValue(
            const double Value,
            const Properties& rProperties) const = 0;
    };

    class ConstantStiffness: public Stiffness
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(ConstantStiffness);

        ConstantStiffness(const Variable<double>& rVariable) : mpVariable(&rVariable) {};


        double GetValue(
            const double Value,
            const Properties& rProperties) const override;

    private:
        Variable<double> const * mpVariable;
    };

    class NonLinearStiffness: public Stiffness
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(NonLinearStiffness);

        NonLinearStiffness(
            const Variable<double>& rVariableX,
            const Variable<double>& rVariableY)
            : mpVariableX(&rVariableX),
              mpVariableY(&rVariableY)
        {};


        double GetValue(
            const double Value,
            const Properties& rProperties) const override;

    private:
        Variable<double> const * mpVariableX;

        Variable<double> const * mpVariableY;
    };

    ///@}
    ///@name Private member variables
    ///@{

    array_1d<Stiffness::UniquePointer, 6> mStiffnessGetters;

    ///@}
    ///@name Private Operations
    ///@{

    void GetValuesVector(
        BoundedVector<double, 12>& rValues,
        const int Step) const;

    void CalculateStiffnessValues(
        array_1d<double, 6>& rStiffnessValues,
        const BoundedVector<double, 12>& rValues,
        const double Length) const;

    void CalculateStiffnessMatrix(
        BoundedMatrix<double, 12, 12>& rStiffnessMatrix,
        const array_1d<double, 6>& rStiffnessValues,
        const double Length) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

}; // Class BushingElement

///@}

} // namespace Kratos.
