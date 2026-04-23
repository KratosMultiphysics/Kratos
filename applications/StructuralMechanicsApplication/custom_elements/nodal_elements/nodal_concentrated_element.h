// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"

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

/**
 * @class NodalConcentratedElement
 * @ingroup StructuralMechanicsApplication
 * @brief Concentrated nodal for 3D and 2D points
 * @details The element can consider both the displacement and rotational stiffness, and both the mass and the inertia
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) NodalConcentratedElement
    : public Element
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the geometry
    using GeometryType = Geometry<Node>;

    /// Definition of the base type
    using BaseType = Element;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Counted pointer of NodalConcentratedElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( NodalConcentratedElement);

    ///@}

     /**
     * @brief Flags related to the element computation
     */
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_DISPLACEMENT_STIFFNESS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_NODAL_MASS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_ROTATIONAL_STIFFNESS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_NODAL_INERTIA );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_DAMPING_RATIO );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_ROTATIONAL_DAMPING_RATIO );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RAYLEIGH_DAMPING );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    NodalConcentratedElement(IndexType NewId, GeometryType::Pointer pGeometry);

    NodalConcentratedElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    NodalConcentratedElement(NodalConcentratedElement const& rOther);

    /// Destructor.
    ~NodalConcentratedElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    NodalConcentratedElement& operator=(NodalConcentratedElement const& rOther);

    ///@}
    ///@name Operations
    ///@{

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
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    //************* GETTING METHODS

    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;


    //************* COMPUTING  METHODS

    /**
     * @brief It is called to initialize the element if the element needs to perform any operation before any calculation is done
     * @details The elemental variables will be initialized and set using this method
     * @param rCurrentProcessInfo: the current process info instance
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This is called during the assembling process in order to calculate all elemental contributions to the global system  matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This calculates just the RHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This calculates just the LHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable rDestinationVariable (double version)
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(
        const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        const Variable<double >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable (array_1d<double, 3>) version rDestinationVariable.
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double, 3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    Flags mElementalFlags; /// Elemental flags

    ///@}
    ///@name Protected Operators
    ///@{

    NodalConcentratedElement() : Element()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes the actual size of the system of equations
     * @return This method returns the size of the system of equations
     */
    std::size_t ComputeSizeOfSystem() const;

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
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
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class NodalConcentratedElement

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
