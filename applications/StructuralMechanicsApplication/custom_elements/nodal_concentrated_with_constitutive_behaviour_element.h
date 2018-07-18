// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_NODAL_CONCENTRATED_WITH_CONSTITUTIVE_BEHAVIOUR_ELEMENT_H_INCLUDED )
#define  KRATOS_NODAL_CONCENTRATED_WITH_CONSTITUTIVE_BEHAVIOUR_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/nodal_concentrated_element.h"

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
 * @class NodalConcentratedWithConstitutiveBehaviourElement
 * @ingroup StructuralMechanicsApplication
 * @brief Concentrated nodal for 3D and 2D points with constitutive behaviour
 * @details The element can consider both the displacement and rotational stiffness, and both the mass and the inertia
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) NodalConcentratedWithConstitutiveBehaviourElement
    : public NodalConcentratedElement
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the node type
    typedef Node<3> NodeType;

    /// Definition of the geometry
    typedef Geometry<NodeType> GeometryType;

    /// Definition of the base type
    typedef NodalConcentratedElement BaseType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Counted pointer of NodalConcentratedWithConstitutiveBehaviourElement
    KRATOS_CLASS_POINTER_DEFINITION( NodalConcentratedWithConstitutiveBehaviourElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    NodalConcentratedWithConstitutiveBehaviourElement(IndexType NewId, GeometryType::Pointer pGeometry, ConstitutiveLaw::Pointer pConstitutiveLaw, bool UseRayleighDamping = false, const bool ComputeActiveNodeFlag = true);

    NodalConcentratedWithConstitutiveBehaviourElement(IndexType NewId, GeometryType::Pointer pGeometry, bool UseRayleighDamping = false, const bool ComputeActiveNodeFlag = true);

    NodalConcentratedWithConstitutiveBehaviourElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, bool UseRayleighDamping = false, const bool ComputeActiveNodeFlag = true);

    ///Copy constructor
    NodalConcentratedWithConstitutiveBehaviourElement(NodalConcentratedWithConstitutiveBehaviourElement const& rOther);

    /// Destructor.
    ~NodalConcentratedWithConstitutiveBehaviourElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    NodalConcentratedWithConstitutiveBehaviourElement& operator=(NodalConcentratedWithConstitutiveBehaviourElement const& rOther);

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

    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize() override;

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;


    //************* COMPUTING  METHODS

    /**
     * This is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo
        ) override;


    /**
     * This calculates just the RHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateRightHandSide( 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * This calculates just the LHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateLeftHandSide( 
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix, 
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix, 
        ProcessInfo& rCurrentProcessInfo
        ) override;

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
    ///@}

    ConstitutiveLaw::Pointer mpConstitutiveLaw;

    ///@name Protected Operators
    ///@{
    NodalConcentratedWithConstitutiveBehaviourElement() : BaseType()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{

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

}; // Class NodalConcentratedWithConstitutiveBehaviourElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_NODAL_CONCENTRATED_WITH_CONSTITUTIVE_BEHAVIOUR_ELEMENT_H_INCLUDED  defined
