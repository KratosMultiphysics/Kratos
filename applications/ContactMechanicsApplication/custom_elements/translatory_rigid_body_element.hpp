//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_TRANSLATORY_RIGID_BODY_ELEMENT_H_INCLUDED )
#define  KRATOS_TRANSLATORY_RIGID_BODY_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/rigid_body_element.hpp"

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

/// Rigid Body Element for 3D space dimension

/**
 * Nodal Variables: DISPLACEMENT, STEP_DISPLACEMENT, VELOCITY, ACCELERATION, ROTATION, STEP_ROTATION, DELTA_ROTATION, ANGULAR_VELOCITY, ANGULAR_ACCELERATION
 * Nodal Dofs: DISPLACEMENT, ROTATION
 */

class KRATOS_API(CONTACT_MECHANICS_APPLICATION) TranslatoryRigidBodyElement
    :public RigidBodyElement
{
public:

    ///@name Type Definitions
    ///@{
   ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw                          ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer      ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure         StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod            IntegrationMethod;
    ///Type definition for beam utilities
    typedef BeamMathUtils<double>                      BeamMathUtilsType;
    ///Type definition for quaternion
    typedef Quaternion<double>                            QuaternionType;
    ///Type for nodes
    typedef Node<3>                                             NodeType;
    ///Type for nodes container
    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;


    /// Counted pointer of TranslatoryRigidBodyElement
    KRATOS_CLASS_POINTER_DEFINITION( TranslatoryRigidBodyElement );

    ///@}

public:

    ///@name Life Cycle
    ///@{

    /// Serializer constructor
    TranslatoryRigidBodyElement() {};

    /// Default constructors
    TranslatoryRigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry);

    TranslatoryRigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    TranslatoryRigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, NodesContainerType::Pointer pNodes);

    ///Copy constructor
    TranslatoryRigidBodyElement(TranslatoryRigidBodyElement const& rOther);

    /// Destructor.
    virtual ~TranslatoryRigidBodyElement();


    ///@}
    ///@name Operators
    ///@{


   /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;


    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone (IndexType NewId, NodesArrayType const& ThisNodes) const override;

    //************* GETTING METHODS


    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues, int Step = 0) override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;


    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize() override;


    //************* COMPUTING  METHODS

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Translatory Rigid Body Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Translatory Rigid Body Element #" << Id();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      GetGeometry().PrintData(rOStream);
    }
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
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Initialize System Matrices
     */
    virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          Flags& rCalculationFlags) override;



    /**
      * Calculation of the Tangent Intertia Matrix
      */
    virtual void CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix,
					   ElementVariables& rVariables,
					   ProcessInfo& rCurrentProcessInfo) override;

    /**
      * Calculation of the Inertial Forces Vector
      */
    virtual void CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector,
					   ElementVariables& rVariables,
					   ProcessInfo& rCurrentProcessInfo) override;


    /**
      * Update rigid body nodes and positions
      */
    virtual void UpdateRigidBodyNodes(ProcessInfo& rCurrentProcessInfo) override;



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
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override;

    virtual void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}


}; // Class TranslatoryRigidBodyElement

} // namespace Kratos.
#endif //  KRATOS_TRANSLATORY_RIGID_BODY_ELEMENT_H_INCLUDED defined

