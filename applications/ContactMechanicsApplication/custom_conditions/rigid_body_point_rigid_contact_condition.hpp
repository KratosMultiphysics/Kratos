//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_RIGID_BODY_POINT_RIGID_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_RIGID_BODY_POINT_RIGID_CONTACT_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/point_rigid_contact_condition.hpp"

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

/// Rigid Body Point Rigid Contact Condition for 3D and 2D geometries. (base class)

/**
 * Implements a Contact Point Load definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) RigidBodyPointRigidContactCondition
    : public PointRigidContactCondition
{
public:

    ///@name Type Definitions

    ///Tensor order 1 definition
    //typedef BoundedVector<double, 3>     PointType;
    typedef array_1d<double, 3>             PointType;

    ///@{
    // Counted pointer of RigidBodyPointRigidContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( RigidBodyPointRigidContactCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RigidBodyPointRigidContactCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    RigidBodyPointRigidContactCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    RigidBodyPointRigidContactCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall );

    /// Copy constructor
    RigidBodyPointRigidContactCondition( RigidBodyPointRigidContactCondition const& rOther);

    /// Destructor
    virtual ~RigidBodyPointRigidContactCondition();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId,
			      NodesArrayType const& ThisNodes,
			      PropertiesType::Pointer pProperties ) const override;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId,
			     NodesArrayType const& ThisNodes) const override;


    //************* GETTING METHODS

    /**
     * Sets on rConditionDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rConditionDofList,
		    ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult,
			  ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues,
			 int Step = 0 ) override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues,
				   int Step = 0 ) override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues,
				    int Step = 0 ) override;


    /**
     * this function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * @param rRHSVector: input variable containing the RHS vector to be assembled
     * @param rRHSVariable: variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable: variable in the database to which the rRHSvector will be assembled
      * @param rCurrentProcessInfo: the current process info instance
     */
    void AddExplicitContribution(const VectorType& rRHSVector,
					 const Variable<VectorType>& rRHSVariable,
					 Variable<array_1d<double,3> >& rDestinationVariable,
					 const ProcessInfo& rCurrentProcessInfo) override;

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

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
    RigidBodyPointRigidContactCondition() {};


    /**
     * Pointer to the spatial bounding box defining the rigid wall
     */
    SpatialBoundingBox::Pointer mpRigidWall;


    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Initialize System Matrices
     */

    void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					  VectorType& rRightHandSideVector,
					  Flags& rCalculationFlags) override;

    /**
     * Calculate Condition Kinematics
     */
    void CalculateKinematics(ConditionVariables& rVariables,
				     const ProcessInfo& rCurrentProcessInfo,
				     const double& rPointNumber) override;


    /**
     * Calculation of the Load Stiffness Matrix which usually is subtracted to the global stiffness matrix
     */
    void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     ConditionVariables& rVariables,
				     double& rIntegrationWeight) override;

    virtual void CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix,
				     ConditionVariables& rVariables,
				     double& rIntegrationWeight);

    /**
     * Calculation of the External Forces Vector for a force or pressure vector
     */
    void CalculateAndAddContactForces(Vector& rRightHandSideVector,
					      ConditionVariables& rVariables,
					      double& rIntegrationWeight ) override;



    virtual void CalculateAndAddNormalContactForce(Vector& rRightHandSideVector, ConditionVariables& rVariables, double& rIntegrationWeight);


    virtual void CalculateAndAddTangentContactForce(Vector& rRightHandSideVector, ConditionVariables& rVariables, double& rIntegrationWeight);


    double& CalculateNormalForceModulus( double& rNormalForceModulus, ConditionVariables& rVariables );

    double CalculateCoulombsFrictionLaw( double& rTangentForceModulus, double& rNormalForceModulus, ConditionVariables& rVariables );

    double CalculateFrictionCoefficient( const double& rTangentRelativeMovement, const double& rDeltaTime );


    /**
     * Calculation of the Contact Force Factors
     */
    virtual void CalculateContactFactors(ConditionVariables &rContact);


    /**
     * Calculation of an SkewSymmetricTensor from a vector
     */
     void VectorToSkewSymmetricTensor( const Vector& rVector,
				       Matrix& rSkewSymmetricTensor );


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

    WeakPointerVector<Element > mMasterElements;

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }


}; // class RigidBodyPointRigidContactCondition.

} // namespace Kratos.

#endif // KRATOS_RIGID_BODY_POINT_RIGID_CONTACT_CONDITION_H_INCLUDED defined
