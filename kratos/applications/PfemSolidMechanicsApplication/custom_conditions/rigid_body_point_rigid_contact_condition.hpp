//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:           September 2014 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RIGID_BODY_POINT_RIGID_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_RIGID_BODY_POINT_RIGID_CONTACT_CONDITION_H_INCLUDED



// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "custom_modelers/spatial_bounding_box.hpp"

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
class RigidBodyPointRigidContactCondition
    : public Condition
{
public:

    ///@name Type Definitions

    ///Tensor order 1 definition
    typedef Vector VectorType;

    ///@{
    // Counted pointer of RigidBodyPointRigidContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( RigidBodyPointRigidContactCondition );
    ///@}

protected:

    /**
     * Parameters to be used in the Condition as they are. 
     */

   typedef struct
    {
        VectorType Normal;        //normal direction
        VectorType Tangent;       //tangent direction
	   
    } SurfaceVector;

    typedef struct
    {
        double Normal;        //normal component 
        double Tangent;       //tangent component
        	   
    } SurfaceScalar;

    typedef struct
    {
      Flags           Options;               //calculation options
      
      //Geometrical gaps:
      SurfaceScalar   Gap;                   //normal and tangential gap

      //Friction:
      double          FrictionCoefficient;   //total friction coeffitient mu

      SurfaceScalar   Penalty;               //Penalty Parameter normal and tangent

      //Geometric variables
      SurfaceVector   Surface;               //normal and tangent vector to the surface
      

       //for axisymmetric use only
      double  CurrentRadius;
      double  ReferenceRadius;


    } GeneralVariables;


public:


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
			      PropertiesType::Pointer pProperties ) const;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, 
			     NodesArrayType const& ThisNodes) const;


    //************* STARTING - ENDING  METHODS


    /**
     * Called at the beginning of each iteration
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);


    /**
     * Called at the end of each iteration
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);


    /**
     * Called at the beginning of each solution step
     */
    virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);


    /**
     * Called at the end of each solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);


    //************* GETTING METHODS

    /**
     * Sets on rConditionDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rConditionDofList,
		    ProcessInfo& rCurrentProcessInfo );

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult,
			  ProcessInfo& rCurrentProcessInfo );

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues,
			 int Step = 0 );

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues,
				   int Step = 0 );

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues,
				    int Step = 0 );


    //************* COMPUTING  METHODS

    /**
     * this function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * @param rRHSVector: input variable containing the RHS vector to be assembled
     * @param rRHSVariable: variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable: variable in the database to which the rRHSvector will be assembled 
      * @param rCurrentProcessInfo: the current process info instance
     */      
    virtual void AddExplicitContribution(const VectorType& rRHSVector, 
					 const Variable<VectorType>& rRHSVariable, 
					 Variable<array_1d<double,3> >& rDestinationVariable, 
					 const ProcessInfo& rCurrentProcessInfo);

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );


    /**
     * Rigid Contact Condition has a spatial bounding box defining the rigid wall properites
     * it must be set and it is a member pointer variable
     */
    void SetRigidWall(SpatialBoundingBox::Pointer pRigidWall);


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
     * Clear Nodal Forces
     */
    void ClearNodalForces ();


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

    WeakPointerVector< Element > mMasterElements;

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

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }


}; // class RigidBodyPointRigidContactCondition.

} // namespace Kratos.

#endif // KRATOS_RIGID_BODY_POINT_RIGID_CONTACT_CONDITION_H_INCLUDED defined 
