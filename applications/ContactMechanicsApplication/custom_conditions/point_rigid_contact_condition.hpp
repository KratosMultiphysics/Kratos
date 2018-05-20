//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_POINT_RIGID_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_RIGID_CONTACT_CONDITION_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "custom_bounding/spatial_bounding_box.hpp"
#include "custom_friction/friction_law.hpp"

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

/// Point Rigid Contact Condition for 3D and 2D geometries. (base class)

/**
 * Implements a Contact Point Load definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) PointRigidContactCondition
    : public Condition
{
public:

    ///@name Type Definitions

    ///Tensor order 1 definition
    //typedef BoundedVector<double, 3>     PointType;
    typedef array_1d<double, 3>             PointType;

    ///@{
    // Counted pointer of PointRigidContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( PointRigidContactCondition );
    ///@}

protected:

    /**
     * Flags related to the condition computation
     */

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );


    /**
     * Parameters to be used in the Condition as they are.
     */

    typedef struct
    {
        PointType Normal;        //normal direction
        PointType Tangent;       //tangent direction

    } SurfaceVector;

    typedef struct
    {
        double Normal;        //normal component
        double Tangent;       //tangent component

    } SurfaceScalar;

    struct ConditionVariables
    {
      Flags           Options;               //calculation options

      //Geometrical gaps:
      SurfaceScalar   Gap;                   //normal and tangential gap
      double          ContributoryFactor;    //distance to neighbour points

      //Friction:
      bool            Slip;                  //slip = true / stick = false
      PointType       RelativeDisplacement;  //relative point displacement
      double          FrictionCoefficient;   //total friction coeffitient mu
      double          DeltaTime;             //time step

      SurfaceScalar   Penalty;               //Penalty Parameter normal and tangent

      //Geometric variables
      SurfaceVector   Surface;               //normal and tangent vector to the surface

      // Elasto-plastic constitutive matrix (axisym and PS) (Normal and Tangent stiffness)
      SurfaceScalar   TangentMatrix;

      //Contact stress
      Vector ContactStressVector;

      //for axisymmetric use only
      double  CurrentRadius;
      double  ReferenceRadius;

      //for rigid body
      Matrix  SkewSymDistance;     //compute the skewsymmmetric tensor of the distance

      void Initialize()
      {
	Gap.Normal  = 0;
	Gap.Tangent = 0;

	ContributoryFactor = 0;

	RelativeDisplacement.resize(3);
	noalias(RelativeDisplacement) = ZeroVector(3);

	Slip = false;
	DeltaTime = 0;
	FrictionCoefficient = 0;
	Penalty.Normal  = 0;
	Penalty.Tangent = 0;

	Surface.Normal.resize(3);
	noalias(Surface.Normal) = ZeroVector(3);
	Surface.Tangent.resize(3);
	noalias(Surface.Tangent) = ZeroVector(3);

	SkewSymDistance.resize(3,3);
	noalias(SkewSymDistance) =  ZeroMatrix(3,3);

	TangentMatrix.Normal = 0;
	TangentMatrix.Tangent = 0;

	CurrentRadius = 0;
	ReferenceRadius = 0;
      }

    };


    /**
     * This struct is used in the component wise calculation only
     * is defined here and is used to declare a member variable in the component wise condition
     * private pointers can only be accessed by means of set and get functions
     * this allows to set and not copy the local system variables
     */

    struct LocalSystemComponents
    {
    private:

      //for calculation local system with compacted LHS and RHS
      MatrixType *mpLeftHandSideMatrix;
      VectorType *mpRightHandSideVector;

      //for calculation local system with LHS and RHS components
      std::vector<MatrixType> *mpLeftHandSideMatrices;
      std::vector<VectorType> *mpRightHandSideVectors;

      //LHS variable components
      const std::vector< Variable< MatrixType > > *mpLeftHandSideVariables;

      //RHS variable components
      const std::vector< Variable< VectorType > > *mpRightHandSideVariables;


    public:

      //calculation flags
      Flags  CalculationFlags;

      /**
       * sets the value of a specified pointer variable
       */
      void SetLeftHandSideMatrix( MatrixType& rLeftHandSideMatrix ) { mpLeftHandSideMatrix = &rLeftHandSideMatrix; };
      void SetLeftHandSideMatrices( std::vector<MatrixType>& rLeftHandSideMatrices ) { mpLeftHandSideMatrices = &rLeftHandSideMatrices; };
      void SetLeftHandSideVariables(const std::vector< Variable< MatrixType > >& rLeftHandSideVariables ) { mpLeftHandSideVariables = &rLeftHandSideVariables; };

      void SetRightHandSideVector( VectorType& rRightHandSideVector ) { mpRightHandSideVector = &rRightHandSideVector; };
      void SetRightHandSideVectors( std::vector<VectorType>& rRightHandSideVectors ) { mpRightHandSideVectors = &rRightHandSideVectors; };
      void SetRightHandSideVariables(const std::vector< Variable< VectorType > >& rRightHandSideVariables ) { mpRightHandSideVariables = &rRightHandSideVariables; };


      /**
       * returns the value of a specified pointer variable
       */
      MatrixType& GetLeftHandSideMatrix() { return *mpLeftHandSideMatrix; };
      std::vector<MatrixType>& GetLeftHandSideMatrices() { return *mpLeftHandSideMatrices; };
      const std::vector< Variable< MatrixType > >& GetLeftHandSideVariables() { return *mpLeftHandSideVariables; };

      VectorType& GetRightHandSideVector() { return *mpRightHandSideVector; };
      std::vector<VectorType>& GetRightHandSideVectors() { return *mpRightHandSideVectors; };
      const std::vector< Variable< VectorType > >& GetRightHandSideVariables() { return *mpRightHandSideVariables; };

    };


public:


    ///@name Life Cycle
    ///@{

    /// Serialization constructor
    PointRigidContactCondition(){};

    /// Default constructor
    PointRigidContactCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    PointRigidContactCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );


    PointRigidContactCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall );

    /// Copy constructor
    PointRigidContactCondition( PointRigidContactCondition const& rOther);

    /// Destructor
    virtual ~PointRigidContactCondition();

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


    //************* STARTING - ENDING  METHODS

    /**
     * Called at the beginning of each iteration
     */
    void Initialize() override;

    /**
     * Called at the end of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the beginning of each iteration
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the end of each iteration
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the end of each solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;


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


    //************* COMPUTING  METHODS

    /**
     * this is called during the assembling process in order
     * to calculate all condition contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the condition left hand side matrix
     * @param rRightHandSideVector: the condition right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
			      VectorType& rRightHandSideVector,
			      ProcessInfo& rCurrentProcessInfo ) override;


    /**
     * this function provides a more general interface to the condition.
     * it is designed so that rLHSvariables and rRHSvariables are passed TO the condition
     * thus telling what is the desired output
     * @param rLeftHandSideMatrices: container with the output left hand side matrices
     * @param rLHSVariables: paramter describing the expected LHSs
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateLocalSystem(std::vector< MatrixType >& rLeftHandSideMatrices,
			      const std::vector< Variable< MatrixType > >& rLHSVariables,
			      std::vector< VectorType >& rRightHandSideVectors,
			      const std::vector< Variable< VectorType > >& rRHSVariables,
			      ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the condition right hand side vector only
      * @param rRightHandSideVector: the condition right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo ) override;


    /**
     * this function provides a more general interface to the condition.
     * it is designed so that rRHSvariables are passed TO the condition
     * thus telling what is the desired output
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateRightHandSide(std::vector< VectorType >& rRightHandSideVectors,
				const std::vector< Variable< VectorType > >& rRHSVariables,
				ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the condition mass matrix
      * @param rMassMatrix: the condition mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo ) override;

    /**
      * this is called during the assembling process in order
      * to calculate the condition damping matrix
      * @param rDampingMatrix: the condition damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo ) override;


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

    /**
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;


    /**
     * Pointer to the spatial bounding box defining the rigid wall
     */
    SpatialBoundingBox::Pointer mpRigidWall;


    /**
     * Pointer to the friction law
     */
    FrictionLaw::Pointer mpFrictionLaw;


    /**
     * Contact stress tensor
     */
    Vector mContactStressVector;

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


    /**
     * Initialize System Matrices
     */
    virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					  VectorType& rRightHandSideVector,
					  Flags& rCalculationFlags);

    /**
     * Initialize General Variables
     */
    virtual void InitializeConditionVariables(ConditionVariables& rVariables,
					    const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(ConditionVariables& rVariables,
				     const ProcessInfo& rCurrentProcessInfo,
				     const double& rPointNumber);

    /**
     * Calculation of the Integration Weight
     */
    virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);


    /**
     * Calculates the condition contributions
     */
    virtual void CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
					  ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculation and addition of the matrices of the LHS
     */
    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    ConditionVariables& rVariables,
                                    double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */
    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    ConditionVariables& rVariables,
				    double& rIntegrationWeight);


    /**
     * Calculation of the Load Stiffness Matrix which usually is subtracted to the global stiffness matrix
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     ConditionVariables& rVariables,
				     double& rIntegrationWeight);


    /**
     * Calculation of the Contact Forces Vector
     */
    virtual void CalculateAndAddContactForces(Vector& rRightHandSideVector,
					      ConditionVariables& rVariables,
					      double& rIntegrationWeight );

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
	  //rSerializer.save("mpRigidWall",mpRigidWall); //rebuild in contact search and in restart
	rSerializer.save("mpFrictionLaw",mpFrictionLaw);
	rSerializer.save("mContactStressVector",mContactStressVector);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
	  //rSerializer.load("mpRigidWall",mpRigidWall);
	rSerializer.load("mpFrictionLaw",mpFrictionLaw);
	rSerializer.load("mContactStressVector",mContactStressVector);
    }

    ///@}

}; // class PointRigidContactCondition.

} // namespace Kratos.

#endif // KRATOS_POINT_RIGID_CONTACT_CONDITION_H_INCLUDED defined
