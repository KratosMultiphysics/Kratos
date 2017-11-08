//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_ROTATION_EMC_SCHEME )
#define  KRATOS_RESIDUAL_BASED_ROTATION_EMC_SCHEME

// System includes

// External includes

// Project includes
#include "custom_strategies/schemes/residual_based_rotation_newmark_scheme.hpp"

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
  
  // Covariant implicit time stepping algorithm: the classical Newmark algorithm of nonlinear elastodynamics and a canonical extension of the Newmark formulas to the orthogonal group SO(3) for the rotational part. (EMC version)

  template<class TSparseSpace, class TDenseSpace >
  class ResidualBasedRotationEMCScheme: public ResidualBasedRotationNewmarkScheme<TSparseSpace,TDenseSpace>
  {
  public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedRotationEMCScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
        
    typedef typename BaseType::TDataType                         TDataType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename Element::DofsVectorType                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                   NodesArrayType;

    typedef ModelPart::ElementsContainerType             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

    typedef typename BaseType::Pointer                     BaseTypePointer;

    typedef BeamMathUtils<double>                        BeamMathUtilsType;

    typedef Quaternion<double>                              QuaternionType;

    typedef ResidualBasedRotationNewmarkScheme<TSparseSpace,TDenseSpace>   DerivedType;
    
  public:
  
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    ResidualBasedRotationEMCScheme(double rDynamic = 1, double rAlpha = 0.0)
      :DerivedType()
    {


      
    }

    
    ///Copy constructor
    ResidualBasedRotationEMCScheme(ResidualBasedRotationEMCScheme& rOther)
      :DerivedType(rOther)
    {
    }


    /// Clone
    BaseTypePointer Clone() override
    {
      return BaseTypePointer( new ResidualBasedRotationEMCScheme(*this) );
    }

    
    /// Destructor
    virtual ~ResidualBasedRotationEMCScheme() override {}


    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /**
     * It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart: The model of the problem to solve
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     *
     */
    
    void InitializeSolutionStep(ModelPart& rModelPart,
				TSystemMatrixType& A,
				TSystemVectorType& Dx,
				TSystemVectorType& b) override
    {
        KRATOS_TRY;

        ProcessInfo CurrentProcessInfo= rModelPart.GetProcessInfo();

        Scheme<TSparseSpace,TDenseSpace>::InitializeSolutionStep(rModelPart, A, Dx, b);

        double DeltaTime = CurrentProcessInfo[DELTA_TIME];

        if (DeltaTime < 1.0e-24)
        {
            KRATOS_ERROR << " ERROR: detected delta_time = 0 in the Solution Scheme DELTA_TIME. PLEASE : check if the time step is created correctly for the current model part ";
        }

	if( this->mDynamic.static_dynamic != 0 ){
	  CurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] = true;  
	  CurrentProcessInfo[EQUILIBRIUM_POINT] = 0.5;
	}
	else{
	  CurrentProcessInfo[EQUILIBRIUM_POINT] = 1;
	}

        //initializing Newmark constants
	this->mDynamic.deltatime = DeltaTime;
	this->mDynamic.c0 = 2.0 / DeltaTime;
	this->mDynamic.c1 = 1.0 / DeltaTime;

        KRATOS_CATCH( "" );
    }
    
    
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


    //*********************************************************************************
    // Predict Linear Movements
    //*********************************************************************************

    virtual void PredictLinearMovements(ModelPart::NodeType& rNode) override
    {
      KRATOS_TRY
	
      //Predicting: NewDisplacement = PreviousDisplacement + PreviousVelocity * DeltaTime;
	  
      array_1d<double, 3 > & CurrentDisplacement      = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double, 3 > & CurrentVelocity          = rNode.FastGetSolutionStepValue(VELOCITY,     0);
      array_1d<double, 3 > & CurrentAcceleration      = rNode.FastGetSolutionStepValue(ACCELERATION, 0);
	  
      array_1d<double, 3 > & CurrentStepDisplacement  = rNode.FastGetSolutionStepValue(STEP_DISPLACEMENT);
	  
      array_1d<double, 3 > & ReferenceDisplacement    = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double, 3 > & ReferenceVelocity        = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double, 3 > & ReferenceAcceleration    = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(CurrentStepDisplacement) = CurrentDisplacement - ReferenceDisplacement;

      //ATTENTION::: the prediction is performed only on free nodes
      this->PredictStepDisplacement( rNode, CurrentStepDisplacement, CurrentVelocity, ReferenceVelocity, CurrentAcceleration, ReferenceAcceleration);	       

      // Updating time derivatives ::: Please note that displacements and its time derivatives can not be consistently fixed separately
      noalias(CurrentDisplacement) = ReferenceDisplacement + CurrentStepDisplacement;

      this->UpdateVelocity     (CurrentVelocity, CurrentStepDisplacement, ReferenceVelocity);
      this->UpdateAcceleration (CurrentAcceleration, CurrentVelocity, ReferenceVelocity);    
      
      KRATOS_CATCH( "" )
    }

    
    //*********************************************************************************
    // Update Linear Movements
    //*********************************************************************************

    virtual void UpdateLinearMovements(ModelPart::NodeType& rNode) override
    {
      KRATOS_TRY
	
      array_1d<double, 3 > DeltaDisplacement;
	  
      // Displacement at iteration i+1
      array_1d<double, 3 > & CurrentDisplacement      = rNode.FastGetSolutionStepValue(DISPLACEMENT,      0);
      // Displacement at iteration i
      array_1d<double, 3 > & ReferenceDisplacement    = rNode.FastGetSolutionStepValue(STEP_DISPLACEMENT, 1); 

      noalias(DeltaDisplacement) = CurrentDisplacement-ReferenceDisplacement;

      array_1d<double, 3 > & CurrentStepDisplacement  = rNode.FastGetSolutionStepValue(STEP_DISPLACEMENT);

      noalias(CurrentStepDisplacement)  = CurrentDisplacement - rNode.FastGetSolutionStepValue(DISPLACEMENT,1);

      array_1d<double, 3 > & CurrentVelocity      = rNode.FastGetSolutionStepValue(VELOCITY,      0);
      array_1d<double, 3 > & CurrentAcceleration  = rNode.FastGetSolutionStepValue(ACCELERATION,  0);
             
      array_1d<double, 3 > & ReferenceVelocity      = rNode.FastGetSolutionStepValue(VELOCITY,     1);

      this->UpdateVelocity     (CurrentVelocity, CurrentStepDisplacement, ReferenceVelocity);
      this->UpdateAcceleration (CurrentAcceleration, CurrentVelocity, ReferenceVelocity);    

      KRATOS_CATCH( "" )
    }


    //*********************************************************************************
    // Predict Angular Movements
    //*********************************************************************************

    virtual void PredictAngularMovements(ModelPart::NodeType& rNode) override
    {
      KRATOS_TRY
	
      array_1d<double, 3 >& CurrentRotation              = rNode.FastGetSolutionStepValue(ROTATION,             0);
      array_1d<double, 3 >& CurrentStepRotation          = rNode.FastGetSolutionStepValue(STEP_ROTATION,        0);
      array_1d<double, 3 >& CurrentAngularVelocity       = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY,     0);
      array_1d<double, 3 >& CurrentAngularAcceleration   = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION, 0);


      array_1d<double, 3 >  ReferenceAngularVelocity     = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY,      1);
      array_1d<double, 3 >  ReferenceAngularAcceleration = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION,  1);

      //ATTENTION::: the prediction is performed only on free nodes
      this->PredictStepRotation( rNode, CurrentStepRotation, CurrentAngularVelocity, ReferenceAngularVelocity, CurrentAngularAcceleration, ReferenceAngularAcceleration);
	    
      // Updating time derivatives ::: Please note that displacements and its time derivatives can not be consistently fixed separately
      noalias(CurrentRotation) += CurrentStepRotation;

      this->UpdateAngularVelocity     (CurrentAngularVelocity, CurrentStepRotation, ReferenceAngularVelocity);	    
      this->UpdateAngularAcceleration (CurrentAngularAcceleration, CurrentStepRotation, CurrentAngularVelocity, ReferenceAngularVelocity);

      KRATOS_CATCH( "" )
    }
    

    //*********************************************************************************
    // Update Angular Movements
    //*********************************************************************************

    virtual void UpdateAngularMovements(ModelPart::NodeType& rNode) override
    {
      KRATOS_TRY
	
      // Rotation at iteration i+1
      array_1d<double, 3 > & CurrentRotation      = rNode.FastGetSolutionStepValue(ROTATION,       0);	  
      // Rotation at iteration i
      array_1d<double, 3 > & PreviousRotation     = rNode.FastGetSolutionStepValue(STEP_ROTATION,  1);
      //Reference Rotation
      array_1d<double, 3 > & ReferenceRotation    = rNode.FastGetSolutionStepValue(ROTATION,       1);	  
      //StepRotation
      array_1d<double, 3 > & CurrentStepRotation  = rNode.FastGetSolutionStepValue(STEP_ROTATION,  0);	
      //DeltaRotation
      array_1d<double, 3 > & CurrentDeltaRotation = rNode.FastGetSolutionStepValue(DELTA_ROTATION, 0);


      // updated delta rotation: (dofs are added, so here the increment is undone)
      noalias(CurrentDeltaRotation) = CurrentRotation-PreviousRotation;      

      if( this->mDynamic.static_dynamic ){

	// updated delta rotation: (dofs are added, so here the increment is undone)
	Matrix CayleyDeltaRotation(3,3);
	noalias(CayleyDeltaRotation) = ZeroMatrix(3,3);
	Matrix CayleyStepRotation(3,3);
	noalias(CayleyStepRotation) = ZeroMatrix(3,3);
	Matrix RotationMatrix(3,3);
	noalias(RotationMatrix) = ZeroMatrix(3,3);

	BeamMathUtilsType::CayleyTransform( CurrentDeltaRotation, CayleyDeltaRotation );
	BeamMathUtilsType::CayleyTransform( CurrentStepRotation, CayleyStepRotation );

	noalias(RotationMatrix) = prod( CayleyDeltaRotation, CayleyStepRotation );

	BeamMathUtilsType::InverseCayleyTransform( RotationMatrix , CurrentStepRotation );


	// updated compound rotation:

	QuaternionType CurrentRotationQuaternion = QuaternionType::FromRotationVector( ReferenceRotation );
	CurrentRotationQuaternion.ToRotationMatrix(RotationMatrix);

	RotationMatrix = prod( CayleyDeltaRotation, RotationMatrix );

	CurrentRotationQuaternion = QuaternionType::FromRotationMatrix( RotationMatrix );
	CurrentRotationQuaternion.ToRotationVector(CurrentRotation);
      }
      else{

	QuaternionType DeltaRotationQuaternion = QuaternionType::FromRotationVector( CurrentDeltaRotation );
	
	QuaternionType StepRotationQuaternion = QuaternionType::FromRotationVector( CurrentStepRotation );
	
	StepRotationQuaternion = DeltaRotationQuaternion * StepRotationQuaternion;
	
	StepRotationQuaternion.ToRotationVector( CurrentStepRotation );
	
	// updated compound rotation:
	QuaternionType RotationQuaternion = QuaternionType::FromRotationVector( PreviousRotation );
	
	RotationQuaternion = DeltaRotationQuaternion * RotationQuaternion;
	
	RotationQuaternion.ToRotationVector( CurrentRotation );
	
      }		
	
      array_1d<double, 3 > & CurrentAngularVelocity        = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY,     0);
      array_1d<double, 3 > & CurrentAngularAcceleration    = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION, 0);
	  
      array_1d<double, 3 > & ReferenceAngularVelocity      = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY,     1);

      this->UpdateAngularVelocity     (CurrentAngularVelocity, CurrentStepRotation, ReferenceAngularVelocity);	    

      this->UpdateAngularAcceleration (CurrentAngularAcceleration, CurrentStepRotation, CurrentAngularVelocity, ReferenceAngularVelocity);

      KRATOS_CATCH( "" )
    }

    
    //*********************************************************************************
    //Predicting Step Displacement variable
    //*********************************************************************************

    virtual void PredictStepDisplacement(ModelPart::NodeType& rNode,
					 array_1d<double, 3 > & CurrentStepDisplacement,
					 const array_1d<double, 3 > & CurrentVelocity,
					 const array_1d<double, 3 > & PreviousVelocity,
					 const array_1d<double, 3 > & CurrentAcceleration,
					 const array_1d<double, 3 > & PreviousAcceleration) override

    {
      KRATOS_TRY

      if (rNode.IsFixed(ACCELERATION_X))
	{
	  CurrentStepDisplacement[0] = ( CurrentAcceleration[0] - PreviousAcceleration[0] ) / this->mDynamic.c1;
	}
      else if (rNode.IsFixed(VELOCITY_X))
	{
	  CurrentStepDisplacement[0] = ( CurrentVelocity[0] - PreviousVelocity[0] ) / this->mDynamic.c0;
	}
      else if (rNode.IsFixed(DISPLACEMENT_X) == false)
	{
	  //CurrentStepDisplacement[0] = this->mDynamic.deltatime * PreviousVelocity[0] + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[0];
	  CurrentStepDisplacement[0] = 0;
	}
      
      if (rNode.IsFixed(ACCELERATION_Y))
	{
	  CurrentStepDisplacement[1] = ( CurrentAcceleration[1] - PreviousAcceleration[1] ) / this->mDynamic.c1;
	}
      else if (rNode.IsFixed(VELOCITY_Y))
	{
	  CurrentStepDisplacement[1] = ( CurrentVelocity[1] - PreviousVelocity[1] ) / this->mDynamic.c0;
	}
      else if (rNode.IsFixed(DISPLACEMENT_Y) == false)
	{
	  //CurrentStepDisplacement[1] = this->mDynamic.deltatime * PreviousVelocity[1] + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[1];
	  CurrentStepDisplacement[1] = 0;
	}
      
      // For 3D cases
      if (rNode.HasDofFor(DISPLACEMENT_Z))
	{
	  if (rNode.IsFixed(ACCELERATION_Z))
	    {
	      CurrentStepDisplacement[2] = ( CurrentAcceleration[2] - PreviousAcceleration[2] ) / this->mDynamic.c1;
	    }
	  else if (rNode.IsFixed(VELOCITY_Z))
	    {
	      CurrentStepDisplacement[2] = ( CurrentVelocity[2] - PreviousVelocity[2] ) / this->mDynamic.c0;
	    }
	  else if (rNode.IsFixed(DISPLACEMENT_Z) == false)
	    {
	      CurrentStepDisplacement[2] = 0;
	    }
      
	}

      KRATOS_CATCH( "" )	      
    }
 

    //*********************************************************************************
    //Predicting step rotation variable
    //*********************************************************************************

    virtual void PredictStepRotation(ModelPart::NodeType& rNode,
				     array_1d<double, 3 > & CurrentStepRotation,
				     const array_1d<double, 3 > & CurrentVelocity,
				     const array_1d<double, 3 > & PreviousVelocity,
				     const array_1d<double, 3 > & CurrentAcceleration,
				     const array_1d<double, 3 > & PreviousAcceleration) override
    {
      KRATOS_TRY
	
      // For 3D cases
      if (rNode.HasDofFor(ROTATION_X))
	{     
	  if (rNode.IsFixed(ANGULAR_ACCELERATION_X))
	    {
	      CurrentStepRotation[0] = ( CurrentAcceleration[0] - PreviousAcceleration[0] ) / this->mDynamic.c1;
	    }
	  else if (rNode.IsFixed(ANGULAR_VELOCITY_X))
	    {
	      CurrentStepRotation[0] = ( CurrentVelocity[0] - PreviousVelocity[0] ) / this->mDynamic.c0;
	    }
	  else if (rNode.IsFixed(ROTATION_X) == false)
	    {
	      //CurrentStepRotation[0] = this->mDynamic.deltatime * PreviousVelocity[0] + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[0];
	      CurrentStepRotation[0] = 0;
	    }
	}

      if (rNode.HasDofFor(ROTATION_Y))
	{

	  if (rNode.IsFixed(ANGULAR_ACCELERATION_Y))
	    {
	      CurrentStepRotation[1] = ( CurrentAcceleration[1] - PreviousAcceleration[1] ) / this->mDynamic.c1;
	    }
	  else if (rNode.IsFixed(ANGULAR_VELOCITY_Y))
	    {
	      CurrentStepRotation[1] = ( CurrentVelocity[1] - PreviousVelocity[1] ) / this->mDynamic.c0;
	    }
	  else if (rNode.IsFixed(ROTATION_Y) == false)
	    {
	      CurrentStepRotation[1] = 0;
	    }
	}
      
      if (rNode.IsFixed(ANGULAR_ACCELERATION_Z))
	{
	  CurrentStepRotation[2] = ( CurrentAcceleration[2] - PreviousAcceleration[2] ) / this->mDynamic.c1;
	}
      else if (rNode.IsFixed(ANGULAR_VELOCITY_Z))
	{
	  CurrentStepRotation[2] = ( CurrentVelocity[2] - PreviousVelocity[2] ) / this->mDynamic.c0;
	}
      else if (rNode.IsFixed(ROTATION_Z) == false)
	{
	  CurrentStepRotation[2] = 0;
	}
      
      KRATOS_CATCH( "" )                             
    }


    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

    inline void UpdateVelocity(array_1d<double, 3 > & CurrentVelocity,
			       const array_1d<double, 3 > & StepDisplacement,
			       const array_1d<double, 3 > & PreviousVelocity)
    { 
      noalias(CurrentVelocity) += ( this->mDynamic.c0 * StepDisplacement - PreviousVelocity ) * this->mDynamic.static_dynamic;  
    }
    
    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************
    
    inline void UpdateAcceleration(array_1d<double, 3 > & CurrentAcceleration,
				   const array_1d<double, 3 > & CurrentVelocity,
				   const array_1d<double, 3 > & PreviousVelocity)
    {   
      noalias(CurrentAcceleration) += ( this->mDynamic.c1 * CurrentVelocity - PreviousVelocity ) * this->mDynamic.static_dynamic;
    }
    

    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

 
    inline void UpdateAngularVelocity(array_1d<double, 3 > & CurrentAngularVelocity,
				      const array_1d<double, 3 > & StepRotation,
                                      const array_1d<double, 3 > & PreviousAngularVelocity)
    {   
      Matrix CayleyStepRotation(3,3);
      noalias(CayleyStepRotation) = ZeroMatrix(3,3);

      BeamMathUtilsType::CayleyTransform( StepRotation, CayleyStepRotation );
      CurrentAngularVelocity = PreviousAngularVelocity; 
      CurrentAngularVelocity = prod(CayleyStepRotation, CurrentAngularVelocity);
      CurrentAngularVelocity *= (-1) * this->mDynamic.static_dynamic;
      
      noalias(CurrentAngularVelocity) += ( this->mDynamic.c0 * StepRotation ) * this->mDynamic.static_dynamic;
    }


    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************

    
    inline void UpdateAngularAcceleration(array_1d<double, 3 > & CurrentAngularAcceleration,
					  const array_1d<double, 3 > & StepRotation,
					  const array_1d<double, 3 > & CurrentAngularVelocity,
                                          const array_1d<double, 3 > & PreviousAngularVelocity)
    {

      Matrix CayleyStepRotation(3,3);
      noalias(CayleyStepRotation) = ZeroMatrix(3,3);

      BeamMathUtilsType::CayleyTransform( StepRotation, CayleyStepRotation );
      CurrentAngularAcceleration = PreviousAngularVelocity; 
      CurrentAngularAcceleration = prod(CayleyStepRotation, CurrentAngularVelocity);
      CurrentAngularAcceleration*= (-1) * this->mDynamic.c1 * this->mDynamic.static_dynamic;

      
      noalias(CurrentAngularAcceleration) += ( this->mDynamic.c1 * CurrentAngularVelocity ) * this->mDynamic.static_dynamic;      
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@{

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
    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
  };  /* Class ResidualBasedRotationEMCScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}  
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ROTATION_EMC_SCHEME  defined */


