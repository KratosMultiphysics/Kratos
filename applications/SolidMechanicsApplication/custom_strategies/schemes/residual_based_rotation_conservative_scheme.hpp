//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2015 $
//   Revision:            $Revision:                  0.0 $
//
// 

#if !defined(KRATOS_RESIDUAL_BASED_ROTATION_CONSERVATIVE_SCHEME )
#define  KRATOS_RESIDUAL_BASED_ROTATION_CONSERVATIVE_SCHEME

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "containers/array_1d.h"

#include "utilities/beam_math_utilities.hpp"


#include "solid_mechanics_application_variables.h"

namespace Kratos
{
/**@name Kratos Globals */
/*@{ */
/*@} */
/**@name Type Definitions */
/*@{ */
/*@} */
/**@name  Enum's */
/*@{ */
/*@} */
/**@name  Functions */
/*@{ */
/*@} */
/**@name Kratos Classes */
/*@{ */
/*@} */

// Energy-momentum conserving scheme in dynamics

template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedRotationConservativeScheme: public Scheme<TSparseSpace,TDenseSpace>
{
protected:


    struct GeneralDynamics
    {
      double c0;
      double c1;
      double deltatime;
      double static_dynamic;

    };


    struct  GeneralMatrices
    {

        std::vector< Matrix > M;     //first derivative matrix  (usually mass matrix)
        std::vector< Matrix > D;     //second derivative matrix (usually damping matrix)

    };

    struct GeneralVectors
    {
        std::vector< Vector > fm;  //fist derivatives vector
        std::vector< Vector > fd;  //second derivative vector
    };



public:


    /**@name Type Definitions */

    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedRotationConservativeScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

    typedef typename BaseType::TDataType                         TDataType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename Element::DofsVectorType                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

    typedef BeamMathUtils<double>                        BeamMathUtilsType;

    typedef Quaternion<double>                              QuaternionType;

    /*@} */

    /**
     * Constructor.
     * The bossak method
     */
    ResidualBasedRotationConservativeScheme(double rDynamic = 1)
        :Scheme<TSparseSpace,TDenseSpace>()
    {
        //For pure stable Newmark Scheme
        mDynamic.static_dynamic= rDynamic;


        //std::cout << " MECHANICAL SCHEME: The Romero and Armero Conservative Scheme "<<std::endl;


        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();

        mMatrix.M.resize(NumThreads);
        mMatrix.D.resize(NumThreads);

        mVector.fm.resize(NumThreads);
	mVector.fd.resize(NumThreads);

    }

    /** Destructor.
     */
    virtual ~ResidualBasedRotationConservativeScheme
    () {}

    /**@name Operators */

    /*@} */
    /**@name Operators
     */
    /*@{ */

    /**
       Performing the update of the solution.
    */


    //***************************************************************************
    /**
     * incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param r_model_part
     * @param rDofSet set of all primary variables
     * @param A	LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b )
    {
        KRATOS_TRY


        //store previous iteration rotation in previous step rotation
	for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                i != r_model_part.NodesEnd(); ++i)
	  {
	    if( i->IsNot(SLAVE) && i->IsNot(RIGID) ){


	      //DELTA_ROTATION variable used to store previous iteration rotation
	      array_1d<double, 3 > & CurrentRotation      = (i)->FastGetSolutionStepValue(ROTATION);
	      array_1d<double, 3 > & ReferenceRotation    = (i)->FastGetSolutionStepValue(DELTA_ROTATION, 1);
    
	      // Rotation at iteration i
	      ReferenceRotation    = CurrentRotation; 

	    }
	  }


        //std::cout << " Update " << std::endl;
        //update of displacement (by DOF)
        for (typename DofsArrayType::iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
        {
            if (i_dof->IsFree() )
            {
                i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
            }
        }

	//LINEAR VELOCITIES AND ACCELERATIONS

        //updating time derivatives (nodally for efficiency)
	for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                i != r_model_part.NodesEnd(); ++i)
	  {

	    if( i->IsNot(SLAVE) && i->IsNot(RIGID) ){
	    
	      array_1d<double, 3 > & CurrentStepDisplacement  = (i)->FastGetSolutionStepValue(STEP_DISPLACEMENT);

	      CurrentStepDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT) - (i)->FastGetSolutionStepValue(DISPLACEMENT,1);

	      array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY);
	      array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION);
             
              array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY,1);
	      array_1d<double, 3 > & PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);

              UpdateVelocity     ((*i), CurrentVelocity, PreviousVelocity, CurrentStepDisplacement);
              UpdateAcceleration ((*i), CurrentAcceleration, PreviousAcceleration, CurrentVelocity, PreviousVelocity);
                            
	    }
	  }

        //updating time derivatives (nodally for efficiency)

	//UPDATE COMPOUND AND DELTA ROTATIONS BY QUATERNION COMPOSITION and ANGULAR VELOCITIES AND ACCELERATIONS
	Vector CurrentRotationVector      = ZeroVector(3);
	Vector PreviousRotationVector     = ZeroVector(3);
	Vector CurrentStepRotationVector  = ZeroVector(3);
	Vector DeltaRotationVector        = ZeroVector(3);


	QuaternionType DeltaRotationQuaternion;
	QuaternionType CurrentRotationQuaternion;
	QuaternionType ReferenceRotationQuaternion;


	for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
	     i != r_model_part.NodesEnd(); ++i)
	  {
	    if( i->IsNot(SLAVE) && i->IsNot(RIGID) ){

	      // Rotation at iteration i+1
	      array_1d<double, 3 > & CurrentRotation     = (i)->FastGetSolutionStepValue(ROTATION);
	      // Rotation at step n
	      array_1d<double, 3 > & PreviousRotation    = (i)->FastGetSolutionStepValue(ROTATION, 1);
	      // Rotation at iteration i
	      array_1d<double, 3 > & ReferenceRotation   = (i)->FastGetSolutionStepValue(DELTA_ROTATION, 1);
	      // StepRotation
	      array_1d<double, 3 > & CurrentStepRotation = (i)->FastGetSolutionStepValue(STEP_ROTATION);
	      // DeltaRotation
	      array_1d<double, 3 > & CurrentDeltaRotation= (i)->FastGetSolutionStepValue(DELTA_ROTATION);

	      // std::cout<<(i)->Id()<<" CurrentRotation "<<CurrentRotation<<std::endl;
	      // std::cout<<(i)->Id()<<" ReferenceRotation "<<ReferenceRotation<<std::endl;
	      // std::cout<<(i)->Id()<<" CurrentStepRotation "<<CurrentStepRotation<<std::endl;

	      CurrentDeltaRotation = CurrentRotation-ReferenceRotation;
	    
	      // std::cout<<" CurrentDeltaRotation "<<CurrentDeltaRotation<<std::endl;

	      for( unsigned int j=0; j<3; j++)
		{
		  DeltaRotationVector[j]       = CurrentDeltaRotation[j];  //delta rotation computed in the iteration

		  CurrentRotationVector[j]     = ReferenceRotation[j];     //previous iteration total rotation

		  PreviousRotationVector[j]    = PreviousRotation[j];      //previous step total rotation

		  CurrentStepRotationVector[j] = CurrentStepRotation[j];   //previous iteration step rotation 
		}
	   

	      if( mDynamic.static_dynamic ){
		
		// Option I:

		//A:
		// Matrix CayleyDeltaRotation = ZeroMatrix(3,3);
		// Matrix CurrentRotationMatrix = ZeroMatrix(3,3);
		// Matrix PreviousRotationMatrix = ZeroMatrix(3,3);

		// BeamMathUtilsType::CayleyTransform( DeltaRotationVector, CayleyDeltaRotation );
		
		// CurrentRotationQuaternion = QuaternionType::FromRotationVector( CurrentRotationVector );
		// CurrentRotationQuaternion.ToRotationMatrix(CurrentRotationMatrix);

		// CurrentRotationMatrix = prod( CayleyDeltaRotation, CurrentRotationMatrix );

		// CurrentRotationQuaternion = QuaternionType::FromRotationMatrix( CurrentRotationMatrix );
		// CurrentRotationQuaternion.ToRotationVector(CurrentRotationVector); //assign

		// //----------

		// Matrix CayleyStepRotation = ZeroMatrix(3,3);
		// ReferenceRotationQuaternion = QuaternionType::FromRotationVector( PreviousRotationVector );
		// ReferenceRotationQuaternion.ToRotationMatrix(PreviousRotationMatrix);

		// CayleyStepRotation =  prod( CurrentRotationMatrix, trans(PreviousRotationMatrix) ); 

		// BeamMathUtilsType::InverseCayleyTransform( CayleyStepRotation, CurrentStepRotationVector ); //assign


		// Option II: 
		Matrix CayleyDeltaRotation   = ZeroMatrix(3,3);
		Matrix CayleyStepRotation    = ZeroMatrix(3,3);

		Matrix StepRotationMatrix    = ZeroMatrix(3,3);

		BeamMathUtilsType::CayleyTransform( DeltaRotationVector, CayleyDeltaRotation );
		BeamMathUtilsType::CayleyTransform( CurrentStepRotationVector,  CayleyStepRotation );
		// std::cout<<" CayleyDeltaRotation "<<CayleyDeltaRotation<<std::endl;
		// std::cout<<" CayleyStepRotation "<<CayleyStepRotation<<std::endl;

		StepRotationMatrix = prod( CayleyDeltaRotation, CayleyStepRotation );

		BeamMathUtilsType::InverseCayleyTransform( StepRotationMatrix , CurrentStepRotationVector ); //assign

		//------

		// A:
		Matrix CurrentRotationMatrix = ZeroMatrix(3,3);

		CurrentRotationQuaternion = QuaternionType::FromRotationVector( CurrentRotationVector );
		CurrentRotationQuaternion.ToRotationMatrix(CurrentRotationMatrix);

		CurrentRotationMatrix = prod( CayleyDeltaRotation, CurrentRotationMatrix );

		CurrentRotationQuaternion = QuaternionType::FromRotationMatrix( CurrentRotationMatrix );
		CurrentRotationQuaternion.ToRotationVector(CurrentRotationVector); //assign

		// std::cout<<(i)->Id()<<" CurrentStepRotationVector "<<CurrentStepRotationVector<<std::endl;
		
		//B:
		// Matrix CurrentRotationMatrix  = ZeroMatrix(3,3);
		// Matrix PreviousRotationMatrix = ZeroMatrix(3,3);

		// CayleyStepRotation = ZeroMatrix(3,3);
		// BeamMathUtilsType::CayleyTransform( CurrentStepRotationVector,  CayleyStepRotation );

		// CurrentRotationQuaternion = QuaternionType::FromRotationVector( PreviousRotationVector );
		// CurrentRotationQuaternion.ToRotationMatrix(PreviousRotationMatrix);

		// CurrentRotationMatrix = prod( CayleyStepRotation, PreviousRotationMatrix );

		// CurrentRotationQuaternion = QuaternionType::FromRotationMatrix( CurrentRotationMatrix );
		// CurrentRotationQuaternion.ToRotationVector(CurrentRotationVector); //assign

	      }
	      else{

		DeltaRotationQuaternion     = QuaternionType::FromRotationVector( DeltaRotationVector );

		// update step rotations:
		ReferenceRotationQuaternion = QuaternionType::FromRotationVector( CurrentStepRotationVector );
		
		CurrentRotationQuaternion   = DeltaRotationQuaternion * ReferenceRotationQuaternion;
		
		CurrentRotationQuaternion.ToRotationVector( CurrentStepRotationVector );
		
		// ---------------
				
		// update compound rotations:
		ReferenceRotationQuaternion = QuaternionType::FromRotationVector( CurrentRotationVector );
		
		CurrentRotationQuaternion   = DeltaRotationQuaternion * ReferenceRotationQuaternion;
		
		CurrentRotationQuaternion.ToRotationVector( CurrentRotationVector );
		
		// ---------------

	      }

	      for( unsigned int j=0; j<3; j++)
		{
		  //update rotation and step_rotation
		  CurrentRotation[j]      = CurrentRotationVector[j];
		  CurrentStepRotation[j]  = CurrentStepRotationVector[j];
		}	    


              array_1d<double, 3 > & CurrentAngularVelocity      = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
              array_1d<double, 3 > & CurrentAngularAcceleration  = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

	      array_1d<double, 3 > & PreviousAngularVelocity      = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY,1);
	      array_1d<double, 3 > & PreviousAngularAcceleration      = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION,1);

              
	      UpdateAngularVelocity     ((*i), CurrentAngularVelocity, PreviousAngularVelocity, CurrentStepRotation);
	      UpdateAngularAcceleration ((*i), CurrentAngularAcceleration, PreviousAngularAcceleration, CurrentAngularVelocity, PreviousAngularVelocity, CurrentStepRotation);

	    }

	  }

	
	//UPDATE SLAVE NODES
	SlaveNodesUpdate(r_model_part);


        KRATOS_CATCH( "" )
    }


    //***************************************************************************
    //***************************************************************************
    // this must be implemented in the ContactMechanicsApplication

    void SlaveNodesUpdate(ModelPart& r_model_part)
    {
        KRATOS_TRY

	// Matrix SkewSymVariable = ZeroMatrix(3,3);
        // Vector RadiusVector    = ZeroVector(3);
	// Vector Variable        = ZeroVector(3);
	// Vector AngularVariable = ZeroVector(3);
	// array_1d<double,3>     VariableArray;

        // array_1d<double,3> Radius;

	// for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
	//      i != r_model_part.NodesEnd(); ++i)
	//   {
	//     if( (i)->Is(SLAVE) ){

	//       Element& MasterElement = (i)->GetValue(MASTER_ELEMENTS).back();

	//       Node<3>& rCenterOfGravity = MasterElement.GetGeometry()[0];

	//       array_1d<double, 3 >&  Center              = rCenterOfGravity.GetInitialPosition();
	//       array_1d<double, 3 >&  Displacement        = rCenterOfGravity.FastGetSolutionStepValue(DISPLACEMENT);
	//       array_1d<double, 3 >&  Rotation            = rCenterOfGravity.FastGetSolutionStepValue(ROTATION);

	//       //Get rotation matrix
	//       QuaternionType TotalQuaternion = QuaternionType::FromRotationVector<array_1d<double,3> >(Rotation);

	//       Radius = (i)->GetInitialPosition() - Center;

	//       Matrix RotationMatrix;
	//       TotalQuaternion.ToRotationMatrix(RotationMatrix);
	      
	//       for(int j=0; j<3; j++)
	// 	RadiusVector[j] = Radius[j];
	      
	//       RadiusVector = prod( RotationMatrix, RadiusVector );

	//       for(int j=0; j<3; j++)
	// 	Radius[j] = RadiusVector[j];


	//       array_1d<double, 3 >&  NodeDisplacement      = (i)->FastGetSolutionStepValue(DISPLACEMENT);
	//       noalias(NodeDisplacement)     = ( (Center + Displacement)  + Radius ) - (i)->GetInitialPosition();


	//       //TotalQuaternion.RotateVector3<array_1d<double,3> >(Radius);
	//       if( (i)->IsNot(RIGID) && (i)->IsNot(MASTER) ){

	// 	//std::cout<<"  [ MasterElement "<<MasterElement.Id()<<" ]"<<std::endl;

	// 	array_1d<double, 3 >&  StepRotation        = rCenterOfGravity.FastGetSolutionStepValue(STEP_ROTATION);
	// 	array_1d<double, 3 >&  DeltaRotation       = rCenterOfGravity.FastGetSolutionStepValue(DELTA_ROTATION);
	// 	array_1d<double, 3 >&  Velocity            = rCenterOfGravity.FastGetSolutionStepValue(VELOCITY);
	// 	array_1d<double, 3 >&  Acceleration        = rCenterOfGravity.FastGetSolutionStepValue(ACCELERATION);
	// 	array_1d<double, 3 >&  AngularVelocity     = rCenterOfGravity.FastGetSolutionStepValue(ANGULAR_VELOCITY);
	// 	array_1d<double, 3 >&  AngularAcceleration = rCenterOfGravity.FastGetSolutionStepValue(ANGULAR_ACCELERATION);

	// 	//std::cout<<"  [ Rotation:"<<Rotation<<",StepRotation:"<<StepRotation<<",DeltaRotation:"<<DeltaRotation<<"]"<<std::endl;
	// 	//std::cout<<"  [ Velocity:"<<Velocity<<",Acceleration:"<<Acceleration<<",Displacement:"<<Displacement<<"]"<<std::endl;

	      
	// 	array_1d<double, 3 >&  NodeStepDisplacement  = (i)->FastGetSolutionStepValue(STEP_DISPLACEMENT);
	// 	array_1d<double, 3 >&  NodeRotation          = (i)->FastGetSolutionStepValue(ROTATION);
	// 	array_1d<double, 3 >&  NodeStepRotation      = (i)->FastGetSolutionStepValue(STEP_ROTATION);
	// 	array_1d<double, 3 >&  NodeDeltaRotation     = (i)->FastGetSolutionStepValue(DELTA_ROTATION);

	// 	noalias(NodeStepDisplacement) = NodeDisplacement - (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
	// 	noalias(NodeRotation)         = Rotation;
	// 	noalias(NodeStepRotation)     = StepRotation;
	// 	noalias(NodeDeltaRotation)    = DeltaRotation;    

	// 	for(int j=0; j<3; j++)
	// 	  RadiusVector[j] = Radius[j];

	// 	//********************
	// 	for(int j=0; j<3; j++)
	// 	  Variable[j] = AngularVelocity[j];

	// 	//compute the skewsymmmetric tensor of the angular velocity
	// 	BeamMathUtilsType::VectorToSkewSymmetricTensor(Variable, SkewSymVariable);
	      
	// 	//compute the contribution of the angular velocity to the velocity v = Wxr
	// 	Variable = prod(SkewSymVariable,RadiusVector);
     
	// 	for(int j=0; j<3; j++)
	// 	  VariableArray[j] = Variable[j];

	// 	(i)->FastGetSolutionStepValue(VELOCITY)               = Velocity + VariableArray;


	// 	//********************
	      
	// 	//centripetal acceleration:
	// 	for(int j=0; j<3; j++)
	// 	  AngularVariable[j] = AngularVelocity[j];
		
	// 	//compute the skewsymmmetric tensor of the angular velocity
	// 	BeamMathUtilsType::VectorToSkewSymmetricTensor(AngularVariable, SkewSymVariable);

	// 	AngularVariable = prod(SkewSymVariable,Variable); //ac = Wx(Wxr)


	// 	for(int j=0; j<3; j++)
	// 	  Variable[j] = AngularAcceleration[j];

	// 	//compute the skewsymmmetric tensor of the angular acceleration
	// 	BeamMathUtilsType::VectorToSkewSymmetricTensor(Variable, SkewSymVariable);
	      
	// 	//compute the contribution of the angular velocity to the velocity a = Axr
	// 	Variable = prod(SkewSymVariable,RadiusVector);

	// 	for(int j=0; j<3; j++)
	// 	  VariableArray[j] = Variable[j] + AngularVariable[j];

	// 	(i)->FastGetSolutionStepValue(ACCELERATION)           = Acceleration + VariableArray;


	// 	//********************
	// 	(i)->FastGetSolutionStepValue(ANGULAR_VELOCITY)       = AngularVelocity;
	// 	(i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION)   = AngularAcceleration;

	      
	//        	//std::cout<<"  [ Finalize Rigid Body Link Point : [Id:"<<(i)->Id()<<"] "<<std::endl;
	//        	//std::cout<<"  [ Rotation:"<<NodeRotation<<" / Angular Acceleration"<<AngularAcceleration<<" ] "<<std::endl;   
		
     
	//       }
		
	//     }
	//   }

	KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    //predicts the solution for the current step:
    // x = xold 

    void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY


 	//LINEAR VELOCITIES AND ACCELERATIONS

        //updating time derivatives (nodally for efficiency)
	array_1d<double, 3 >  DeltaDisplacement;

	for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                i != r_model_part.NodesEnd(); ++i)
	  {

	    if( i->IsNot(SLAVE) && i->IsNot(RIGID) ){
	    
	      array_1d<double, 3 > & CurrentDisplacement      = (i)->FastGetSolutionStepValue(DISPLACEMENT);
	      array_1d<double, 3 > & CurrentStepDisplacement  = (i)->FastGetSolutionStepValue(STEP_DISPLACEMENT);

	      unsigned int previous = 2;
              array_1d<double, 3 > & ReferenceDisplacement    = (i)->FastGetSolutionStepValue(DISPLACEMENT,previous);

	      noalias(CurrentStepDisplacement) = CurrentDisplacement - ReferenceDisplacement;

	      PredictStepDisplacement( (*i), CurrentStepDisplacement);

	      PredictDisplacement(CurrentDisplacement, CurrentStepDisplacement);

	      array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY);
	      array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION);
             
              array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY,previous);
	      array_1d<double, 3 > & PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,previous);

              UpdateVelocity     ((*i), CurrentVelocity, PreviousVelocity, CurrentStepDisplacement );
              UpdateAcceleration ((*i), CurrentAcceleration, PreviousAcceleration, CurrentVelocity, PreviousVelocity);
                            
	    }
	  }

        //updating time derivatives (nodally for efficiency)

	//UPDATE COMPOUND AND DELTA ROTATIONS BY QUATERNION COMPOSITION and ANGULAR VELOCITIES AND ACCELERATIONS
	Vector CurrentRotationVector         = ZeroVector(3);
	Vector PreviousRotationVector        = ZeroVector(3);
	Vector CurrentStepRotationVector     = ZeroVector(3);
	Vector DeltaRotationVector           = ZeroVector(3);


	QuaternionType DeltaRotationQuaternion;
	QuaternionType CurrentRotationQuaternion;
	QuaternionType ReferenceRotationQuaternion;


	for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
	     i != r_model_part.NodesEnd(); ++i)
	  {
	    if( i->IsNot(SLAVE) && i->IsNot(RIGID) ){

	      // Rotation
	      array_1d<double, 3 > & CurrentRotation     = (i)->FastGetSolutionStepValue(ROTATION);
	      // StepRotation
	      array_1d<double, 3 > & CurrentStepRotation = (i)->FastGetSolutionStepValue(STEP_ROTATION);

	      PredictStepRotation( (*i), CurrentStepRotation );

	      PredictRotation(CurrentRotation, CurrentStepRotation);   


              array_1d<double, 3 > & CurrentAngularAcceleration  = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION);
              array_1d<double, 3 > & CurrentAngularVelocity      = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY);

	      unsigned int previous = 2;
	      array_1d<double, 3 >&  PreviousAngularVelocity     = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY,previous);
	      array_1d<double, 3 >&  PreviousAngularAcceleration = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION,previous);

              
	      UpdateAngularVelocity     ((*i), CurrentAngularVelocity, PreviousAngularVelocity, CurrentStepRotation);
	      UpdateAngularAcceleration ((*i), CurrentAngularAcceleration, PreviousAngularAcceleration, CurrentAngularVelocity, PreviousAngularVelocity, CurrentStepRotation);

	    }

	  }

	
	//UPDATE SLAVE NODES
	SlaveNodesUpdate(r_model_part);


        KRATOS_CATCH( "" )
    }


    //***************************************************************************
    //***************************************************************************
    /**
    this is the place to initialize the elements.
    This is intended to be called just once when the strategy is initialized
     */
    void InitializeElements(ModelPart& rModelPart)
    {
        KRATOS_TRY

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin() + ElementPartition[k];
            ElementsArrayType::iterator ElemEnd = rModelPart.Elements().begin() + ElementPartition[k + 1];

            for (ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)
            {
                itElem->Initialize(); //function to initialize the element
            }

        }

        this->mElementsAreInitialized = true;

        KRATOS_CATCH( "" )
    }


    /**
    this is the place to initialize the conditions.
    This is intended to be called just once when the strategy is initialized
    */
    void InitializeConditions(ModelPart& rModelPart)
    {
        KRATOS_TRY

        if(this->mElementsAreInitialized==false)
            KRATOS_THROW_ERROR( std::logic_error, "Before initilizing Conditions, initialize Elements FIRST", "" )

        int NumThreads = OpenMPUtils::GetNumThreads();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Conditions().size(), NumThreads, ConditionPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin() + ConditionPartition[k];
            ConditionsArrayType::iterator CondEnd = rModelPart.Conditions().begin() + ConditionPartition[k + 1];

            for (ConditionsArrayType::iterator itCond = CondBegin; itCond != CondEnd; itCond++)
            {
                itCond->Initialize(); //function to initialize the condition
            }

        }

        this->mConditionsAreInitialized = true;
        KRATOS_CATCH( "" )
    }

    /**
     * initializes time step solution
     * only for reasons if the time step solution is restarted
     * @param r_model_part
     * @param A	LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        ProcessInfo& CurrentProcessInfo= r_model_part.GetProcessInfo();


        Scheme<TSparseSpace,TDenseSpace>::InitializeSolutionStep(r_model_part,A,Dx,b);

	
        double DeltaTime = CurrentProcessInfo[DELTA_TIME];

        if (DeltaTime == 0)
            KRATOS_THROW_ERROR(std::logic_error, "detected delta_time = 0 in the Solution Scheme ... check if the time step is created correctly for the current model part", "" )

	//--------------------------
	
	if( mDynamic.static_dynamic != 0 ){
	  CurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] = true;  
	  CurrentProcessInfo[EQUILIBRIUM_POINT] = 0.5;
	}
	else{
	  CurrentProcessInfo[EQUILIBRIUM_POINT] = 1;
	}

        //initializing Newmark constants
	mDynamic.deltatime = DeltaTime;
	mDynamic.c0 = 2.0 / DeltaTime;
	mDynamic.c1 = 1.0 / DeltaTime;

        KRATOS_CATCH( "" )
    }
    //***************************************************************************
    //***************************************************************************

    /**
    function called once at the end of a solution step, after convergence is reached if
    an iterative process is needed
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        //finalizes solution step for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ElementsArrayType::iterator ElementsBegin = rElements.begin() + ElementPartition[k];
            ElementsArrayType::iterator ElementsEnd = rElements.begin() + ElementPartition[k + 1];

            for (ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                itElem->FinalizeSolutionStep(CurrentProcessInfo);
            }
        }

        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ConditionsArrayType::iterator ConditionsBegin = rConditions.begin() + ConditionPartition[k];
            ConditionsArrayType::iterator ConditionsEnd = rConditions.begin() + ConditionPartition[k + 1];

            for (ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; itCond++)
            {
                itCond->FinalizeSolutionStep(CurrentProcessInfo);
            }
        }
        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void InitializeNonLinIteration(ModelPart& r_model_part,
                                   TSystemMatrixType& A,
                                   TSystemVectorType& Dx,
                                   TSystemVectorType& b)
    {
        KRATOS_TRY
        ElementsArrayType& pElements = r_model_part.Elements();
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for (ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            (it) -> InitializeNonLinearIteration(CurrentProcessInfo);
        }

        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            (it) -> InitializeNonLinearIteration(CurrentProcessInfo);
        }
        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void InitializeNonLinearIteration(Condition::Pointer rCurrentCondition,
                                      ProcessInfo& CurrentProcessInfo)
    {
        (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }


    //***************************************************************************
    //***************************************************************************

    void InitializeNonLinearIteration(Element::Pointer rCurrentElement,
                                      ProcessInfo& CurrentProcessInfo)
    {
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }



    //***************************************************************************
    //***************************************************************************

    /** this function is designed to be called in the builder and solver to introduce*/

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);


        if(mDynamic.static_dynamic !=0)
        {
	    
 	  (rCurrentElement) -> CalculateSecondDerivativesContributions(mMatrix.M[thread],mVector.fm[thread],CurrentProcessInfo);

	  (rCurrentElement) -> CalculateFirstDerivativesContributions(mMatrix.D[thread],mVector.fd[thread],CurrentProcessInfo);

        }


        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);


        if(mDynamic.static_dynamic !=0)
        {
	  
	  AddDynamicsToLHS(LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

	  AddDynamicsToRHS(RHS_Contribution, mVector.fd[thread], mVector.fm[thread], CurrentProcessInfo);

        }

        //AssembleTimeSpaceLHS(rCurrentElement, LHS_Contribution, DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {

        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        if(mDynamic.static_dynamic !=0)
        {

	  (rCurrentElement) -> CalculateSecondDerivativesRHS(mVector.fm[thread],CurrentProcessInfo);

	  (rCurrentElement) -> CalculateFirstDerivativesRHS(mVector.fd[thread],CurrentProcessInfo);

        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mDynamic.static_dynamic !=0)
        {

	  AddDynamicsToRHS(RHS_Contribution, mVector.fd[thread], mVector.fm[thread], CurrentProcessInfo);

        }

        KRATOS_CATCH( "" )

    }

    //***************************************************************************
    //***************************************************************************

    /** functions totally analogous to the precedent but applied to
          the "condition" objects
    */
    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {


        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        if(mDynamic.static_dynamic !=0)
        {

	  (rCurrentCondition) -> CalculateSecondDerivativesContributions(mMatrix.M[thread],mVector.fm[thread],CurrentProcessInfo);
	  
	  (rCurrentCondition) -> CalculateFirstDerivativesContributions(mMatrix.D[thread],mVector.fd[thread],CurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mDynamic.static_dynamic !=0)
        {

	  AddDynamicsToLHS(LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

	  AddDynamicsToRHS(RHS_Contribution, mVector.fd[thread], mVector.fm[thread], CurrentProcessInfo);

        }

        //AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }


    //***************************************************************************
    //***************************************************************************

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current condition
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        if(mDynamic.static_dynamic !=0)
        {

	  (rCurrentCondition) -> CalculateSecondDerivativesRHS(mVector.fm[thread],CurrentProcessInfo);

	  (rCurrentCondition) -> CalculateFirstDerivativesRHS(mVector.fd[thread],CurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

        //adding the dynamic contributions (static is already included)

        if(mDynamic.static_dynamic !=0)
        {
	  
	  AddDynamicsToRHS (RHS_Contribution, mVector.fd[thread], mVector.fm[thread], CurrentProcessInfo);

        }

        KRATOS_CATCH( "" )
    }


    //***************************************************************************
    //***************************************************************************

    /** Function that returns the list of Degrees of freedom to be
    assembled in the system for a Given Element
     */
    void GetElementalDofList(
        Element::Pointer rCurrentElement,
        Element::DofsVectorType& ElementalDofList,
        ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement->GetDofList(ElementalDofList, CurrentProcessInfo);
    }

    //***************************************************************************
    //***************************************************************************

    /** Function that returns the list of Degrees of freedom to be
    assembled in the system for a Given Element
     */
    void GetConditionDofList(
        Condition::Pointer rCurrentCondition,
        Element::DofsVectorType& ConditionDofList,
        ProcessInfo& CurrentProcessInfo)
    {
        rCurrentCondition->GetDofList(ConditionDofList, CurrentProcessInfo);
    }

    //***************************************************************************
    //***************************************************************************

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    virtual int Check(ModelPart& r_model_part)
    {
        KRATOS_TRY

        int err = Scheme<TSparseSpace, TDenseSpace>::Check(r_model_part);
        if(err!=0) return err;

        //check for variables keys
        //verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )

        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin();
                it!=r_model_part.NodesEnd(); it++)
        {
            if (it->SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if (it->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if (it->SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
        }

        //check that dofs exist
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin();
                it!=r_model_part.NodesEnd(); it++)
        {
            if(it->HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_X dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Y dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Z dof on node ",it->Id() )
        }


        //check for minimum value of the buffer index
        //verify buffer size
        if (r_model_part.GetBufferSize() < 2)
            KRATOS_THROW_ERROR( std::logic_error, "insufficient buffer size. Buffer size should be greater than 2. Current size is", r_model_part.GetBufferSize() )


        return 0;
        KRATOS_CATCH( "" )
    }


    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */

protected:


    GeneralDynamics     mDynamic;

    GeneralMatrices     mMatrix;
    GeneralVectors      mVector;


    /*@} */
    /**@name Protected Operators*/
    /*@{ */

   //*********************************************************************************
    //Predicting Step Displacement 
    //*********************************************************************************

    inline void PredictStepDisplacement(ModelPart::NodeType& Node,
					array_1d<double, 3 > & CurrentStepDisplacement)
    {
      
      //CurrentStepDisplacement.clear();

      if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
      	{
      	  CurrentStepDisplacement[0]  = 0;	  
      	}
      
      if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
      	{
      	  CurrentStepDisplacement[1]  = 0;	 
      	}
      
      if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
      	{
      	  CurrentStepDisplacement[2]  = 0;	  
      	}

    }

    //*********************************************************************************
    //Predicting displacement variable
    //*********************************************************************************

    inline void PredictDisplacement(array_1d<double, 3 > & CurrentDisplacement,
				    const array_1d<double, 3 > & StepDisplacement)
    {
      noalias(CurrentDisplacement) += StepDisplacement; //needed:: to capture the correct displacements
    }

   //*********************************************************************************
    //Predicting Step Rotation
    //*********************************************************************************

    inline void PredictStepRotation(ModelPart::NodeType& Node,
				    array_1d<double, 3 > & CurrentStepRotation)
    {

      //CurrentStepRotation.clear(); 

      if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
      	{
      	  CurrentStepRotation[0]  = 0;	  
      	}
      
      if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
      	{
      	  CurrentStepRotation[1]  = 0;	 
      	}
      
      if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
      	{
      	  CurrentStepRotation[2]  = 0;	  
      	}

    }

   //*********************************************************************************
    //Predicting rotation variable
    //*********************************************************************************

    inline void PredictRotation(array_1d<double, 3 > & CurrentRotation,
				const array_1d<double, 3 > & StepRotation)
    {
      noalias(CurrentRotation) += StepRotation; //needed for imposed rotations
    }

    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

    inline void UpdateVelocity(ModelPart::NodeType& Node,                                      
			       array_1d<double, 3 > & CurrentVelocity,
			       const array_1d<double, 3 > & PreviousVelocity,
			       const array_1d<double, 3 > & StepDisplacement)
    { 
      
      if(Node.Is(SLAVE)) return;

      if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false) 
      	{
      	  CurrentVelocity[0] = mDynamic.c0 * StepDisplacement[0] - PreviousVelocity[0];
      	}

      if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
      	{
      	  CurrentVelocity[1] = mDynamic.c0 * StepDisplacement[1] - PreviousVelocity[1];
      	}

      if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
      	{
      	  CurrentVelocity[2] = mDynamic.c0 * StepDisplacement[2] - PreviousVelocity[2];
      }
      
      CurrentVelocity *= mDynamic.static_dynamic;

     
    }
    
    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************

    
    inline void UpdateAcceleration(ModelPart::NodeType& Node,
				   array_1d<double, 3 > & CurrentAcceleration,
				   array_1d<double, 3 > & PreviousAcceleration,
				   const array_1d<double, 3 > & CurrentVelocity,
				   const array_1d<double, 3 > & PreviousVelocity)
    {   

      if(Node.Is(SLAVE)) return;        
      

      if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
      	{
      	  CurrentAcceleration[0] =  mDynamic.c1 * (CurrentVelocity[0] - PreviousVelocity[0]);
      	}
 
      if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
      	{
      	  CurrentAcceleration[1] =  mDynamic.c1 * (CurrentVelocity[1] - PreviousVelocity[1]);
      	}

      if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
      	{
      	  CurrentAcceleration[2] =  mDynamic.c1 * (CurrentVelocity[2] - PreviousVelocity[2]);
      	}

      // if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
      // 	{
      // 	  CurrentAcceleration[0] =  mDynamic.c1 * CurrentVelocity[0] - PreviousAcceleration[0];
      // 	}
 
      // if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
      // 	{
      // 	  CurrentAcceleration[1] =  mDynamic.c1 * CurrentVelocity[1] - PreviousAcceleration[1];
      // 	}

      // if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
      // 	{
      // 	  CurrentAcceleration[2] =  mDynamic.c1 * CurrentVelocity[2] - PreviousAcceleration[2];
      // 	}

      CurrentAcceleration *= mDynamic.static_dynamic;

    }
    
    


    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

 
    inline void UpdateAngularVelocity(ModelPart::NodeType& Node,
                                      array_1d<double, 3 > & CurrentAngularVelocity,
				      const array_1d<double, 3 > & PreviousAngularVelocity,
				      const array_1d<double, 3 > & StepRotation)
    { 

      if(Node.Is(SLAVE)) return;        
      
      array_1d<double, 3 > AngularVelocity = PreviousAngularVelocity;

      ///-------------------

      Vector CurrentStepRotationVector = ZeroVector(3);
      Vector PreviousAngularVelocityVector = ZeroVector(3);
      

      for( unsigned int j=0; j<3; j++)
      	{
      	  CurrentStepRotationVector[j]     = StepRotation[j];   //previous iteration step rotation 
      	  PreviousAngularVelocityVector[j] = PreviousAngularVelocity[j];
      	}
      
      Matrix CayleyStepRotation  = ZeroMatrix(3,3);

      BeamMathUtilsType::CayleyTransform( CurrentStepRotationVector, CayleyStepRotation );      
 
      PreviousAngularVelocityVector = prod(CayleyStepRotation, PreviousAngularVelocityVector);

      for( unsigned int j=0; j<3; j++)
      	{
      	  AngularVelocity[j] = PreviousAngularVelocityVector[j];
      	}

     ///-------------------

      if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
      	{
      	  CurrentAngularVelocity[0] =  mDynamic.c0 * StepRotation[0] - AngularVelocity[0];
      	}

      if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
      	{
      	  CurrentAngularVelocity[1] =  mDynamic.c0 * StepRotation[1] - AngularVelocity[1];
      	}

      if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
      	{
      	  CurrentAngularVelocity[2] =  mDynamic.c0 * StepRotation[2] - AngularVelocity[2];
	}

      CurrentAngularVelocity *= mDynamic.static_dynamic;

      // PROBLEM IN THE INITIAL CONDITIONS: the first aproach is very bad.
      // if( AngularVelocity[0] == 0 && AngularVelocity[1] == 0 && AngularVelocity[2] == 0 )
      //  	CurrentAngularVelocity *= 0.5;

      //std::cout<<" CurrentAngularVelocity "<<CurrentAngularVelocity<<std::endl;

    }


    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************

    
    inline void UpdateAngularAcceleration(ModelPart::NodeType& Node,
                                          array_1d<double, 3 > & CurrentAngularAcceleration,
                                          array_1d<double, 3 > & PreviousAngularAcceleration,					  
                                          const array_1d<double, 3 > & CurrentAngularVelocity,
                                          const array_1d<double, 3 > & PreviousAngularVelocity,
					  const array_1d<double, 3 > & StepRotation)
    {           


      if(Node.Is(SLAVE)) return;        
      
      array_1d<double, 3 > AngularVelocity = PreviousAngularVelocity;

      ///-------------------

      Vector CurrentStepRotationVector = ZeroVector(3);
      Vector PreviousAngularVelocityVector = ZeroVector(3);
      

      for( unsigned int j=0; j<3; j++)
      	{
      	  CurrentStepRotationVector[j] = StepRotation[j];   //previous iteration step rotation 
      	  PreviousAngularVelocityVector[j] = PreviousAngularVelocity[j];
      	}
      
      Matrix CayleyStepRotation  = ZeroMatrix(3,3);

      BeamMathUtilsType::CayleyTransform( CurrentStepRotationVector, CayleyStepRotation );

      PreviousAngularVelocityVector = prod(CayleyStepRotation, PreviousAngularVelocityVector);

     
      for( unsigned int j=0; j<3; j++)
      	{
      	  AngularVelocity[j] = PreviousAngularVelocityVector[j];
      	}

     ///-------------------

      if(Node.Is(SLAVE)) return;        

      if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
      	{
      	  CurrentAngularAcceleration[0] =  mDynamic.c1 * ( CurrentAngularVelocity[0] - AngularVelocity[0] );
      	}

      if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
      	{
      	  CurrentAngularAcceleration[1] =  mDynamic.c1 * ( CurrentAngularVelocity[1] - AngularVelocity[1] );
      	}

      if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
      	{
      	  CurrentAngularAcceleration[2] =  mDynamic.c1 * ( CurrentAngularVelocity[2] - AngularVelocity[2] );
      	}

      // if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularAcceleration[0] =  mDynamic.c1 * CurrentAngularVelocity[0] - AngularAcceleration[0];
      // 	}

      // if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularAcceleration[1] =  mDynamic.c1 * CurrentAngularVelocity[1] - AngularAcceleration[1];
      // 	}

      // if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularAcceleration[2] =  mDynamic.c1 * CurrentAngularVelocity[2] - AngularAcceleration[2];
      // 	}

      CurrentAngularAcceleration *= mDynamic.static_dynamic;
    }


    //Elements and Conditions:
    //****************************************************************************

    /**
    Atangent = M*c0 + D*c1 + K

     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {

        // adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += M * mDynamic.static_dynamic;
        }

        //adding  damping contribution
        if (D.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += D * mDynamic.static_dynamic;
        }

    }

    /**
    bdyn = b - M*a - D*v

     */
     void AddDynamicsToRHS(
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemVectorType& fd,
        LocalSystemVectorType& fm,
        ProcessInfo& CurrentProcessInfo)
    {

        //adding inertia contribution
        if (fm.size() != 0)
        {
            noalias(RHS_Contribution) -=  mDynamic.static_dynamic * fm;
        }

        //adding damping contribution
        if (fd.size() != 0)
        {
	    noalias(RHS_Contribution) -=  mDynamic.static_dynamic * fd;
        }


    }


    //Elements:
    //****************************************************************************

    /**
    bdyn = b - M*a - D*v

     */
    void AddDynamicsToRHS(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
	    rCurrentElement->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= mDynamic.static_dynamic ;

            rCurrentElement->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mVector.ap[thread] * mDynamic.static_dynamic;

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

        //adding damping contribution
        if (D.size1() != 0)
        {
            rCurrentElement->GetFirstDerivativesVector(mVector.v[thread], 0);

            (mVector.v[thread]) *= mDynamic.static_dynamic ;

            noalias(RHS_Contribution) -= prod(D, mVector.v[thread]);
        }


    }



    //Conditions:
    //****************************************************************************


    void AddDynamicsToRHS(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentCondition->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *=  mDynamic.static_dynamic;

            rCurrentCondition->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) +=  mVector.ap[thread] * mDynamic.static_dynamic;

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

        //adding damping contribution
        if (D.size1() != 0)
        {
            rCurrentCondition->GetFirstDerivativesVector(mVector.v[thread], 0);

            (mVector.v[thread]) *= mDynamic.static_dynamic ;

            noalias(RHS_Contribution) -= prod(D, mVector.v [thread]);
        }

    }



    /*@} */
    /**@name Protected member Variables */
    /*@{ */
    /*@} */
    /**@name Protected Operators*/
    /*@{ */
    /*@} */
    /**@name Protected Operations*/
    /*@{ */
    /*@} */
    /**@name Protected  Access */
    /*@{ */
    /*@} */
    /**@name Protected Inquiry */
    /*@{ */
    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */
private:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */
    //DofsVectorType mElementalDofList;
    /*@} */
    /**@name Private Operators*/
    /*@{ */
    /*@} */
    /**@name Private Operations*/
    /*@{ */
    /*@} */
    /**@name Private  Access */
    /*@{ */
    /*@} */
    /**@name Private Inquiry */
    /*@{ */
    /*@} */
    /**@name Unaccessible methods */
    /*@{ */
}; /* Class Scheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ROTATION_CONSERVATIVE_SCHEME  defined */


