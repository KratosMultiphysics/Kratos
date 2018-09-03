//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:             October 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_ROTATION_SCHEME)
#define KRATOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_ROTATION_SCHEME

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
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

// Covariant implicit time stepping algorithm: the classical Newmark algorithm of nonlinear elastodynamics and a canonical extension of the Newmark formulas to the orthogonal group SO(3) for the rotational part.

template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedBossakDisplacementRotationScheme: public Scheme<TSparseSpace,TDenseSpace>
{
protected:


    struct GeneralAlphaMethod
    {
        double f;  //alpha Hilbert
        double m;  //alpha Bosssak

    };

    struct NewmarkMethod
    {

        double beta;
        double gamma;
        double deltatime;

        //system constants
        double c0;
        double c1;
        double c2;
        double c3;
        double c4;
        double c5;
        double c6;
        double c7;

        //static-dynamic parameter
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
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedBossakDisplacementRotationScheme );

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
    ResidualBasedBossakDisplacementRotationScheme(double rDynamic = 1, double mAlphaM = 0.0)
        :Scheme<TSparseSpace,TDenseSpace>()
    {
        //For pure stable Newmark Scheme
        // mNewmark.beta=  0.25;
        // mNewmark.gamma= 0.5;

	//For Bossak Scheme
	mAlpha.f= 0;
        mAlpha.m= mAlphaM; //0.0 to -0.3

        mNewmark.beta= (1.0+mAlpha.f-mAlpha.m)*(1.0+mAlpha.f-mAlpha.m)*0.25;
        mNewmark.gamma= 0.5+mAlpha.f-mAlpha.m;


        mNewmark.static_dynamic= rDynamic;


        //std::cout << " MECHANICAL SCHEME: The Newmark Simo Time Integration Scheme [ beta= "<<mNewmark.beta<<" gamma= "<<mNewmark.gamma<<"]"<<std::endl;


        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();

        mMatrix.M.resize(NumThreads);
        mMatrix.D.resize(NumThreads);

        mVector.fm.resize(NumThreads);
	mVector.fd.resize(NumThreads);

    }

    /** Destructor.
     */
    virtual ~ResidualBasedBossakDisplacementRotationScheme
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

	      //STEP_DISPLACEMENT variable used to store previous iteration displacement
	      array_1d<double, 3 > & CurrentDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);
	      array_1d<double, 3 > & ReferenceDisplacement= (i)->FastGetSolutionStepValue(STEP_DISPLACEMENT, 1);

	      // Displacement at iteration i
	      ReferenceDisplacement = CurrentDisplacement;

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
	array_1d<double, 3 >  DeltaDisplacement;

	for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                i != r_model_part.NodesEnd(); ++i)
	  {

	    if( i->IsNot(SLAVE) && i->IsNot(RIGID) ){

	      array_1d<double, 3 > & CurrentDisplacement      = (i)->FastGetSolutionStepValue(DISPLACEMENT);
	      array_1d<double, 3 > & ReferenceDisplacement    = (i)->FastGetSolutionStepValue(STEP_DISPLACEMENT, 1); //the CurrentDisplacement of the last iteration was stored here

	      noalias(DeltaDisplacement) = CurrentDisplacement-ReferenceDisplacement;

	      array_1d<double, 3 > & CurrentStepDisplacement  = (i)->FastGetSolutionStepValue(STEP_DISPLACEMENT);

	      noalias(CurrentStepDisplacement)  = CurrentDisplacement - (i)->FastGetSolutionStepValue(DISPLACEMENT,1);

	      //CurrentStepDisplacement += DeltaDisplacement;

	      array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY);
	      array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION);

              array_1d<double, 3 > & PreviousVelocity      = (i)->FastGetSolutionStepValue(VELOCITY,1);
              array_1d<double, 3 > & PreviousAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION,1);

              UpdateAcceleration ((*i), CurrentAcceleration, DeltaDisplacement, CurrentStepDisplacement, PreviousAcceleration, PreviousVelocity);
              UpdateVelocity     ((*i), CurrentVelocity, DeltaDisplacement, CurrentAcceleration, PreviousAcceleration, PreviousVelocity);

	    }
	  }

        //updating time derivatives (nodally for efficiency)

	//UPDATE COMPOUND AND DELTA ROTATIONS BY QUATERNION COMPOSITION and ANGULAR VELOCITIES AND ACCELERATIONS
	Vector TotalRotationVector = ZeroVector(3);
	Vector StepRotationVector  = ZeroVector(3);
	Vector DeltaRotationVector = ZeroVector(3);

	Vector LinearDeltaRotationVector = ZeroVector(3);

	QuaternionType DeltaRotationQuaternion;
	QuaternionType CurrentRotationQuaternion;
	QuaternionType ReferenceRotationQuaternion;


	for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
	     i != r_model_part.NodesEnd(); ++i)
	  {
	    if( i->IsNot(SLAVE) && i->IsNot(RIGID) ){

	      // Rotation at iteration i+1
	      array_1d<double, 3 > & CurrentRotation     = (i)->FastGetSolutionStepValue(ROTATION);
	      // Rotation at iteration i
	      array_1d<double, 3 > & ReferenceRotation   = (i)->FastGetSolutionStepValue(DELTA_ROTATION, 1);

	      //StepRotation
	      array_1d<double, 3 > & CurrentStepRotation = (i)->FastGetSolutionStepValue(STEP_ROTATION);

	      //DeltaRotation
	      array_1d<double, 3 > & CurrentDeltaRotation= (i)->FastGetSolutionStepValue(DELTA_ROTATION);

	      CurrentDeltaRotation = CurrentRotation-ReferenceRotation;

	      // if( (i)->Id()==26 )
	      //   {
	      // 	  std::cout<<" CurrentDeltaRotation "<<CurrentDeltaRotation<<std::endl;
	      // 	  std::cout<<" ReferenceDeltaRotation "<<ReferenceRotation<<" CurrentRotation "<<CurrentRotation<<std::endl;
	      // 	}


	      for( unsigned int j=0; j<3; j++)
		{
		  DeltaRotationVector[j]  = CurrentDeltaRotation[j]; //delta rotation computed in the iteration

		  TotalRotationVector[j]  = ReferenceRotation[j];    //previous iteration total rotation

		  StepRotationVector[j]   = CurrentStepRotation[j];  //previous iteration step rotation
		}


	      DeltaRotationQuaternion = QuaternionType::FromRotationVector( DeltaRotationVector );

	      // updated step rotations:
	      ReferenceRotationQuaternion = QuaternionType::FromRotationVector( StepRotationVector );

	      CurrentRotationQuaternion = DeltaRotationQuaternion * ReferenceRotationQuaternion;

	      CurrentRotationQuaternion.ToRotationVector( StepRotationVector );

	      // ---------------

	      // std::cout<<" Step    Rotation B: "<<StepRotationVector<<std::endl;

	      // updated compound rotations:
	      ReferenceRotationQuaternion = QuaternionType::FromRotationVector( TotalRotationVector );

	      CurrentRotationQuaternion = DeltaRotationQuaternion * ReferenceRotationQuaternion;

	      CurrentRotationQuaternion.ToRotationVector( TotalRotationVector );

	      // ---------------

	      for( unsigned int j=0; j<3; j++)
		{
		  LinearDeltaRotationVector[j] = StepRotationVector[j] - CurrentStepRotation[j];
		}

	      // update delta rotation composed:
	      Matrix RotationMatrix;
	      (ReferenceRotationQuaternion.conjugate()).ToRotationMatrix(RotationMatrix);
	      LinearDeltaRotationVector = prod( RotationMatrix, LinearDeltaRotationVector );

	      //(ReferenceRotationQuaternion.conjugate()).RotateVector3(LinearDeltaRotationVector); //exp[-X_n]*(Rotation_i+1-Rotation_i)
	      //(CurrentRotationQuaternion.conjugate()).RotateVector3(LinearDeltaRotationVector); //exp[-X_n+1]*(Rotation_i+1-Rotation_i)

	      array_1d<double, 3 > LinearDeltaRotation;

	      for( unsigned int j=0; j<3; j++)
		{
		  //update rotation and step_rotation
		  CurrentRotation[j]      = TotalRotationVector[j];
		  CurrentStepRotation[j]  = StepRotationVector[j];

		  //array_1d to pass in the acceleration and velocity updates
		  LinearDeltaRotation[j]  = LinearDeltaRotationVector[j];
		}


	      array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
	      array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION);


              array_1d<double, 3 > & PreviousAngularAcceleration  = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION,1);
              array_1d<double, 3 > & PreviousAngularVelocity      = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY,1);


              UpdateAngularAcceleration ((*i), CurrentAcceleration, LinearDeltaRotation, CurrentStepRotation, PreviousAngularAcceleration, PreviousAngularVelocity);
	      UpdateAngularVelocity     ((*i), CurrentVelocity, LinearDeltaRotation, CurrentAcceleration, PreviousAngularAcceleration, PreviousAngularVelocity);


	    }

	  }


	//UPDATE SLAVE NODES
	SlaveNodesUpdate(r_model_part);


        KRATOS_CATCH( "" )
    }




    //***************************************************************************
    //***************************************************************************

    void SlaveNodesUpdate(ModelPart& r_model_part)
    {
        KRATOS_TRY

	Matrix SkewSymVariable = ZeroMatrix(3,3);
        Vector RadiusVector    = ZeroVector(3);
	Vector Variable        = ZeroVector(3);
	Vector AngularVariable = ZeroVector(3);
	array_1d<double,3>     VariableArray;

        array_1d<double,3> Radius;

	for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
	     i != r_model_part.NodesEnd(); ++i)
	  {
	    if( (i)->Is(SLAVE) && i->IsNot(RIGID) ){

	      Element& MasterElement = (i)->GetValue(MASTER_ELEMENTS).back();

	      Node<3>& rCenterOfGravity = MasterElement.GetGeometry()[0];

	      array_1d<double, 3 >&  Center              = rCenterOfGravity.GetInitialPosition();
	      array_1d<double, 3 >&  Displacement        = rCenterOfGravity.FastGetSolutionStepValue(DISPLACEMENT);
	      array_1d<double, 3 >&  Rotation            = rCenterOfGravity.FastGetSolutionStepValue(ROTATION);
	      array_1d<double, 3 >&  StepRotation        = rCenterOfGravity.FastGetSolutionStepValue(STEP_ROTATION);
	      array_1d<double, 3 >&  DeltaRotation       = rCenterOfGravity.FastGetSolutionStepValue(DELTA_ROTATION);

	      array_1d<double, 3 >&  Velocity            = rCenterOfGravity.FastGetSolutionStepValue(VELOCITY);
	      array_1d<double, 3 >&  Acceleration        = rCenterOfGravity.FastGetSolutionStepValue(ACCELERATION);
	      array_1d<double, 3 >&  AngularVelocity     = rCenterOfGravity.FastGetSolutionStepValue(ANGULAR_VELOCITY);
	      array_1d<double, 3 >&  AngularAcceleration = rCenterOfGravity.FastGetSolutionStepValue(ANGULAR_ACCELERATION);

	      //std::cout<<" [  MasterElement "<<MasterElement.Id() ];
	      //std::cout<<" [ Rotation:"<<Rotation<<",StepRotation:"<<StepRotation<<",DeltaRotation:"<<DeltaRotation<<"]"<<std::endl;
	      //std::cout<<" [ Velocity:"<<Velocity<<",Acceleration:"<<Acceleration<<",Displacement:"<<Displacement<<",DeltaDisplacement"<<Displacement-rCenterOfGravity.FastGetSolutionStepValue(DISPLACEMENT,1)<<"]"<<std::endl;

	      //Get rotation matrix
	      QuaternionType TotalQuaternion = QuaternionType::FromRotationVector<array_1d<double,3> >(Rotation);

	      Radius = (i)->GetInitialPosition() - Center;

	      Matrix RotationMatrix;
	      TotalQuaternion.ToRotationMatrix(RotationMatrix);

	      for(int j=0; j<3; j++)
		RadiusVector[j] = Radius[j];

	      RadiusVector = prod( RotationMatrix, RadiusVector );

	      for(int j=0; j<3; j++)
		Radius[j] = RadiusVector[j];

	      //TotalQuaternion.RotateVector3<array_1d<double,3> >(Radius);

	      array_1d<double, 3 >&  NodeDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);
	      array_1d<double, 3 >&  NodeRotation      = (i)->FastGetSolutionStepValue(ROTATION);
	      array_1d<double, 3 >&  NodeStepRotation  = (i)->FastGetSolutionStepValue(STEP_ROTATION);
	      array_1d<double, 3 >&  NodeDeltaRotation = (i)->FastGetSolutionStepValue(DELTA_ROTATION);

	      noalias(NodeDisplacement)  = ( (Center + Displacement)  + Radius ) - (i)->GetInitialPosition();
	      noalias(NodeRotation)      = Rotation;
	      noalias(NodeStepRotation)  = StepRotation;
	      noalias(NodeDeltaRotation) = DeltaRotation;

	      for(int j=0; j<3; j++)
		RadiusVector[j] = Radius[j];

	      //********************

	      for(int j=0; j<3; j++)
		Variable[j] = AngularVelocity[j];

	      //compute the skewsymmmetric tensor of the angular velocity
	      BeamMathUtilsType::VectorToSkewSymmetricTensor(Variable, SkewSymVariable);

	      //compute the contribution of the angular velocity to the velocity v = Wxr
	      Variable = prod(SkewSymVariable,RadiusVector);

	      for(int j=0; j<3; j++)
		VariableArray[j] = Variable[j];

	      (i)->FastGetSolutionStepValue(VELOCITY)               = Velocity + VariableArray;

	      //********************

	      //centripetal acceleration:
	      for(int j=0; j<3; j++)
		AngularVariable[j] = AngularVelocity[j];

	      //compute the skewsymmmetric tensor of the angular velocity
	      BeamMathUtilsType::VectorToSkewSymmetricTensor(AngularVariable, SkewSymVariable);

	      AngularVariable = prod(SkewSymVariable,Variable); //ac = Wx(Wxr)

	      for(int j=0; j<3; j++)
		Variable[j] = AngularAcceleration[j];

	      //compute the skewsymmmetric tensor of the angular acceleration
	      BeamMathUtilsType::VectorToSkewSymmetricTensor(Variable, SkewSymVariable);

	      //compute the contribution of the angular velocity to the velocity a = Axr
	      Variable = prod(SkewSymVariable,RadiusVector);

     	      for(int j=0; j<3; j++)
		VariableArray[j] = Variable[j] + AngularVariable[j];

	      (i)->FastGetSolutionStepValue(ACCELERATION)           = Acceleration + VariableArray;


	      //********************
	      (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY)       = AngularVelocity;
	      (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION)   = AngularAcceleration;

	      // 	std::cout<<"  [ Finalize Rigid Body Link Point : [Id:"<<(i)->Id()<<"] "<<std::endl;
	      // 	std::cout<<"  [ Displacement:"<<NodeDisplacement<<" / StepRotation"<<NodeStepRotation<<" ] "<<std::endl;
	      // 	std::cout<<"  [ Rotation:"<<NodeRotation<<" / Angular Acceleration"<<AngularAcceleration<<" ] "<<std::endl;


	    }

	  }

	KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    //predicts the solution for the current step:
    // x = xold + vold * Dt

    void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {

        KRATOS_TRY

        //std::cout << " Prediction " << std::endl;

        //double DeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];

	//DISPLACEMENTS

        for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                i != r_model_part.NodesEnd(); ++i)
	  {

	    if( i->IsNot(SLAVE) && i->IsNot(RIGID) ){
	      //predicting displacement = PreviousDisplacement + PreviousVelocity * DeltaTime;
	      //ATTENTION::: the prediction is performed only on free nodes


	      array_1d<double, 3 > & CurrentDisplacement      = (i)->FastGetSolutionStepValue(DISPLACEMENT);
	      array_1d<double, 3 > & CurrentVelocity          = (i)->FastGetSolutionStepValue(VELOCITY);
	      array_1d<double, 3 > & CurrentAcceleration      = (i)->FastGetSolutionStepValue(ACCELERATION);

	      array_1d<double, 3 > & CurrentStepDisplacement  = (i)->FastGetSolutionStepValue(STEP_DISPLACEMENT);


              //DISPLACEMENT AND STEP DISPLACEMENT
	      unsigned int previous = 1;
              array_1d<double, 3 > & ReferenceDisplacement    = (i)->FastGetSolutionStepValue(DISPLACEMENT,previous);
	      array_1d<double, 3 > & PreviousVelocity         = (i)->FastGetSolutionStepValue(VELOCITY,previous);
	      array_1d<double, 3 > & PreviousAcceleration     = (i)->FastGetSolutionStepValue(ACCELERATION,previous);

	      noalias(CurrentStepDisplacement) = CurrentDisplacement - ReferenceDisplacement;

              //only used in a special prediction of the displacement
	      array_1d<double, 3 > & CurrentAngularVelocity   = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY);

	      PredictStepDisplacement( (*i),CurrentStepDisplacement, CurrentVelocity, PreviousAcceleration, CurrentAcceleration, CurrentAngularVelocity);

	      if (i->HasDofFor(LAGRANGE_MULTIPLIER_NORMAL))
		{
		  double& PreviousLM    = (i)->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_NORMAL, 1);
		  double& CurrentLM     = (i)->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_NORMAL);

		  CurrentLM = PreviousLM;
		}

	      //DISPLACEMENT, LINEAR VELOCITIES AND ACCELERATIONS

	      PredictDisplacement(CurrentDisplacement, CurrentStepDisplacement);

	      //updating time derivatives ::: please note that displacements and its time derivatives
	      //can not be consistently fixed separately

	      PredictAcceleration ( (*i), CurrentAcceleration, CurrentStepDisplacement, PreviousAcceleration, PreviousVelocity);
	      PredictVelocity     ( (*i), CurrentVelocity, CurrentStepDisplacement, CurrentAcceleration, PreviousAcceleration, PreviousVelocity);

	      //std::cout<<" Prediction: Velocity "<<CurrentVelocity<<" Accel "<<CurrentAcceleration<<" Disp "<<CurrentDisplacement<<" StepDisp "<<CurrentStepDisplacement<<std::endl;
	    }

	  }


	//ROTATIONS

        for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                i != r_model_part.NodesEnd(); ++i)
	  {

	    if( i->IsNot(SLAVE) && i->IsNot(RIGID) ){

	      //ATTENTION::: the prediction is performed only on free nodes
	      array_1d<double, 3 >& CurrentRotation              = (i)->FastGetSolutionStepValue(ROTATION);
	      array_1d<double, 3 >& CurrentStepRotation          = (i)->FastGetSolutionStepValue(STEP_ROTATION);
	      array_1d<double, 3 >& CurrentAngularVelocity       = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
	      array_1d<double, 3 >& CurrentAngularAcceleration   = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

	      //array_1d<double, 3 >& ImposedRotation              = (i)->FastGetSolutionStepValue(IMPOSED_ROTATION);

	      //ROTATIONS AND STEP ROTATIONS
	      unsigned int previous = 1;
	      array_1d<double, 3 >  PreviousAngularVelocity     = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY,previous);
	      array_1d<double, 3 >  PreviousAngularAcceleration = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION,previous);

	      PredictStepRotation( (*i), CurrentStepRotation, CurrentAngularVelocity, PreviousAngularAcceleration, CurrentAngularAcceleration);

	      //ROTATIONS, ANGULAR VELOCITIES AND ACCELERATIONS

	      //updating time derivatives ::: please note that rotations and its time derivatives
	      //can not be consistently fixed separately

	      PredictRotation(CurrentRotation, CurrentStepRotation);


	      PredictAngularAcceleration ( (*i), CurrentAngularAcceleration, CurrentStepRotation, PreviousAngularAcceleration,  PreviousAngularVelocity);
	      PredictAngularVelocity     ( (*i), CurrentAngularVelocity, CurrentStepRotation, CurrentAngularAcceleration, PreviousAngularAcceleration, PreviousAngularVelocity);

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
        //std::cout<<" mechanical elements are initialized "<<std::endl;

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

	//Set Newmark coefficients

	if( mNewmark.static_dynamic != 0 ){
	  CurrentProcessInfo[NEWMARK_BETA]  = mNewmark.beta;
	  CurrentProcessInfo[NEWMARK_GAMMA] = mNewmark.gamma;
	  CurrentProcessInfo[BOSSAK_ALPHA]  = mAlpha.m;
	  CurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] = true;
	}

        //initializing Newmark constants
	mNewmark.deltatime = DeltaTime;

        mNewmark.c0 = ( mNewmark.gamma / ( DeltaTime * mNewmark.beta ) );
        mNewmark.c1 = ( 1.0 / ( DeltaTime * DeltaTime * mNewmark.beta ) );

        mNewmark.c2 = ( DeltaTime * ( 1.0 - mNewmark.gamma ) );
        mNewmark.c3 = ( DeltaTime * mNewmark.gamma );
        mNewmark.c4 = ( DeltaTime / mNewmark.beta );
        mNewmark.c5 = ( DeltaTime * DeltaTime * ( 0.5 - mNewmark.beta ) / mNewmark.beta );

	//extra
	mNewmark.c6 = ( 1.0 / (mNewmark.beta * DeltaTime) );
	mNewmark.c7 = ( 0.5 / (mNewmark.beta) - 1.0 );

        //std::cout<<" Newmark Variables "<<mNewmark.c0<<" "<<mNewmark.c1<<" "<<mNewmark.c2<<" "<<mNewmark.c3<<" "<<mNewmark.c4<<std::endl;


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


        if(mNewmark.static_dynamic !=0)
        {

 	  (rCurrentElement) -> CalculateSecondDerivativesContributions(mMatrix.M[thread],mVector.fm[thread],CurrentProcessInfo);

	  (rCurrentElement) -> CalculateFirstDerivativesContributions(mMatrix.D[thread],mVector.fd[thread],CurrentProcessInfo);

        }


        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);


        if(mNewmark.static_dynamic !=0)
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

        if(mNewmark.static_dynamic !=0)
        {

	  (rCurrentElement) -> CalculateSecondDerivativesRHS(mVector.fm[thread],CurrentProcessInfo);

	  (rCurrentElement) -> CalculateFirstDerivativesRHS(mVector.fd[thread],CurrentProcessInfo);

        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
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

        if(mNewmark.static_dynamic !=0)
        {

	  (rCurrentCondition) -> CalculateSecondDerivativesContributions(mMatrix.M[thread],mVector.fm[thread],CurrentProcessInfo);

	  (rCurrentCondition) -> CalculateFirstDerivativesContributions(mMatrix.D[thread],mVector.fd[thread],CurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
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

        if(mNewmark.static_dynamic !=0)
        {

	  (rCurrentCondition) -> CalculateSecondDerivativesRHS(mVector.fm[thread],CurrentProcessInfo);

	  (rCurrentCondition) -> CalculateFirstDerivativesRHS(mVector.fd[thread],CurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

        //adding the dynamic contributions (static is already included)

        if(mNewmark.static_dynamic !=0)
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

    GeneralAlphaMethod  mAlpha;
    NewmarkMethod       mNewmark;

    GeneralMatrices     mMatrix;
    GeneralVectors      mVector;


    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    //*********************************************************************************
    //Predicting Step Displacement variable
    //*********************************************************************************

    inline void PredictStepDisplacement(ModelPart::NodeType& Node,
					array_1d<double, 3 > & CurrentStepDisplacement,
					const array_1d<double, 3 > & PreviousVelocity,
					const array_1d<double, 3 > & PreviousAcceleration,
					const array_1d<double, 3 > & CurrentAcceleration,
					const array_1d<double, 3 > & CurrentAngularVelocity)
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

	//noalias(CurrentStepDisplacement) = ( mNewmark.c4 * PreviousVelocity + mNewmark.c5 * PreviousAcceleration + CurrentAcceleration) * mNewmark.beta *  mNewmark.static_dynamic;

        //noalias(CurrentStepDisplacement) = ( mNewmark.deltatime * PreviousVelocity ) *  mNewmark.static_dynamic;

	// // Angular velocity correction:

       	// QuaternionType AngularVelocityQuaternion;

	// Vector AngularVelocityVector = ZeroVector(3);

	// for( unsigned int j=0; j<3; j++)
	//   {
	//     AngularVelocityVector[j]  =  mNewmark.c4 * mNewmark.beta * CurrentAngularVelocity[j];
	//   }

	// //updated total rotations:
	// AngularVelocityQuaternion = QuaternionType::FromRotationVector( CurrentAngularVelocity );

	// Vector CurrentStepDisplacementVector   = ZeroVector(3);
	// for( unsigned int j=0; j<3; j++)
	//   {
	//     CurrentStepDisplacementVector[j]  = CurrentStepDisplacement[j];
	//   }

	// AngularVelocityQuaternion.RotateVector3(CurrentStepDisplacementVector);

	// for( unsigned int j=0; j<3; j++)
	//   {
	//     CurrentStepDisplacement[j]  = CurrentStepDisplacementVector[j];
	//   }
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
    //Predicting first time Derivative
    //*********************************************************************************


    inline void PredictVelocity(ModelPart::NodeType& Node,
				array_1d<double, 3 > & CurrentVelocity,
				const array_1d<double, 3 > & StepDisplacement,
				const array_1d<double, 3 > & CurrentAcceleration,
				const array_1d<double, 3 > & PreviousAcceleration,
				const array_1d<double, 3 > & PreviousVelocity)

    {

      // BELYTSCHKO:
      // if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
      // 	{

      // 	  CurrentVelocity[0]  =   ( mNewmark.c0 * StepDisplacement[0]
      // 	  			    - ( (mNewmark.gamma / mNewmark.beta) - 1.0  ) * PreviousVelocity[0]
      // 	  			    - ( mNewmark.deltatime * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2 ) ) * PreviousAcceleration[0]) * mNewmark.static_dynamic;
      // 	}

      // if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
      // 	{
      // 	  CurrentVelocity[1]  =  (  mNewmark.c0 * StepDisplacement[1]
      // 	  			    - ( (mNewmark.gamma / mNewmark.beta) - 1.0  ) * PreviousVelocity[1]
      // 	  			    - ( mNewmark.deltatime * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2 ) ) * PreviousAcceleration[1]) * mNewmark.static_dynamic;
      // 	}

      // if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
      // 	{
      // 	  CurrentVelocity[2]  =  (  mNewmark.c0 * StepDisplacement[2]
      // 	  			   - ( (mNewmark.gamma / mNewmark.beta) - 1.0  ) * PreviousVelocity[2]
      // 	  			   - ( mNewmark.deltatime * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2 ) ) * PreviousAcceleration[2]) * mNewmark.static_dynamic;
      // 	}


      //SIMO:
      if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
      	{
      	  CurrentVelocity[0]  =  PreviousVelocity[0] + (mNewmark.c2 * PreviousAcceleration[0]
      			           + mNewmark.c3 * CurrentAcceleration[0] ) * mNewmark.static_dynamic;
      	}
      if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
      	{
      	  CurrentVelocity[1]  =  PreviousVelocity[1] + (mNewmark.c2 * PreviousAcceleration[1]
      			           + mNewmark.c3 * CurrentAcceleration[1] ) * mNewmark.static_dynamic;
      	}
      if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
      	{
      	  CurrentVelocity[2]  =  PreviousVelocity[2] + (mNewmark.c2 * PreviousAcceleration[2]
      			           + mNewmark.c3 * CurrentAcceleration[2] ) * mNewmark.static_dynamic;
      	}

    }
    //*********************************************************************************
    //Predicting second time Derivative
    //*********************************************************************************

    inline void PredictAcceleration(ModelPart::NodeType& Node,
				    array_1d<double, 3 > & CurrentAcceleration,
				    const array_1d<double, 3 > & StepDisplacement,
				    const array_1d<double, 3 > & PreviousAcceleration,
				    const array_1d<double, 3 > & PreviousVelocity)
    {
      //if(Node.Is(SLAVE)) return;

      // BELYTSCHKO::
      if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
      	{
      	  CurrentAcceleration[0] =  ( mNewmark.c1 * StepDisplacement[0] - mNewmark.c6 * PreviousVelocity[0] - mNewmark.c7 * PreviousAcceleration[0]) * mNewmark.static_dynamic;
      	}

      if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
      	{
      	  CurrentAcceleration[1] =  ( mNewmark.c1 * StepDisplacement[1] - mNewmark.c6 * PreviousVelocity[1] - mNewmark.c7 * PreviousAcceleration[1]) * mNewmark.static_dynamic;
      	}

      if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
      	{
      	  CurrentAcceleration[2] =  ( mNewmark.c1 * StepDisplacement[2] - mNewmark.c6 * PreviousVelocity[2] - mNewmark.c7 * PreviousAcceleration[2]) * mNewmark.static_dynamic;
      	}

      // SIMO:
      // if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
      // 	{
      // 	  //CurrentAcceleration[0] =  (- mNewmark.c4 * PreviousVelocity[0] + mNewmark.c5 * PreviousAcceleration[0]) * mNewmark.static_dynamic;

      // 	  //CurrentAcceleration[0] = PreviousAcceleration[0]; //not stable
      // 	  CurrentAcceleration[0] = 0; //ok
      // 	}

      // if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
      // 	{
      // 	  //CurrentAcceleration[1] =  (- mNewmark.c4 * PreviousVelocity[1] + mNewmark.c5 * PreviousAcceleration[1]) * mNewmark.static_dynamic;

      // 	  //CurrentAcceleration[1] = PreviousAcceleration[1]; //not stable
      // 	  CurrentAcceleration[1] = 0; //ok
      // 	}

      // if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
      // 	{
      // 	  //CurrentAcceleration[2] =  (- mNewmark.c4 * PreviousVelocity[2] + mNewmark.c5 * PreviousAcceleration[2]) * mNewmark.static_dynamic;

      // 	  //CurrentAcceleration[2] = PreviousAcceleration[2]; //not stable
      // 	  CurrentAcceleration[2] = 0; //ok
      // 	}

    }

    //*********************************************************************************
    //Predicting step rotation variable
    //*********************************************************************************

    inline void PredictStepRotation(ModelPart::NodeType& Node,
				    array_1d<double, 3 > & CurrentStepRotation,
				    const array_1d<double, 3 > & PreviousVelocity,
				    const array_1d<double, 3 > & PreviousAcceleration,
				    const array_1d<double, 3 > & CurrentAcceleration)
    {

      //CurrentStepRotation.clear(); //not needed :rigid bodies remain without rotation
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

      //noalias(CurrentStepRotation) = ( mNewmark.c4 * PreviousVelocity + mNewmark.c5 * PreviousAcceleration + CurrentAcceleration) * mNewmark.beta *  mNewmark.static_dynamic;

      //noalias(CurrentStepRotation) = ( mNewmark.deltatime * PreviousVelocity ) *  mNewmark.static_dynamic;


    }



    //*********************************************************************************
    //Predicting rotation variable
    //*********************************************************************************

    inline void PredictRotation(array_1d<double, 3 > & CurrentRotation,
				const array_1d<double, 3 > & StepRotation)
    {

        noalias(CurrentRotation) += StepRotation; //needed for imposed rotations

	// QuaternionType StepRotationQuaternion;
	// QuaternionType CurrentRotationQuaternion;
	// QuaternionType ReferenceRotationQuaternion;

	// Vector PreviousRotationVector  = ZeroVector(3);
	// Vector StepRotationVector = ZeroVector(3);

	// for( unsigned int j=0; j<3; j++)
	//   {
	//     PreviousRotationVector[j]  = PreviousRotation[j];
	//     StepRotationVector[j]      = StepRotation[j];
	//   }

	// // updated total rotations:
	// StepRotationQuaternion      = QuaternionType::FromRotationVector( StepRotationVector );
	// ReferenceRotationQuaternion = QuaternionType::FromRotationVector( PreviousRotationVector );

	// CurrentRotationQuaternion = StepRotationQuaternion * ReferenceRotationQuaternion;
	// CurrentRotationQuaternion.ToRotationVector( PreviousRotationVector );

	// for( unsigned int j=0; j<3; j++)
	//   {
	//     CurrentRotation[j] = PreviousRotationVector[j];
	//   }

    }

    //*********************************************************************************
    //Predicting first time Derivative
    //*********************************************************************************

    inline void PredictAngularVelocity(ModelPart::NodeType& Node,
				       array_1d<double, 3 > & CurrentAngularVelocity,
				       const array_1d<double, 3 > & StepRotation,
				       const array_1d<double, 3 > & CurrentAngularAcceleration,
				       const array_1d<double, 3 > & PreviousAngularAcceleration,
				       const array_1d<double, 3 > & PreviousAngularVelocity)
    {

      if(Node.Is(SLAVE)) return;

      // BELYTSCHKO:
      // if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
      // 	{

      // 	  CurrentAngularVelocity[0]  =   ( mNewmark.c0 * StepRotation[0]
      // 	  			    - ( (mNewmark.gamma / mNewmark.beta) - 1.0  ) * PreviousAngularVelocity[0]
      // 	  			    - ( mNewmark.deltatime * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2 ) ) * PreviousAngularAcceleration[0]) * mNewmark.static_dynamic;
      // 	}

      // if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularVelocity[1]  =  (  mNewmark.c0 * StepRotation[1]
      // 	  			    - ( (mNewmark.gamma / mNewmark.beta) - 1.0  ) * PreviousAngularVelocity[1]
      // 	  			    - ( mNewmark.deltatime * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2 ) ) * PreviousAngularAcceleration[1]) * mNewmark.static_dynamic;
      // 	}

      // if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularVelocity[2]  =  (  mNewmark.c0 * StepRotation[2]
      // 	  			   - ( (mNewmark.gamma / mNewmark.beta) - 1.0  ) * PreviousAngularVelocity[2]
      // 	  			   - ( mNewmark.deltatime * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2 ) ) * PreviousAngularAcceleration[2]) * mNewmark.static_dynamic;
      // 	}


      // SIMO:
      if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
	{
	  CurrentAngularVelocity[0]  =  PreviousAngularVelocity[0] + (mNewmark.c2 * PreviousAngularAcceleration[0]
						   + mNewmark.c3 * CurrentAngularAcceleration[0] ) * mNewmark.static_dynamic;

	}
      if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
	{
	  CurrentAngularVelocity[1]  =  PreviousAngularVelocity[1] + (mNewmark.c2 * PreviousAngularAcceleration[1]
						   + mNewmark.c3 * CurrentAngularAcceleration[1] ) * mNewmark.static_dynamic;

	}
      if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
	{
	  CurrentAngularVelocity[2]  =  PreviousAngularVelocity[2] + (mNewmark.c2 * PreviousAngularAcceleration[2]
							+ mNewmark.c3 * CurrentAngularAcceleration[2] ) * mNewmark.static_dynamic;

	}
    }

    //*********************************************************************************
    //Predicting second time Derivative
    //*********************************************************************************

    inline void PredictAngularAcceleration(ModelPart::NodeType& Node,
					   array_1d<double, 3 > & CurrentAngularAcceleration,
					   const array_1d<double, 3 > & StepRotation,
					   const array_1d<double, 3 > & PreviousAngularAcceleration,
					   const array_1d<double, 3 > & PreviousAngularVelocity)
    {
      if(Node.Is(SLAVE)) return;

      // BELYTSCHKO:
      if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
      	{
      	  CurrentAngularAcceleration[0] =  ( mNewmark.c1 * StepRotation[0] - mNewmark.c6 * PreviousAngularVelocity[0] - mNewmark.c7 * PreviousAngularAcceleration[0]) * mNewmark.static_dynamic;
      	}

       if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
      	{
      	  CurrentAngularAcceleration[1] =  ( mNewmark.c1 * StepRotation[1] - mNewmark.c6 * PreviousAngularVelocity[1] - mNewmark.c7 * PreviousAngularAcceleration[1]) * mNewmark.static_dynamic;

      	}

       if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
      	{
      	  CurrentAngularAcceleration[2] =  ( mNewmark.c1 * StepRotation[2] - mNewmark.c6 * PreviousAngularVelocity[2] - mNewmark.c7 * PreviousAngularAcceleration[2]) * mNewmark.static_dynamic;
      	}


      //SIMO
      // if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
      // 	{
      // 	  //CurrentAngularAcceleration[0] =  (- mNewmark.c4 * PreviousAngularVelocity[0] + mNewmark.c5 * PreviousAngularAcceleration[0]) * mNewmark.static_dynamic;

      // 	  //CurrentAngularAcceleration[0] = PreviousAngularAcceleration[0]; //not stable
      // 	  CurrentAngularAcceleration[0] = 0; //ok
      // 	}

      // if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
      // 	{
      // 	  //CurrentAngularAcceleration[1] =  (- mNewmark.c4 * PreviousAngularVelocity[1] + mNewmark.c5 * PreviousAngularAcceleration[1]) * mNewmark.static_dynamic;

      // 	  //CurrentAngularAcceleration[1] = PreviousAngularAcceleration[1]; //not stable
      // 	  CurrentAngularAcceleration[1] = 0; //ok
      // 	}

      // if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
      // 	{
      // 	  //CurrentAngularAcceleration[2] =  (- mNewmark.c4 * PreviousAngularVelocity[2] + mNewmark.c5 * PreviousAngularAcceleration[2]) * mNewmark.static_dynamic;

      // 	  //CurrentAngularAcceleration[2] = PreviousAngularAcceleration[2]; //not stable
      // 	  CurrentAngularAcceleration[2] = 0; //ok
      // 	}


    }


    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

    inline void UpdateVelocity(ModelPart::NodeType& Node,
			       array_1d<double, 3 > & CurrentVelocity,
			       const array_1d<double, 3 > & DeltaDisplacement,
			       const array_1d<double, 3 > & CurrentAcceleration,
			       const array_1d<double, 3 > & PreviousAcceleration,
			       const array_1d<double, 3 > & PreviousVelocity)
    {

      if(Node.Is(SLAVE)) return;

      // BELYTSCHKO:
      if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
	{
	  CurrentVelocity[0] =  PreviousVelocity[0] + mNewmark.c2 * PreviousAcceleration[0]
                                + mNewmark.c3 * CurrentAcceleration[0];
	}

      if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
	{
	  CurrentVelocity[1] =  PreviousVelocity[1] +  mNewmark.c2 * PreviousAcceleration[1]
                                + mNewmark.c3 * CurrentAcceleration[1];
	}

      if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
	{
	  CurrentVelocity[2] =  PreviousVelocity[2] + mNewmark.c2 * PreviousAcceleration[2]
                                + mNewmark.c3 * CurrentAcceleration[2];
      }

      // SIMO:
      // if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
      // 	{
      // 	  CurrentVelocity[0] += mNewmark.c0 * DeltaDisplacement[0];
      // 	}

      // if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
      // 	{
      // 	  CurrentVelocity[1] += mNewmark.c0 * DeltaDisplacement[1];
      // 	}

      // if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
      // 	{
      // 	  CurrentVelocity[2] += mNewmark.c0 * DeltaDisplacement[2];
      // }


    }

    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************


    inline void UpdateAcceleration(ModelPart::NodeType& Node,
				   array_1d<double, 3 > & CurrentAcceleration,
				   const array_1d<double, 3 > & DeltaDisplacement,
				   const array_1d<double, 3 > & StepDisplacement,
				   const array_1d<double, 3 > & PreviousAcceleration,
				   const array_1d<double, 3 > & PreviousVelocity)
    {

      if(Node.Is(SLAVE)) return;

      // BELYTSCHKO:
      if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
	{
	  CurrentAcceleration[0] =  ( mNewmark.c1 * StepDisplacement[0]- mNewmark.c6 * PreviousVelocity[0] - mNewmark.c7 * PreviousAcceleration[0]) * mNewmark.static_dynamic;
	}

      if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
	{
	  CurrentAcceleration[1] =  ( mNewmark.c1 * StepDisplacement[1]- mNewmark.c6 * PreviousVelocity[1] - mNewmark.c7 * PreviousAcceleration[1]) * mNewmark.static_dynamic;
	}

      if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
	{
	  CurrentAcceleration[2] =  ( mNewmark.c1 * StepDisplacement[2]- mNewmark.c6 * PreviousVelocity[2] - mNewmark.c7 * PreviousAcceleration[2]) * mNewmark.static_dynamic;
	}

      // SIMO:
      // if ((Node.pGetDof(DISPLACEMENT_X))->IsFixed() == false)
      // 	{
      // 	  CurrentAcceleration[0] +=  mNewmark.c1 * DeltaDisplacement[0];
      // 	}

      // if ((Node.pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
      // 	{
      // 	  CurrentAcceleration[1] +=  mNewmark.c1 * DeltaDisplacement[1];
      // 	}

      // if ((Node.pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
      // 	{
      // 	  CurrentAcceleration[2] +=  mNewmark.c1 * DeltaDisplacement[2];
      // }

    }




    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************


    inline void UpdateAngularVelocity(ModelPart::NodeType& Node,
                                      array_1d<double, 3 > & CurrentAngularVelocity,
				      const array_1d<double, 3 > & DeltaRotation,
				      const array_1d<double, 3 > & CurrentAngularAcceleration,
                                      const array_1d<double, 3 > & PreviousAngularAcceleration,
                                      const array_1d<double, 3 > & PreviousAngularVelocity)
    {

      if(Node.Is(SLAVE)) return;

      // BELYTSCHKO:
      if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
	{
	  CurrentAngularVelocity[0] = ( PreviousAngularVelocity[0] + (1-mNewmark.gamma) * mNewmark.deltatime * PreviousAngularAcceleration[0]
                                + mNewmark.gamma * mNewmark.deltatime * CurrentAngularAcceleration[0]) * mNewmark.static_dynamic;
	}

      if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
	{
	  CurrentAngularVelocity[1] = ( PreviousAngularVelocity[1] + (1-mNewmark.gamma) * mNewmark.deltatime * PreviousAngularAcceleration[1]
                                + mNewmark.gamma * mNewmark.deltatime * CurrentAngularAcceleration[1]) * mNewmark.static_dynamic;
	}

      if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
	{
	  CurrentAngularVelocity[2] = ( PreviousAngularVelocity[2] + (1-mNewmark.gamma) * mNewmark.deltatime * PreviousAngularAcceleration[2]
                                + mNewmark.gamma * mNewmark.deltatime * CurrentAngularAcceleration[2]) * mNewmark.static_dynamic;
      }

      // SIMO:
      // if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularVelocity[0] +=  mNewmark.c0 * DeltaRotation[0];
      // 	}

      // if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularVelocity[1] +=  mNewmark.c0 * DeltaRotation[1];
      // 	}

      // if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularVelocity[2] +=  mNewmark.c0 * DeltaRotation[2];
      // }
    }


    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************


    inline void UpdateAngularAcceleration(ModelPart::NodeType& Node,
                                          array_1d<double, 3 > & CurrentAngularAcceleration,
					  const array_1d<double, 3 > & DeltaRotation,
					  const array_1d<double, 3 > & StepRotation,
                                          const array_1d<double, 3 > & PreviousAngularAcceleration,
                                          const array_1d<double, 3 > & PreviousAngularVelocity)
    {
      if(Node.Is(SLAVE)) return;

      //BELYTSCHKO:
      if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
      	{
      	  CurrentAngularAcceleration[0] =  ( mNewmark.c1 * StepRotation[0]- mNewmark.c6 * PreviousAngularVelocity[0] - mNewmark.c7 * PreviousAngularAcceleration[0]) * mNewmark.static_dynamic;
      	}

      if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
      	{
      	  CurrentAngularAcceleration[1] =  ( mNewmark.c1 * StepRotation[1]- mNewmark.c6 * PreviousAngularVelocity[1] - mNewmark.c7 * PreviousAngularAcceleration[1]) * mNewmark.static_dynamic;
      	}

      if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
      	{
      	  CurrentAngularAcceleration[2] =  ( mNewmark.c1 * StepRotation[2]- mNewmark.c6 * PreviousAngularVelocity[2] - mNewmark.c7 * PreviousAngularAcceleration[2]) * mNewmark.static_dynamic;
      	}

      // SIMO:
      // if ((Node.pGetDof(ROTATION_X))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularAcceleration[0] +=  mNewmark.c1 * DeltaRotation[0];
      // 	}

      // if ((Node.pGetDof(ROTATION_Y))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularAcceleration[1] +=  mNewmark.c1 * DeltaRotation[1];
      // 	}

      // if ((Node.pGetDof(ROTATION_Z))->IsFixed() == false)
      // 	{
      // 	  CurrentAngularAcceleration[2] +=  mNewmark.c1 * DeltaRotation[2];
      // }
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
            noalias(LHS_Contribution) += M * mNewmark.static_dynamic;

            //std::cout<<" Mass Matrix "<<M<<" coeficient "<<mNewmark.c0<<std::endl;
        }

        //adding  damping contribution
        if (D.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += D * mNewmark.static_dynamic;

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
            noalias(RHS_Contribution) -=  mNewmark.static_dynamic * fm;
        }

        //adding damping contribution
        if (fd.size() != 0)
        {
	    noalias(RHS_Contribution) -=  mNewmark.static_dynamic * fd;
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

            (mVector.a[thread]) *= mNewmark.static_dynamic ;

            rCurrentElement->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mVector.ap[thread] * mNewmark.static_dynamic;

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
            //KRATOS_WATCH( prod(M, macc[thread] ) )

        }

        //adding damping contribution
        if (D.size1() != 0)
        {
            rCurrentElement->GetFirstDerivativesVector(mVector.v[thread], 0);

            (mVector.v[thread]) *= mNewmark.static_dynamic ;

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

            (mVector.a[thread]) *=  mNewmark.static_dynamic;

            rCurrentCondition->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) +=  mVector.ap[thread] * mNewmark.static_dynamic;

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

        //adding damping contribution
        //damping contribution
        if (D.size1() != 0)
        {
            rCurrentCondition->GetFirstDerivativesVector(mVector.v[thread], 0);

            (mVector.v[thread]) *= mNewmark.static_dynamic ;

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

#endif /* KRATOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_ROTATION_SCHEME  defined */


