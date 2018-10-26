//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_EXPLICIT_HAMILTON_SCHEME_H_INCLUDED)
#define  KRATOS_EXPLICIT_HAMILTON_SCHEME_H_INCLUDED

/* System includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "utilities/quaternion.h"

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

  template<class TSparseSpace,
	   class TDenseSpace //= DenseSpace<double>
	   >
  class ExplicitHamiltonScheme : public Scheme<TSparseSpace,TDenseSpace>
  {

  public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( ExplicitHamiltonScheme );

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType ElementsArrayType;

    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    typedef ModelPart::NodesContainerType NodesArrayType;

    typedef Quaternion<double> QuaternionType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ExplicitHamiltonScheme(const double  rMaximumDeltaTime,
			   const double  rDeltaTimeFraction,
			   const double  rDeltaTimePredictionLevel,
			   const bool    rRayleighDamping)
      : Scheme<TSparseSpace,TDenseSpace>()
    {

      mDeltaTime.PredictionLevel  = rDeltaTimePredictionLevel;

      mDeltaTime.Maximum          = rMaximumDeltaTime;
      mDeltaTime.Fraction         = rDeltaTimeFraction;

      mRayleighDamping            = rRayleighDamping;

      //Allocate auxiliary memory
      int NumThreads = OpenMPUtils::GetNumThreads();


      mMatrix.D.resize(NumThreads);
      //mMatrix.M.resize(NumThreads);

      mVector.v.resize(NumThreads);
      //mVector.a.resize(NumThreads);
      //mVector.ap.resize(NumThreads);


    }


    /** Destructor.
     */
    virtual ~ExplicitHamiltonScheme() {}


    /*@} */
    /**@name Operators
     */
    /*@{ */

    virtual int Check(ModelPart& r_model_part)
    {
      KRATOS_TRY

        BaseType::Check(r_model_part);

      if(r_model_part.GetBufferSize() < 2)
        {
	  KRATOS_THROW_ERROR(std::logic_error, "Insufficient buffer size for Central Difference Scheme. It has to be 2", "")
	    }

      KRATOS_CATCH("")

      return 0;
    }


    virtual void Initialize(ModelPart& r_model_part)
    {
      KRATOS_TRY

        mModelPart = r_model_part;

      if(mDeltaTime.PredictionLevel>0)
        {
	  this->CalculateDeltaTime(r_model_part);
        }

      //initialize scheme variables
      this->InitializeExplicitScheme(r_model_part);

      ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();
      rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] = true;

      mTime.Previous       = 0.0;
      mTime.PreviousMiddle = 0.0;

      mSchemeIsInitialized = true;

      KRATOS_CATCH("")
	}

    //***************************************************************************

    void InitializeSolutionStep(ModelPart& r_model_part,
				TSystemMatrixType& A,
				TSystemVectorType& Dx,
				TSystemVectorType& b
				)
    {
      KRATOS_TRY

      BaseType::InitializeSolutionStep(r_model_part,A,Dx,b);

      if(mDeltaTime.PredictionLevel>1)
	{
	  this->CalculateDeltaTime(r_model_part);
	}


      //initialize residual
      //this->InitializeResidual(r_model_part);

      //initialize movements
      this->InitializeMovements(r_model_part);

      //initialize update flags
      mUpdatePositionFlag = true;
      mUpdateRotationFlag = true;
      mUpdateMomentumFlag = true;


      KRATOS_CATCH("")
    }

    //**************************************************************************

    void InitializeResidual( ModelPart& r_model_part )

    {
        KRATOS_TRY

        NodesArrayType& pNodes   = r_model_part.Nodes();

        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

        #pragma omp parallel for
        for(int k=0; k<number_of_threads; k++)
        {
            typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
            typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

            for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {
	      (i)->FastGetSolutionStepValue(FORCE_RESIDUAL).clear();
	      (i)->FastGetSolutionStepValue(MOMENT_RESIDUAL).clear();
            }
        }

        KRATOS_CATCH("")
    }


    //**************************************************************************

    void InitializeMovements( ModelPart& r_model_part )

    {
        KRATOS_TRY

        NodesArrayType& pNodes   = r_model_part.Nodes();

        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

        #pragma omp parallel for
        for(int k=0; k<number_of_threads; k++)
        {
            typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
            typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

            for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {
	      array_1d<double, 3 >& CurrentStepRotation = (i)->FastGetSolutionStepValue(STEP_ROTATION);
	      if( (i->pGetDof(ROTATION_X))->IsFree() )
		CurrentStepRotation[0]  = 0;
	      if( (i->pGetDof(ROTATION_Y))->IsFree() )
		CurrentStepRotation[1]  = 0;
	      if( (i->pGetDof(ROTATION_Z))->IsFree() )
		CurrentStepRotation[2]  = 0;
	    }
        }

        KRATOS_CATCH("")
    }


    //***************************************************************************

    void FinalizeSolutionStep(ModelPart& rModelPart,
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

            for (ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; ++itElem)
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

            for (ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; ++itCond)
            {
                itCond->FinalizeSolutionStep(CurrentProcessInfo);
            }
        }
        KRATOS_CATCH( "" )
    }


    //***************************************************************************

    void CalculateDeltaTime(ModelPart& r_model_part)
    {

      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
      ElementsArrayType& pElements      = r_model_part.Elements();

#ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
#else
      int number_of_threads = 1;

#endif

      vector<unsigned int> element_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

      double safety_factor = 0.65;  //most autors recommend a value near 0.80 (Belytschko - Nonlinear FE.. 2000. chap 6. pag. 315)

      Vector delta_times(number_of_threads);

      double stable_delta_time = 0.00;

      for(int i = 0; i < number_of_threads; i++)
	delta_times[i] = mDeltaTime.Maximum/safety_factor;

#pragma omp parallel for private(stable_delta_time)
      for(int k=0; k<number_of_threads; k++)
        {
	  typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
	  typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
	  for(ElementsArrayType::iterator it=it_begin; it!= it_end; ++it)
            {

	      //get geometric and material properties
	      double length   = it->GetGeometry().Length();
	      double alpha    = it->GetProperties()[RAYLEIGH_ALPHA];
	      double beta     = it->GetProperties()[RAYLEIGH_BETA];
	      double E        = it->GetProperties()[YOUNG_MODULUS];
	      double v        = it->GetProperties()[POISSON_RATIO];
	      double ro       = it->GetProperties()[DENSITY];

	      //compute courant criterion
	      double bulk       = E/(3.0*(1.0-2.0*v));
	      double wavespeed  = sqrt(bulk/ro);
	      double w          = 2.0*wavespeed/length;   //frequency

	      double psi        = 0.5*(alpha/w + beta*w); //critical ratio;

	      stable_delta_time = (2.0/w)*(sqrt(1.0 + psi*psi)-psi);

	      if(stable_delta_time > 0.00){
		if(stable_delta_time < delta_times[k]){
		  delta_times[k] = stable_delta_time;
		}
	      }

            }
        }

      stable_delta_time  = *std::min_element(delta_times.begin(), delta_times.end());

      stable_delta_time *= safety_factor * 0.5; //extra factor added to get an stable delta time

      if(stable_delta_time < mDeltaTime.Maximum){

	rCurrentProcessInfo[DELTA_TIME] = stable_delta_time;

	std::cout<< "  New computed stable Time step = "<< stable_delta_time <<"         "<< std::endl;
	std::cout<< "   " << std::endl;
      }

      KRATOS_CATCH("")
    }

    //***************************************************************************

    void InitializeExplicitScheme(ModelPart& r_model_part)
    {
      KRATOS_TRY

        //Set velocity to zero...

        NodesArrayType& pNodes        = r_model_part.Nodes();

#ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
#else
      int number_of_threads = 1;
#endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

#pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
	{
	  typename NodesArrayType::iterator i_begin = pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end   = pNodes.ptr_begin()+node_partition[k+1];

	  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
	    {
	      i->FastGetSolutionStepValue(POSITION_MOMENTUM).clear();
	      i->FastGetSolutionStepValue(ROTATION_MOMENTUM).clear();

	      i->FastGetSolutionStepValue(FORCE_RESIDUAL).clear();
	      i->FastGetSolutionStepValue(MOMENT_RESIDUAL).clear();

	      i->FastGetSolutionStepValue(DISPLACEMENT).clear();
	      i->FastGetSolutionStepValue(ROTATION).clear();
	    }
	}

      KRATOS_CATCH("")
    }

    /**
       Performing the update of the solution.
    */
    //***************************************************************************
    virtual void Update(ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			)
    {
      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();

      mUpdatePositionFlag   = rCurrentProcessInfo[POSITION_UPDATE_LABEL];
      mUpdateRotationFlag   = rCurrentProcessInfo[ROTATION_UPDATE_LABEL];
      mUpdateMomentumFlag   = rCurrentProcessInfo[MOMENTUM_UPDATE_LABEL];


      if( mUpdatePositionFlag == true ){
	//1) Update nodal positions:
	UpdateNodalPosition( r_model_part );

	mUpdatePositionFlag = false;
      }

      if( mUpdateRotationFlag == true ){
	//1) Update nodal positions:
	UpdateNodalRotation( r_model_part );

	mUpdateRotationFlag = false;
      }

      if( mUpdateMomentumFlag == true ){

	//2) Update momentum equations:
	UpdateNodalMomentum( r_model_part );

	mUpdateMomentumFlag = false;
      }


      KRATOS_CATCH("")
    }


    //***************************************************************************
    //***************************************************************************


    void UpdateNodalPosition(ModelPart& r_model_part)
    {
      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes            = r_model_part.Nodes();

      //Step Update
      mTime.Delta     = rCurrentProcessInfo[DELTA_TIME];

      double  Alpha   = rCurrentProcessInfo[ALPHA_TRAPEZOIDAL_RULE];

      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

      #pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
        {
	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {

	      //is active by default
	      bool node_is_active = true;
	      if( (i)->IsDefined(ACTIVE) )
		node_is_active = (i)->Is(ACTIVE);

	      if(node_is_active)
		{
		  // Update DISPLACEMENTS, LINEAR VELOCITIES AND ACCELERATIONS

		  const double& nodal_mass                    = i->FastGetSolutionStepValue(NODAL_MASS);

		  array_1d<double,3>& force_residual          = i->FastGetSolutionStepValue(FORCE_RESIDUAL);
		  array_1d<double,3>& position_momentum       = i->FastGetSolutionStepValue(POSITION_MOMENTUM);


		  array_1d<double,3>& CurrentDisplacement    = i->FastGetSolutionStepValue(DISPLACEMENT,0);
		  array_1d<double,3>& ReferenceDisplacement  = i->FastGetSolutionStepValue(DISPLACEMENT,1);

		  ReferenceDisplacement = CurrentDisplacement;

		  //Solution of the explicit equation: force_residual is (Fext-Fint) then the sign is changed respect the reference formulation

		  //std::cout<<" nodal mass: ("<<i->Id()<<"): "<<nodal_mass<<std::endl;

		  if( (i->pGetDof(DISPLACEMENT_X))->IsFree() ){

		    CurrentDisplacement[0]  = ReferenceDisplacement[0] + ( mTime.Delta / nodal_mass ) * ( position_momentum[0] + Alpha * mTime.Delta * force_residual[0] );

		  }


		  if( (i->pGetDof(DISPLACEMENT_Y))->IsFree() ){

		    CurrentDisplacement[1]  = ReferenceDisplacement[1] + ( mTime.Delta / nodal_mass ) * ( position_momentum[1] + Alpha * mTime.Delta * force_residual[1] );

		  }


		  if( i->HasDofFor(DISPLACEMENT_Z) ){

		    if( (i->pGetDof(DISPLACEMENT_Z))->IsFree() ){

		      CurrentDisplacement[2]  = ReferenceDisplacement[2] + ( mTime.Delta / nodal_mass ) * ( position_momentum[2] + Alpha * mTime.Delta * force_residual[2] );

		    }
		  }

		  //std::cout<<" Node ["<<i->Id()<<"] (mass:"<<nodal_mass<<",force:"<<force_residual<<",moment: "<<moment_residual<<", position: "<<position_momentum<<"); Alpha:"<<Alpha<<std::endl;
		  //std::cout<<" Node ["<<i->Id()<<"] (Dt/mass:"<<mTime.Delta / nodal_mass<<", position_momentum:"<<position_momentum<<", alpha*Dt*force_residual: "<<Alpha * mTime.Delta * force_residual<<"); Disp:"<<CurrentDisplacement<<", "<<ReferenceDisplacement<<std::endl;

		  //linear velocities and accelerations
		  array_1d<double, 3 >  DeltaDisplacement = CurrentDisplacement-ReferenceDisplacement;

		  array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY, 0);
		  array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

		  array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
		  //array_1d<double, 3 > & PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

		  PreviousVelocity = CurrentVelocity;

		  UpdateVelocity (CurrentVelocity, DeltaDisplacement, mTime.Delta);

		  array_1d<double, 3 >  DeltaVelocity = CurrentVelocity-PreviousVelocity;

		  UpdateAcceleration (CurrentAcceleration, DeltaVelocity, mTime.Delta);

		  //---------------------------//

		}
            }
        }


      KRATOS_CATCH("")

    }


    //***************************************************************************
    //***************************************************************************


    void UpdateNodalRotation(ModelPart& r_model_part)
    {
      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes            = r_model_part.Nodes();

      //Step Update
      mTime.Delta     = rCurrentProcessInfo[DELTA_TIME];

      //double  Alpha   = rCurrentProcessInfo[ALPHA_TRAPEZOIDAL_RULE];


      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

      #pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
        {
	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {

	      //is active by default
	      bool node_is_active = true;
	      if( (i)->IsDefined(ACTIVE) )
		node_is_active = (i)->Is(ACTIVE);

	      if(node_is_active)
		{

		  // Update ROTATIONS, ANGULAR VELOCITIES AND ACCELERATIONS

		  //Total Rotation
		  array_1d<double, 3 > & CurrentRotation     = (i)->FastGetSolutionStepValue(ROTATION, 0);
		  array_1d<double, 3 > & ReferenceRotation   = (i)->FastGetSolutionStepValue(ROTATION, 1);

		  //StepRotation
		  array_1d<double, 3 > & CurrentStepRotation = (i)->FastGetSolutionStepValue(STEP_ROTATION, 0);

		  Vector TotalRotationVector = ZeroVector(3);
		  Vector StepRotationVector  = ZeroVector(3);

		  QuaternionType StepRotationQuaternion;
		  QuaternionType CurrentRotationQuaternion;
		  QuaternionType ReferenceRotationQuaternion;

		  for( unsigned int j=0; j<3; j++)
		    {

		      TotalRotationVector[j] = ReferenceRotation[j];    //previous iteration total rotation
		      StepRotationVector[j]  = CurrentStepRotation[j];  //previous iteration step rotation
		    }

		  StepRotationQuaternion = QuaternionType::FromRotationVector( StepRotationVector );

		  // updated step rotations:
		  ReferenceRotationQuaternion = QuaternionType::FromRotationVector( TotalRotationVector );

		  CurrentRotationQuaternion = StepRotationQuaternion * ReferenceRotationQuaternion;

		  CurrentRotationQuaternion.ToRotationVector( TotalRotationVector );

		  for( unsigned int j=0; j<3; j++)
		    {
		      CurrentRotation[j] = TotalRotationVector[j];
		    }


		  //angular velocities and accelerations
		  array_1d<double, 3 > & CurrentAngularVelocity      = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY, 0);
		  array_1d<double, 3 > & PreviousAngularVelocity     = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY, 1);

		  array_1d<double, 3 > & CurrentAngularAcceleration  = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION, 0);

		  PreviousAngularVelocity = CurrentAngularVelocity;

		  UpdateAngularVelocity (CurrentAngularVelocity, StepRotationVector, mTime.Delta);

		  array_1d<double, 3 >  DeltaAngularVelocity = CurrentAngularVelocity-PreviousAngularVelocity;

		  UpdateAngularAcceleration (CurrentAngularAcceleration, DeltaAngularVelocity, mTime.Delta);

		  //------------------------

		}
            }
        }


      KRATOS_CATCH("")

    }


    //***************************************************************************
    //***************************************************************************


    void UpdateNodalMomentum(ModelPart& r_model_part)
    {
      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes            = r_model_part.Nodes();

      //Step Update
      mTime.Delta     = rCurrentProcessInfo[DELTA_TIME];

      double Alpha    = rCurrentProcessInfo[ALPHA_TRAPEZOIDAL_RULE];

#ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
#else
      int number_of_threads = 1;
#endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

#pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
        {
	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {
	      //is active by default
	      bool node_is_active = true;
	      if( (i)->IsDefined(ACTIVE) )
		node_is_active = (i)->Is(ACTIVE);

	      if(node_is_active)
		{
		  //Current step information "N+1" (before step update).
		  array_1d<double,3>& previous_force_residual  = i->FastGetSolutionStepValue(FORCE_RESIDUAL,1);
		  array_1d<double,3>& previous_moment_residual = i->FastGetSolutionStepValue(MOMENT_RESIDUAL,1);

		  array_1d<double,3>& force_residual           = i->FastGetSolutionStepValue(FORCE_RESIDUAL);
		  array_1d<double,3>& moment_residual          = i->FastGetSolutionStepValue(MOMENT_RESIDUAL);

		  array_1d<double,3>& previous_position_momentum   = i->FastGetSolutionStepValue(POSITION_MOMENTUM,1);
		  array_1d<double,3>& previous_rotation_momentum   = i->FastGetSolutionStepValue(ROTATION_MOMENTUM,1);

		  array_1d<double,3>& position_momentum        = i->FastGetSolutionStepValue(POSITION_MOMENTUM);
		  array_1d<double,3>& rotation_momentum        = i->FastGetSolutionStepValue(ROTATION_MOMENTUM);


		  //Solution of the explicit momentum equations: force_residual is (Fext-Fint) then the sign is changed respect the reference formulation (the same for moment_residual)
		  position_momentum = previous_position_momentum + mTime.Delta * ( (1.0 - Alpha) * previous_force_residual  + Alpha * force_residual );

		  rotation_momentum = previous_rotation_momentum + mTime.Delta * ( (1.0 - Alpha) * previous_moment_residual + Alpha * moment_residual );

		  //std::cout<<" Node ["<<i->Id()<<"] (force_residual:"<<force_residual<<", "<<previous_force_residual<<",moment_residual: "<<moment_residual<<", "<<previous_moment_residual<<"); position_momentum:"<<position_momentum<<std::endl;
		}
            }
        }


      KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************

    /** this function is designed to be called in the builder and solver to introduce*/

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

	(rCurrentElement) -> CalculateSecondDerivativesContributions(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

	//add explicit contribution of the Element Residual (RHS) to nodal Force Residual (nodal RHS)
	(rCurrentElement) -> AddExplicitContribution(LHS_Contribution, TANGENT_MATRIX, TANGENT_LYAPUNOV, rCurrentProcessInfo);

	//add explicit contribution of the Element Residual (RHS) to nodal Moment Residual (nodal RHS)
	(rCurrentElement) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, RESIDUAL_LYAPUNOV, rCurrentProcessInfo);

	//std::cout<<" Element ["<<(rCurrentElement)->Id()<<"]: (Tangent:"<<LHS_Contribution<<") "<<std::endl;
	//std::cout<<" Element ["<<(rCurrentElement)->Id()<<"]: (Residual:"<<RHS_Contribution<<") "<<std::endl;

	//----------------------

	// (rCurrentElement) -> CalculateFirstDerivativesContributions(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

	// //add explicit contribution of the Element Residual (RHS) to nodal Force Residual (nodal RHS)
	// (rCurrentElement) -> AddExplicitContribution(LHS_Contribution, TANGENT_MATRIX, TANGENT_LYAPUNOV, rCurrentProcessInfo);

	// //add explicit contribution of the Element Residual (RHS) to nodal Moment Residual (nodal RHS)
	// (rCurrentElement) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, RESIDUAL_LYAPUNOV, rCurrentProcessInfo);

	// //std::cout<<" Element ["<<(rCurrentElement)->Id()<<"]: (Tangent:"<<LHS_Contribution<<") "<<std::endl;
	// //std::cout<<" Element ["<<(rCurrentElement)->Id()<<"]: (Residual:"<<RHS_Contribution<<") "<<std::endl;

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
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

	// (rCurrentCondition) -> CalculateSecondDerivativesContributions(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

	// //add explicit contribution of the Element Residual (RHS) to nodal Force Residual (nodal RHS)
	// (rCurrentCondition) -> AddExplicitContribution(LHS_Contribution, TANGENT_MATRIX, TANGENT_LYAPUNOV, rCurrentProcessInfo);

	// //add explicit contribution of the Element Residual (RHS) to nodal Moment Residual (nodal RHS)
	// (rCurrentCondition) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, RESIDUAL_LYAPUNOV, rCurrentProcessInfo);

	// //std::cout<<" Condition ["<<(rCurrentCondition)->Id()<<"]: (Tangent:"<<LHS_Contribution<<") "<<std::endl;
	// //std::cout<<" Condition ["<<(rCurrentCondition)->Id()<<"]: (Residual:"<<RHS_Contribution<<") "<<std::endl;

	// //----------------------

	// (rCurrentCondition) -> CalculateFirstDerivativesContributions(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

	// //add explicit contribution of the Element Residual (RHS) to nodal Force Residual (nodal RHS)
	// (rCurrentCondition) -> AddExplicitContribution(LHS_Contribution, TANGENT_MATRIX, TANGENT_LYAPUNOV, rCurrentProcessInfo);

	// //add explicit contribution of the Element Residual (RHS) to nodal Moment Residual (nodal RHS)
	// (rCurrentCondition) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, RESIDUAL_LYAPUNOV, rCurrentProcessInfo);

	// //std::cout<<" Condition ["<<(rCurrentCondition)->Id()<<"]: (Tangent:"<<LHS_Contribution<<") "<<std::endl;
	// //std::cout<<" Condition ["<<(rCurrentCondition)->Id()<<"]: (Residual:"<<RHS_Contribution<<") "<<std::endl;

        KRATOS_CATCH( "" )
    }


    //***************************************************************************
    //***************************************************************************

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
				    LocalSystemVectorType& RHS_Contribution,
				    Element::EquationIdVectorType& EquationId,
				    ProcessInfo& rCurrentProcessInfo)
    {

      KRATOS_TRY


      //basic operations for the element considered
      (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);

      //std::cout<<" Element ["<<(rCurrentElement)->Id()<<"]: (Residual:"<<RHS_Contribution<<") "<<std::endl;

      //add explicit contribution of the Element Residual (RHS) to nodal Force Residual (nodal RHS)
      (rCurrentElement) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);

      //add explicit contribution of the Element Residual (RHS) to nodal Moment Residual (nodal RHS)
      (rCurrentElement) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, MOMENT_RESIDUAL, rCurrentProcessInfo);


      KRATOS_CATCH( "" )
	}


    //***************************************************************************
    //***************************************************************************

    /** functions totally analogous to the precedent but applied to
	the "condition" objects
    */
    virtual void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
						      LocalSystemVectorType& RHS_Contribution,
						      Element::EquationIdVectorType& EquationId,
						      ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      //basic operations for the element considered
      (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);

      //std::cout<<" Condition ["<<(rCurrentCondition)->Id()<<"]: (Residual:"<<RHS_Contribution<<") "<<std::endl;

      //add explicit contribution of the Condition Residual (RHS) to nodal Force Residual (nodal RHS)
      (rCurrentCondition) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);

      //add explicit contribution of the Condition Residual (RHS) to nodal Moment Residual (nodal RHS)
      (rCurrentCondition) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, MOMENT_RESIDUAL, rCurrentProcessInfo);


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


    /*@} */

  protected:
    /**@name Protected static Member Variables */
    /*@{ */


    struct GeneralMatrices
    {
      //std::vector< Matrix > M;     //first derivative matrix  (usually mass matrix)
      std::vector< Matrix > D;     //second derivative matrix (usually damping matrix)
    };

    struct GeneralVectors
    {
      std::vector< Vector > v;    //velocity
      //std::vector< Vector > a;    //acceleration
      //std::vector< Vector > ap;   //previous acceleration
    };


    struct DeltaTimeParameters
    {
      double PredictionLevel; // 0, 1, 2

      double Maximum;  //maximum delta time
      double Fraction; //fraction of the delta time
    };


    struct TimeVariables
    {
      double PreviousMiddle; //n-1/2
      double Previous;       //n
      double Middle;         //n+1/2
      double Current;        //n+1

      double Delta;          //time step
    };


    GeneralMatrices     mMatrix;
    GeneralVectors      mVector;

    bool                mSchemeIsInitialized;

    ModelPart           mModelPart;


    TimeVariables       mTime;
    DeltaTimeParameters mDeltaTime;
    bool                mRayleighDamping;

    bool                mUpdatePositionFlag;
    bool                mUpdateRotationFlag;
    bool                mUpdateMomentumFlag;


    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    /*@} */
    /**@name Protected Operations*/
    /*@{ */

    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

    inline void UpdateVelocity(array_1d<double, 3 >& rCurrentVelocity,
			       const array_1d<double, 3 >& rDeltaDisplacement,
			       const double& rDeltaTime)
    {
      noalias(rCurrentVelocity) = (1.0/rDeltaTime) * rDeltaDisplacement;
    }


    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************

    inline void UpdateAcceleration(array_1d<double, 3 >& rCurrentAcceleration,
                                   const array_1d<double, 3 >& rDeltaVelocity,
				   const double& rDeltaTime)
    {
      noalias(rCurrentAcceleration) = (1.0/rDeltaTime) *  rDeltaVelocity;
    }

    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

    inline void UpdateAngularVelocity(array_1d<double, 3 >& rCurrentVelocity,
				      const array_1d<double, 3 >& rDeltaRotation,
				      const double& rDeltaTime)
    {
      noalias(rCurrentVelocity) = (1.0/rDeltaTime) * rDeltaRotation;
    }


    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************

    inline void UpdateAngularAcceleration(array_1d<double, 3 >& rCurrentAcceleration,
					  const array_1d<double, 3 >& rDeltaVelocity,
					  const double& rDeltaTime)
    {
      noalias(rCurrentAcceleration) = (1.0/rDeltaTime) * rDeltaVelocity;
    }

    /*@} */
    /**@name Protected  Access */
    /*@{ */

    /*@} */
    /**@name Protected Inquiry */
    /*@{ */

    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */

    /*@} */

  private:
    /**@name Static Member Variables */
    /*@{ */

    /*@} */
    /**@name Member Variables */
    /*@{ */

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
    /**@name Un accessible methods */
    /*@{ */

    /*@} */

  }; /* Class Scheme */

  /*@} */

  /**@name Type Definitions */
  /*@{ */


  /*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_HAMILTON_SCHEME_H_INCLUDED  defined */

