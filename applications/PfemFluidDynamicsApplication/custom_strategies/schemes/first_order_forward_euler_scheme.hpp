//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:    MSantasuna JMCarbonell $
//   Last modified by:    $Co-Author:                AFranci $
//   Date:                $Date:                  April 2018 $
//   Revision:            $Revision:                     0.0 $
//  
//

#if !defined(KRATOS_FIRST_ORDER_FORWARD_EULER_SCHEME)
#define  KRATOS_FIRST_ORDER_FORWARD_EULER_SCHEME


// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"

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

  template<class TSparseSpace, class TDenseSpace>
  class FirstOrderForwardEulerScheme : public Scheme<TSparseSpace,TDenseSpace>
  {
  public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( FirstOrderForwardEulerScheme );

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

    
  protected:


    struct DeltaTimeParameters
    {
      double PredictionLevel; // 0, 1, 2
      double Maximum;         //maximum delta time
      double Fraction;        //fraction of the delta time
    }; 


    struct TimeVariables
    {
      double PreviousMiddle; //n-1/2
      double Previous;       //n
      double Middle;         //n+1/2
      double Current;        //n+1

      double Delta;          //time step
    };

  public:
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors    
    FirstOrderForwardEulerScheme(const double  rMaximumDeltaTime,
				 const double  rDeltaTimeFraction,
				 const double  rDeltaTimePredictionLevel,
				 const bool    rRayleighDamping = false)
      : Scheme<TSparseSpace,TDenseSpace>()
    {

      mDeltaTime.PredictionLevel  = rDeltaTimePredictionLevel;

      mDeltaTime.Maximum          = rMaximumDeltaTime;

      mDeltaTime.Fraction         = rDeltaTimeFraction;

      mRayleighDamping            = rRayleighDamping;
	
      // Allocate auxiliary memory
      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      mMatrix.resize(NumThreads);
      mVector.resize(NumThreads);

      mSchemeIsInitialized = false;
    }

    ///Copy constructor
    FirstOrderForwardEulerScheme(FirstOrderForwardEulerScheme& rOther)
    :BaseType(rOther)
    ,mMatrix(rOther.mMatrix)
    ,mVector(rOther.mVector)
    ,mDeltaTime(rOther.mDeltaTime)
    ,mRayleighDamping(rOther.mRayleighDamping)
    ,mSchemeIsInitialized(rOther.mSchemeIsInitialized)	 
    {
    }

    /// Clone
    BaseTypePointer Clone() override
    {
      return BaseTypePointer( new FirstOrderForwardEulerScheme(*this) );
    }

    /// Destructor
    virtual ~FirstOrderForwardEulerScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    virtual void Initialize(ModelPart& rModelPart) override
    {
      KRATOS_TRY

	if( (mDeltaTime.PredictionLevel>0) && (mSchemeIsInitialized==false) )
	  {
	    this->CalculateDeltaTime(rModelPart);
	  }

      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

      //Preparing the time values for the first step (where time = initial_time + dt)
      mTime.Current         = rCurrentProcessInfo[TIME]+rCurrentProcessInfo[DELTA_TIME];
      mTime.Delta           = rCurrentProcessInfo[DELTA_TIME];
      mTime.Middle          = mTime.Current -  0.5*mTime.Delta;
      mTime.Previous        = mTime.Current -  mTime.Delta;
      mTime.PreviousMiddle  = mTime.Current - 1.5*mTime.Delta;

      if(mSchemeIsInitialized==false)
	{
	  this->InitializeExplicitScheme(rModelPart);
	}
      else
	{
	  this->SchemeCustomInitialization(rModelPart);
	}

      mSchemeIsInitialized = true;

      KRATOS_CATCH("")
	}


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
      KRATOS_TRY

	BaseType::InitializeSolutionStep(rModelPart,A,Dx,b);

      // if(mDeltaTime.PredictionLevel>1)
      // 	{
      // 	  this->CalculateDeltaTime(rModelPart);
      // 	}
      this->CalculateDeltaTime(rModelPart);

      this->InitializeResidual(rModelPart);       
	
      KRATOS_CATCH("")
	}

    /**
     * Performing the update of the solution.
     * incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param rModelPart
     * @param rDofSet set of all primary variables
     * @param A	LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void Update(ModelPart& rModelPart,
		DofsArrayType& rDofSet,
		TSystemMatrixType& A,
		TSystemVectorType& Dx,
		TSystemVectorType& b) override
    {
      KRATOS_TRY

	// std::cout<<"Update in forward euler"<<std::endl;
      
	ProcessInfo& rCurrentProcessInfo  = rModelPart.GetProcessInfo();

      //Step Update
      mTime.Current   = rCurrentProcessInfo[TIME];  //the first step is (time = initial_time + delta time )
      mTime.Delta     = rCurrentProcessInfo[DELTA_TIME];
      mTime.Middle    = 0.5 * ( mTime.Previous + mTime.Current );

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);
	
      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {

	  NodesArrayType::iterator itNode = NodeBegin + i;
	  if((itNode)->IsNot(ISOLATED)){
	    // Current step information "N+1" (before step update).
	    const double& nodal_mass                    = (itNode)->FastGetSolutionStepValue(NODAL_MASS);
	    array_1d<double,3>& current_residual        = (itNode)->FastGetSolutionStepValue(FORCE_RESIDUAL);
	    array_1d<double,3>& current_velocity        = (itNode)->FastGetSolutionStepValue(VELOCITY);
	    array_1d<double,3>& current_displacement    = (itNode)->FastGetSolutionStepValue(DISPLACEMENT); 
	    array_1d<double,3>& current_acceleration    = (itNode)->FastGetSolutionStepValue(ACCELERATION);
	  

	    // Solution of the explicit equation:
	    if((itNode)->IsFixed(ACCELERATION_X) == false)
	      current_acceleration[0] = current_residual[0]/nodal_mass;
	  
	    if((itNode)->IsFixed(ACCELERATION_Y) == false)
	      current_acceleration[1] = current_residual[1]/nodal_mass;

	    // For 3D cases
	    if ((itNode)->HasDofFor(DISPLACEMENT_Z))
	      {
		if((itNode)->IsFixed(ACCELERATION_Z) == false)
		  current_acceleration[2] = current_residual[2]/nodal_mass;
	      }
	  
	    if ((itNode)->IsFixed(VELOCITY_X)){
	      current_acceleration[0] = 0.0;
	      // middle_velocity[0]      = current_velocity[0];
	    }
	    else if((itNode)->IsFixed(DISPLACEMENT_X)){
	      current_acceleration[0] = 0.0;
	      // current_velocity[0] = 0.0;
	      // middle_velocity[0]    = 0.0;
	    }

	    if ((itNode)->IsFixed(VELOCITY_Y)){
	      current_acceleration[1] = 0.0;
	      // middle_velocity[1]      = current_velocity[1];
	    }
	    else if ((itNode)->IsFixed(DISPLACEMENT_Y)){
	      current_acceleration[1] = 0.0;
	      // current_velocity[1] = 0.0;
	      // middle_velocity[1]    = 0.0;	    
	    }
	  
	    // For 3D cases
	    if ((itNode)->HasDofFor(DISPLACEMENT_Z))
	      {
		if ((itNode)->IsFixed(VELOCITY_Z)){
		  current_acceleration[2] = 0.0;
		  // middle_velocity[2]      = current_velocity[2];
		}
		else if ((itNode)->IsFixed(DISPLACEMENT_Z)){
		  current_acceleration[2] = 0.0;
		  // current_velocity[2] = 0.0;
		  // middle_velocity[2]    = 0.0;	  	
		}
	      }
	  
	    noalias(current_velocity)     += ( mTime.Current - mTime.Previous ) * current_acceleration;
	    noalias(current_displacement) += mTime.Delta * current_velocity;    
  
	  }else{
	    // std::cout<<"in updateMomentum ISOLATED NODE "<<(itNode)->X()<<" "<<(itNode)->Y()<<std::endl;
	    array_1d<double,3>& current_velocity        = (itNode)->FastGetSolutionStepValue(VELOCITY);
	    array_1d<double,3>& current_displacement    = (itNode)->FastGetSolutionStepValue(DISPLACEMENT);
	    array_1d<double,3>& current_acceleration    = (itNode)->FastGetSolutionStepValue(ACCELERATION);
	    if((itNode)->SolutionStepsDataHas(VOLUME_ACCELERATION)){
	      array_1d<double, 3 >& VolumeAcceleration = (itNode)->FastGetSolutionStepValue(VOLUME_ACCELERATION);
	      current_acceleration = VolumeAcceleration;
	      noalias(current_velocity)     += ( mTime.Current - mTime.Previous ) * current_acceleration;
	      noalias(current_displacement) += mTime.Delta * current_velocity; 
	    }
	  }

	}

      mTime.Previous = mTime.Current;
      mTime.PreviousMiddle = mTime.Middle;

      KRATOS_CATCH("")
	}



    /**
     * Functions that calculates the RHS of a "element" object
     * @param rCurrentElement: The element to compute
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the condition degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
     */

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
				    LocalSystemVectorType& RHS_Contribution,
				    Element::EquationIdVectorType& EquationId,
				    ProcessInfo& rCurrentProcessInfo) override
    {

      KRATOS_TRY

	int thread = OpenMPUtils::ThisThread();

      //basic operations for the element considered
      (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

      if(mRayleighDamping)
	{
	  (rCurrentElement) -> CalculateDampingMatrix(mMatrix[thread], rCurrentProcessInfo);

	  AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix[thread], rCurrentProcessInfo);
	}

      //add explicit contribution of the Element Residual (RHS) to nodal Force Residual (nodal RHS)
      (rCurrentElement) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);

      KRATOS_CATCH( "" )
	}

    /**
     * Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition: The condition to compute
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the condition degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
     */
    
    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
					      LocalSystemVectorType& RHS_Contribution,
					      Element::EquationIdVectorType& EquationId,
					      ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

	//int thread = OpenMPUtils::ThisThread();

	//basic operations for the element considered
	(rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
          
      // if(mRayleighDamping)
      //    {
      // 	   (rCurrentCondition) -> CalculateDampingMatrix(mMatrix[thread], rCurrentProcessInfo);
	   
      // 	   AddDynamicsToRHS (rCurrentCondition, RHS_Contribution, mMatrix[thread], rCurrentProcessInfo);
      //    }


      //add explicit contribution of the Condition Residual (RHS) to nodal Force Residual (nodal RHS)
      (rCurrentCondition) -> AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);


      KRATOS_CATCH( "" )
	}

    
    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart
     * @return 0 all ok
     */    
    virtual int Check(ModelPart& rModelPart) override
    {
      BaseType::Check(rModelPart);
      
      // Check for variables keys
      // Verify that the variables are correctly initialized
      if(DISPLACEMENT.Key() == 0)
        {
	  KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(VELOCITY.Key() == 0)
        {
	  KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(ACCELERATION.Key() == 0)
        {
	  KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(NODAL_MASS.Key() == 0)
        {
	  KRATOS_ERROR << "NODAL_MASS has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(NODAL_ERROR.Key() == 0)
        {
     	  KRATOS_ERROR << "NODAL_ERROR has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(FORCE_RESIDUAL.Key() == 0)
        {
	  KRATOS_ERROR << "FORCE_RESIDUAL has Key zero! (check if the application is correctly registered" << std::endl;
        }

      // Check that variables are correctly allocated
      for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin();
	  it!=rModelPart.NodesEnd(); ++it)
        {
	  if (it->SolutionStepsDataHas(DISPLACEMENT) == false)
            {
	      KRATOS_ERROR << "DISPLACEMENT variable is not allocated for node " << it->Id() << std::endl;
            }
	  if (it->SolutionStepsDataHas(VELOCITY) == false)
            {
	      KRATOS_ERROR << "VELOCITY variable is not allocated for node " << it->Id() << std::endl;
            }
	  if (it->SolutionStepsDataHas(ACCELERATION) == false)
            {
	      KRATOS_ERROR << "ACCELERATION variable is not allocated for node " << it->Id() << std::endl;
            }
 	  if (it->SolutionStepsDataHas(NODAL_MASS) == false)
            {
	      KRATOS_ERROR << "NODAL_MASS variable is not allocated for node " << it->Id() << std::endl;
            }
 	  if (it->SolutionStepsDataHas(NODAL_ERROR) == false)
            {
	      KRATOS_ERROR << "NODAL_ERROR variable is not allocated for node " << it->Id() << std::endl;
            }  
	  if (it->SolutionStepsDataHas(FORCE_RESIDUAL) == false)
            {
	      KRATOS_ERROR << "FORCE_RESIDUAL variable is not allocated for node " << it->Id() << std::endl;
            }
        }
      
      // Check for minimum value of the buffer index
      // Verify buffer size
      if (rModelPart.GetBufferSize() < 2)
        {
	  KRATOS_ERROR << "insufficient buffer size. Buffer size should be greater than 2. Current size is" << rModelPart.GetBufferSize() << std::endl;
        }
      
      return 0;
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
    
    std::vector<Matrix> mMatrix;
    std::vector<Vector> mVector;

    TimeVariables       mTime;
    DeltaTimeParameters mDeltaTime;   

    bool                mRayleighDamping;
    bool                mSchemeIsInitialized;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{


    //*********************************************************************************
    // Initialize residual
    //*********************************************************************************
    
    void InitializeResidual( ModelPart& rModelPart )

    {
      KRATOS_TRY
	
	const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);
	
      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesArrayType::iterator itNode = NodeBegin + i;
  	  (itNode)->FastGetSolutionStepValue(FORCE_RESIDUAL).clear(); //RHS of momentum
  	  (itNode)->FastGetSolutionStepValue(NODAL_ERROR)=0; //RHS of continuity 
        }

      KRATOS_CATCH("")
	}
    
    //*********************************************************************************
    // Custom initialization
    //*********************************************************************************

    void CalculateDeltaTime(ModelPart& rModelPart)
    {

      KRATOS_TRY

	ProcessInfo& rCurrentProcessInfo= rModelPart.GetProcessInfo();

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      //most autors recommend a value near 0.80 (Belytschko - Nonlinear FE.. 2000. chap 6. pag. 315)
      double safety_factor     = 0.5;  
      double stable_delta_time = 0.00;
      Vector delta_times(NumThreads);
     
      for(unsigned int i = 0; i < NumThreads; i++)
	delta_times[i] = mDeltaTime.Maximum/safety_factor;
      
      OpenMPUtils::PartitionVector ElementPartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(), NumThreads, ElementPartition);
	
      const int nelements = static_cast<int>(rModelPart.Elements().size());
      ElementsArrayType::iterator ElementBegin = rModelPart.Elements().begin();

#pragma omp parallel for firstprivate(ElementBegin) private(stable_delta_time)
      for(int i = 0;  i < nelements; i++)
        {
	  ElementsArrayType::iterator itElement = ElementBegin + i;
	  //get geometric and material properties
	  double inRadius   = (itElement)->GetGeometry().Inradius();
	  double waveSpeed  = 1450;

	  stable_delta_time = 0.9*inRadius/waveSpeed;

	  if(stable_delta_time > 0.00)
	    {
	      int thread = OpenMPUtils::ThisThread();
	      if(stable_delta_time < delta_times[thread])
		{
		  delta_times[thread] = stable_delta_time;
		}
	    }

	}

      
      stable_delta_time  = (*std::min_element(delta_times.begin(), delta_times.end()));
      stable_delta_time *= safety_factor;// * 0.5; //extra factor added to get an stable delta time

      double current_delta_time = rCurrentProcessInfo[DELTA_TIME];

      // std::cout<< " ATTENTION!!!!! Stable delta time is "<< stable_delta_time <<" and current delta time is "<< current_delta_time<<std::endl;
      
      if(stable_delta_time<(0.99999*current_delta_time)){
	std::cout<< " ATTENTION!!!!! Stable delta time is "<< stable_delta_time <<" and current delta time is "<< current_delta_time<<std::endl;
      }
      
      if(stable_delta_time < mDeltaTime.Maximum){
       	rCurrentProcessInfo[DELTA_TIME] = stable_delta_time;	  
      }
      else{
       	if( current_delta_time > mDeltaTime.Maximum/safety_factor )
       	  rCurrentProcessInfo[DELTA_TIME] = mDeltaTime.Maximum;
      }
      
      // std::cout<< "  [EXPLICIT PREDICTION LEVEL"<<mDeltaTime.PredictionLevel<<"]:(computed stable time step = "<< stable_delta_time <<" s)"<< std::endl;
      // std::cout<< "  Using  = "<< rCurrentProcessInfo[DELTA_TIME] <<" s as time step DELTA_TIME)"<< std::endl;
        
      KRATOS_CATCH("")
	}


    // void CalculateDeltaTime(ModelPart& rModelPart)
    // {

    //   KRATOS_TRY

    //   ProcessInfo& rCurrentProcessInfo= rModelPart.GetProcessInfo();

    //   const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

    //   //most autors recommend a value near 0.80 (Belytschko - Nonlinear FE.. 2000. chap 6. pag. 315)
    //   double safety_factor     = 0.5;  
    //   double stable_delta_time = 0.00;
    //   Vector delta_times(NumThreads);
     
    //   for(unsigned int i = 0; i < NumThreads; i++)
    // 	delta_times[i] = mDeltaTime.Maximum/safety_factor;
      
    //   OpenMPUtils::PartitionVector ElementPartition;
    //   OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(), NumThreads, ElementPartition);
	
    //   const int nelements = static_cast<int>(rModelPart.Elements().size());
    //   ElementsArrayType::iterator ElementBegin = rModelPart.Elements().begin();

    //   #pragma omp parallel for firstprivate(ElementBegin) private(stable_delta_time)
    //   for(int i = 0;  i < nelements; i++)
    //     {
    // 	  ElementsArrayType::iterator itElement = ElementBegin + i;
    // 	  //get geometric and material properties
    // 	  double length   = (itElement)->GetGeometry().Length();
    // 	  // double alpha    = (itElement)->GetProperties()[RAYLEIGH_ALPHA];
    // 	  // double beta     = (itElement)->GetProperties()[RAYLEIGH_BETA];
    // 	  double alpha = 0;
    // 	  double beta = 0;
    // 	  double E        = (itElement)->GetProperties()[YOUNG_MODULUS];
    // 	  double v        = (itElement)->GetProperties()[POISSON_RATIO];
    // 	  double ro       = (itElement)->GetProperties()[DENSITY];

    // 	  //compute courant criterion
    // 	  double bulk       = E/(3.0*(1.0-2.0*v));               
    // 	  double wavespeed  = sqrt(bulk/ro);
    // 	  double w          = 2.0*wavespeed/length;   //frequency

    // 	  double psi        = 0.5*(alpha/w + beta*w); //critical ratio;
    // 	  stable_delta_time = (2.0/w)*(sqrt(1.0 + psi*psi)-psi);

    // 	  if(stable_delta_time > 0.00)
    // 	    {
    // 	      int thread = OpenMPUtils::ThisThread();
    // 	      if(stable_delta_time < delta_times[thread])
    // 		{
    // 		  delta_times[thread] = stable_delta_time;
    // 		}
    // 	    }

    // 	}

    //   stable_delta_time  = (*std::min_element(delta_times.begin(), delta_times.end()));
    //   stable_delta_time *= safety_factor;// * 0.5; //extra factor added to get an stable delta time

    //   double current_delta_time = rCurrentProcessInfo[DELTA_TIME];
      
    //   if(stable_delta_time < mDeltaTime.Maximum){
    // 	rCurrentProcessInfo[DELTA_TIME] = stable_delta_time;	  
    //   }
    //   else{
    // 	if( current_delta_time > mDeltaTime.Maximum/safety_factor )
    // 	  rCurrentProcessInfo[DELTA_TIME] = mDeltaTime.Maximum;
    //   }

    //   std::cout<< "  [EXPLICIT PREDICTION LEVEL"<<mDeltaTime.PredictionLevel<<"]:(computed stable time step = "<< stable_delta_time <<" s)"<< std::endl;
    //   std::cout<< "  Using  = "<< rCurrentProcessInfo[DELTA_TIME] <<" s as time step DELTA_TIME)"<< std::endl;
        
    //   KRATOS_CATCH("")
    // }

    
    //*********************************************************************************
    // Custom initialization
    //*********************************************************************************

    void InitializeExplicitScheme(ModelPart& rModelPart)
    {
      KRATOS_TRY

	const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);
	
      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesArrayType::iterator itNode = NodeBegin + i;
		  
	  // array_1d<double,3>& middle_velocity       = (itNode)->FastGetSolutionStepValue(MIDDLE_VELOCITY);
	  // array_1d<double,3>& current_velocity      = (itNode)->FastGetSolutionStepValue(VELOCITY);
	  array_1d<double,3>& current_residual      = (itNode)->FastGetSolutionStepValue(FORCE_RESIDUAL);
	  array_1d<double,3>& current_displacement  = (itNode)->FastGetSolutionStepValue(DISPLACEMENT);
          
	  // noalias(middle_velocity) = current_velocity;
	  current_residual.clear();
	  current_displacement.clear();
	}
      
      KRATOS_CATCH("")
	}

    //*********************************************************************************
    // Custom initialization
    //*********************************************************************************
    
    void SchemeCustomInitialization(ModelPart& rModelPart)
    {
      KRATOS_TRY

	const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);
	
      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesArrayType::iterator itNode = NodeBegin + i;
	  
	  // Current step information "N+1" (before step update).
	  const double& nodal_mass                 = (itNode)->FastGetSolutionStepValue(NODAL_MASS);
	  array_1d<double,3>& current_residual     = (itNode)->FastGetSolutionStepValue(FORCE_RESIDUAL);
	  
	  array_1d<double,3>& current_velocity     = (itNode)->FastGetSolutionStepValue(VELOCITY);
	  array_1d<double,3>& current_displacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT);
	  // array_1d<double,3>& middle_velocity      = (itNode)->FastGetSolutionStepValue(MIDDLE_VELOCITY);
	  
	  array_1d<double,3>& current_acceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION);
	  
	  if((itNode)->IsFixed(ACCELERATION_X) == false)
	    current_acceleration[0] = current_residual[0]/nodal_mass;
	  
	  if((itNode)->IsFixed(ACCELERATION_Y) == false)
	    current_acceleration[1] = current_residual[1]/nodal_mass;

	  // For 3D cases
	  if ((itNode)->HasDofFor(DISPLACEMENT_Z))
	    {
	      if((itNode)->IsFixed(ACCELERATION_Z) == false)
		current_acceleration[2] = current_residual[2]/nodal_mass;
	    }

	  if ((itNode)->IsFixed(VELOCITY_X)){
	    current_acceleration[0] = 0.0;
	    // middle_velocity[0]      = current_velocity[0];
	  }
	  else if((itNode)->IsFixed(DISPLACEMENT_X)){
	    current_acceleration[0] = 0.0;
	    current_velocity[0] = 0.0;
	    // middle_velocity[0]    = 0.0;
	  }

	  if ((itNode)->IsFixed(VELOCITY_Y)){
	    current_acceleration[1] = 0.0;
	    // middle_velocity[1]      = current_velocity[1];
	  }
	  else if ((itNode)->IsFixed(DISPLACEMENT_Y)){
	    current_acceleration[1] = 0.0;
	    current_velocity[1] = 0.0;
	    // middle_velocity[1]    = 0.0;	    
	  }
	  
	  // For 3D cases
	  if ((itNode)->HasDofFor(DISPLACEMENT_Z))
	    {
	      if ((itNode)->IsFixed(VELOCITY_Z)){
	  	current_acceleration[2] = 0.0;
	  	// middle_velocity[2]      = current_velocity[2];
	      }
	      else if ((itNode)->IsFixed(DISPLACEMENT_Z)){
		current_acceleration[2] = 0.0;
		current_velocity[2] = 0.0;
		// middle_velocity[2]    = 0.0;	  	
	      }
	    }

	  //Solution of the explicit equation:
	  // noalias(middle_velocity)   = ( mTime.Middle - mTime.Previous ) * current_acceleration;
	  // noalias(current_velocity)  = middle_velocity + ( mTime.Previous - mTime.PreviousMiddle ) * current_acceleration;
	  noalias(current_velocity)  += ( mTime.Current - mTime.Previous ) * current_acceleration;
	  current_displacement.clear();
          
	}

      mTime.Previous = mTime.Current;
      mTime.PreviousMiddle = mTime.Middle;

      KRATOS_CATCH("")
	}
    

    
    /**
     * It adds the dynamic RHS contribution of the condition: b - D*v
     * @param rCurrentElement: The element to compute
     * @param RHS_Contribution: The dynamic contribution for the RHS
     * @param D: The damping matrix
     * @param CurrentProcessInfo: The current process info instance
     */

    void AddDynamicsToRHS(Element::Pointer rCurrentElement,
			  LocalSystemVectorType& RHS_Contribution,
			  LocalSystemMatrixType& D,
			  ProcessInfo& CurrentProcessInfo)
    {
      int thread = OpenMPUtils::ThisThread();

      // Adding damping contribution
      if (D.size1() != 0)
	{
	  // this->GetFirstDerivativesVector(rCurrentElement, mVector[thread]);

	  noalias(RHS_Contribution) -= prod(D, mVector[thread]);
	}
    }

    /**
     * It adds the dynamic RHS contribution of the condition: b - D*v
     * @param rCurrentCondition: The condition to compute
     * @param RHS_Contribution: The dynamic contribution for the RHS
     * @param D: The damping matrix
     * @param CurrentProcessInfo: The current process info instance
     */
   
    void AddDynamicsToRHS(Condition::Pointer rCurrentCondition,
			  LocalSystemVectorType& RHS_Contribution,
			  LocalSystemMatrixType& D,
			  ProcessInfo& CurrentProcessInfo)
    {
      int thread = OpenMPUtils::ThisThread();

      // Adding damping contribution
      if (D.size1() != 0)
	{
	  // this->GetFirstDerivativesVector(rCurrentCondition, mVector[thread]);

	  noalias(RHS_Contribution) -= prod(D, mVector[thread]);
	}      
    }


    /**
     *  Obtain explicit first derivatives (velocity)
     */
    
    // void GetFirstDerivativesVector(Element::Pointer rCurrentElement, Vector& rValues ) //V at time n-1/2 old
    // {

    //   const unsigned int number_of_nodes = rCurrentElement->GetGeometry().size();
    //   const unsigned int dimension       = rCurrentElement->GetGeometry().WorkingSpaceDimension();
    //   unsigned int       element_size    = number_of_nodes * dimension;

    //   if ( rValues.size() != element_size )
    // 	rValues.resize( element_size, false );

    //   for ( unsigned int i = 0; i < number_of_nodes; i++ )
    // 	{
    // 	  unsigned int index = i * dimension;


    // 	  rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue(MIDDLE_VELOCITY);

    // 	  rValues[index]     = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[0];
    // 	  rValues[index + 1] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[1];

    // 	  if ( dimension == 3 )
    // 	    rValues[index + 2] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[2];
    // 	}
    // }

    
    /**
     *  Obtain explicit first derivatives (velocity)
     */
    
    // void GetFirstDerivativesVector(Condition::Pointer rCurrentCondition, Vector& rValues ) //V at time n-1/2 old
    // {

    //   const unsigned int number_of_nodes = rCurrentCondition->GetGeometry().size();
    //   const unsigned int dimension       = rCurrentCondition->GetGeometry().WorkingSpaceDimension();
    //   unsigned int       condition_size  = number_of_nodes * dimension;
      
    //   if ( rValues.size() != condition_size )
    // 	rValues.resize( condition_size, false );

    //   for ( unsigned int i = 0; i < number_of_nodes; i++ )
    // 	{
    // 	  unsigned int index = i * dimension;

    // 	  rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue(MIDDLE_VELOCITY);

    // 	  rValues[index]     = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[0];
    // 	  rValues[index + 1] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[1];

    // 	  if ( dimension == 3 )
    // 	    rValues[index + 2] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[2];
    // 	}    
    // }
    
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
    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
  }; /* Class FirstOrderForwardEulerScheme */
  ///@}
  ///@name Type Definitions
  ///@{
  ///@}
  ///@name Input and output
  ///@{
  ///@}  
}  /* namespace Kratos.*/

#endif /* KRATOS_FIRST_ORDER_FORWARD_EULER_SCHEME  defined */

