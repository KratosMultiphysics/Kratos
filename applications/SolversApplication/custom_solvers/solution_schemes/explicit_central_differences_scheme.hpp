//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:           MSantasusana $
//   Last modified by:    $Co-Author:         JMCarbonell $
//   Date:                $Date:               April 2014 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME_H_INCLUDED)
#define  KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/solution_scheme.hpp"

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
  class ExplicitCentralDifferencesScheme : public SolutionScheme<TSparseSpace,TDenseSpace>
  {
  public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ExplicitCentralDifferencesScheme );

    typedef SolutionScheme<TSparseSpace,TDenseSpace>                             BaseType;
    typedef typename BaseType::SolutionSchemePointerType                  BasePointerType;
    typedef typename BaseType::LocalFlagType                                LocalFlagType;

    typedef typename BaseType::DofsArrayType                                DofsArrayType;
    typedef typename BaseType::SystemMatrixType                          SystemMatrixType;
    typedef typename BaseType::SystemVectorType                          SystemVectorType;
    typedef typename BaseType::LocalSystemVectorType                LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType                LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                              NodesContainerType;
    typedef ModelPart::ElementsContainerType                        ElementsContainerType;
    typedef ModelPart::ConditionsContainerType                    ConditionsContainerType;


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
    ExplicitCentralDifferencesScheme(Flags& rOptions,
                                     const double  rMaximumDeltaTime,
				     const double  rDeltaTimeFraction,
				     const double  rDeltaTimePredictionLevel)
      : BaseType(rOptions)
    {

      mDeltaTime.PredictionLevel  = rDeltaTimePredictionLevel;

      mDeltaTime.Maximum          = rMaximumDeltaTime;

      mDeltaTime.Fraction         = rDeltaTimeFraction;

      // Allocate auxiliary memory
      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      mMatrix.resize(NumThreads);
      mVector.resize(NumThreads);

      this->Set(LocalFlagType::INITIALIZED, false);
    }

    ///Copy constructor
    ExplicitCentralDifferencesScheme(ExplicitCentralDifferencesScheme& rOther)
        :BaseType(rOther)
        ,mMatrix(rOther.mMatrix)
        ,mVector(rOther.mVector)
	,mDeltaTime(rOther.mDeltaTime)
    {
    }

    /// Clone
    BasePointerType Clone() override
    {
        return BasePointerType( new ExplicitCentralDifferencesScheme(*this) );
    }

    /// Destructor
    ~ExplicitCentralDifferencesScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    void Initialize(ModelPart& rModelPart) override
    {
      KRATOS_TRY

          if( (mDeltaTime.PredictionLevel>0) && this->IsNot(LocalFlagType::INITIALIZED) )
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

      if(this->IsNot(LocalFlagType::INITIALIZED))
	{
	  this->InitializeExplicitScheme(rModelPart);
	}
      else
	{
	  this->SchemeCustomInitialization(rModelPart);
	}

      this->Set(LocalFlagType::INITIALIZED, true);

      KRATOS_CATCH("")
    }


    /**
     * It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart: The model of the problem to solve
     * @param rA: LHS matrix
     * @param rDx: Incremental update of primary variables
     * @param rb: RHS Vector
     *
     */
    void InitializeSolutionStep(ModelPart& rModelPart) override
    {
      KRATOS_TRY

      BaseType::InitializeSolutionStep(rModelPart);

      if(mDeltaTime.PredictionLevel>1)
	{
	  this->CalculateDeltaTime(rModelPart);
	}

      this->InitializeResidual(rModelPart);

      KRATOS_CATCH("")
    }

     /**
     * Performing the update of the solution.
     * incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param rModelPart
     * @param rDofSet set of all primary variables
     * @param rDx incremental update of primary variables
     */
    void Update(ModelPart& rModelPart,
		DofsArrayType& rDofSet,
		SystemVectorType& rDx) override
    {
      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo  = rModelPart.GetProcessInfo();

      //Step Update
      mTime.Current   = rCurrentProcessInfo[TIME];  //the first step is (time = initial_time + delta time )
      mTime.Delta     = rCurrentProcessInfo[DELTA_TIME];
      mTime.Middle    = 0.5 * ( mTime.Previous + mTime.Current );

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

      #pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesContainerType::iterator itNode = NodeBegin + i;

	  // Current step information "N+1" (before step update).
	  const double& nodal_mass                    = (itNode)->FastGetSolutionStepValue(NODAL_MASS);
	  array_1d<double,3>& current_residual        = (itNode)->FastGetSolutionStepValue(FORCE_RESIDUAL);

	  array_1d<double,3>& current_velocity        = (itNode)->FastGetSolutionStepValue(VELOCITY);
	  array_1d<double,3>& current_displacement    = (itNode)->FastGetSolutionStepValue(DISPLACEMENT);
	  array_1d<double,3>& middle_velocity         = (itNode)->FastGetSolutionStepValue(MIDDLE_VELOCITY);

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
	    middle_velocity[0]      = current_velocity[0];
	  }
	  else if((itNode)->IsFixed(DISPLACEMENT_X)){
	    current_acceleration[0] = 0.0;
	    middle_velocity[0]    = 0.0;
	  }

	  if ((itNode)->IsFixed(VELOCITY_Y)){
	    current_acceleration[1] = 0.0;
	    middle_velocity[1]      = current_velocity[1];
	  }
	  else if ((itNode)->IsFixed(DISPLACEMENT_Y)){
	    current_acceleration[1] = 0.0;
	    middle_velocity[1]    = 0.0;
	  }

	  // For 3D cases
	  if ((itNode)->HasDofFor(DISPLACEMENT_Z))
	    {
	      if ((itNode)->IsFixed(VELOCITY_Z)){
	  	current_acceleration[2] = 0.0;
	  	middle_velocity[2]      = current_velocity[2];
	      }
	      else if ((itNode)->IsFixed(DISPLACEMENT_Z)){
		current_acceleration[2] = 0.0;
		middle_velocity[2]    = 0.0;
	      }
	    }

	  noalias(current_velocity)      = middle_velocity + ( mTime.Previous - mTime.PreviousMiddle ) * current_acceleration;
	  noalias(middle_velocity)       = current_velocity + ( mTime.Middle - mTime.Previous ) * current_acceleration;
	  noalias(current_displacement) += mTime.Delta * middle_velocity;

	}

      mTime.Previous = mTime.Current;
      mTime.PreviousMiddle = mTime.Middle;

      KRATOS_CATCH("")
    }



    /**
     * Functions that calculates the RHS of a "element" object
     * @param rCurrentElement: The element to compute
     * @param rRHS_Contribution: The RHS vector contribution
     * @param rEquationId: The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
				    LocalSystemVectorType& rRHS_Contribution,
				    Element::EquationIdVectorType& rEquationId,
				    ProcessInfo& rCurrentProcessInfo) override
    {

      KRATOS_TRY

      int thread = OpenMPUtils::ThisThread();

      //basic operations for the element considered
      (rCurrentElement) -> CalculateRightHandSide(rRHS_Contribution,rCurrentProcessInfo);

      if(this->mOptions.Is(LocalFlagType::RAYLEIGH_DAMPING))
	{
	  (rCurrentElement) -> CalculateDampingMatrix(mMatrix[thread], rCurrentProcessInfo);

	  AddDynamicsToRHS (rCurrentElement, rRHS_Contribution, mMatrix[thread], rCurrentProcessInfo);
	}

      //add explicit contribution of the Element Residual (RHS) to nodal Force Residual (nodal RHS)
      (rCurrentElement) -> AddExplicitContribution(rRHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);

      KRATOS_CATCH("")
    }

    /**
     * Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition: The condition to compute
     * @param rRHS_Contribution: The RHS vector contribution
     * @param rEquationId: The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
					      LocalSystemVectorType& rRHS_Contribution,
					      Element::EquationIdVectorType& rEquationId,
					      ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      //int thread = OpenMPUtils::ThisThread();

      //basic operations for the element considered
      (rCurrentCondition) -> CalculateRightHandSide(rRHS_Contribution,rCurrentProcessInfo);

      // if(this->mOptions.Is(RAYLEIGH_DAMPING))
      //    {
      // 	   (rCurrentCondition) -> CalculateDampingMatrix(mMatrix[thread], rCurrentProcessInfo);

      // 	   AddDynamicsToRHS (rCurrentCondition, RHS_Contribution, mMatrix[thread], rCurrentProcessInfo);
      //    }


      //add explicit contribution of the Condition Residual (RHS) to nodal Force Residual (nodal RHS)
      (rCurrentCondition) -> AddExplicitContribution(rRHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);


      KRATOS_CATCH("")
    }


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart
     * @return 0 all ok
     */
    int Check(ModelPart& rModelPart) override
    {
      KRATOS_TRY

      // Perform base base checks
      int ErrorCode = 0;
      ErrorCode  = BaseType::Check(rModelPart);

      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
      KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
      KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
      KRATOS_CHECK_VARIABLE_KEY(NODAL_MASS);
      KRATOS_CHECK_VARIABLE_KEY(MIDDLE_VELOCITY);
      KRATOS_CHECK_VARIABLE_KEY(FORCE_RESIDUAL);

      // Check that variables are correctly allocated
      for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); ++it)
        {
	  // Nodal data
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_MASS,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MIDDLE_VELOCITY,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(FORCE_RESIDUAL,(*it));

	  // Nodal dofs
	  KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,(*it));
	  KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,(*it));
	  if( rModelPart.GetProcessInfo()[SPACE_DIMENSION] == 3 )
	    KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,(*it));
        }

      // Check for minimum value of the buffer index
      // Verify buffer size
      if (rModelPart.GetBufferSize() < 2)
        {
	  KRATOS_ERROR << "insufficient buffer size. Buffer size should be greater than 2. Current size is" << rModelPart.GetBufferSize() << std::endl;
        }

      return ErrorCode;

      KRATOS_CATCH("")
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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void InitializeResidual( ModelPart& rModelPart )

    {
      KRATOS_TRY

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

      #pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesContainerType::iterator itNode = NodeBegin + i;
  	  (itNode)->FastGetSolutionStepValue(FORCE_RESIDUAL).clear();
        }

      KRATOS_CATCH("")
    }

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
      ElementsContainerType::iterator ElementBegin = rModelPart.Elements().begin();

      #pragma omp parallel for firstprivate(ElementBegin) private(stable_delta_time)
      for(int i = 0;  i < nelements; i++)
        {
	  ElementsContainerType::iterator itElement = ElementBegin + i;
	  //get geometric and material properties
	  double length   = (itElement)->GetGeometry().Length();
	  double alpha    = (itElement)->GetProperties()[RAYLEIGH_ALPHA];
	  double beta     = (itElement)->GetProperties()[RAYLEIGH_BETA];
	  double E        = (itElement)->GetProperties()[YOUNG_MODULUS];
	  double v        = (itElement)->GetProperties()[POISSON_RATIO];
	  double ro       = (itElement)->GetProperties()[DENSITY];

	  //compute courant criterion
	  double bulk       = E/(3.0*(1.0-2.0*v));
	  double wavespeed  = sqrt(bulk/ro);
	  double w          = 2.0*wavespeed/length;   //frequency

	  double psi        = 0.5*(alpha/w + beta*w); //critical ratio;
	  stable_delta_time = (2.0/w)*(sqrt(1.0 + psi*psi)-psi);

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

      if(stable_delta_time < mDeltaTime.Maximum){
	rCurrentProcessInfo[DELTA_TIME] = stable_delta_time;
      }
      else{
	if( current_delta_time > mDeltaTime.Maximum/safety_factor )
	  rCurrentProcessInfo[DELTA_TIME] = mDeltaTime.Maximum;
      }

      KRATOS_INFO("Time Step prediction") << "  [EXPLICIT PREDICTION LEVEL"<<mDeltaTime.PredictionLevel<<"]:(computed stable time step = "<< stable_delta_time <<" s)"<< std::endl;
      KRATOS_INFO("Using") << rCurrentProcessInfo[DELTA_TIME] <<" s as time step DELTA_TIME)"<< std::endl;

      KRATOS_CATCH("")
    }


    void InitializeExplicitScheme(ModelPart& rModelPart)
    {
      KRATOS_TRY

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

      #pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesContainerType::iterator itNode = NodeBegin + i;

	  array_1d<double,3>& middle_velocity       = (itNode)->FastGetSolutionStepValue(MIDDLE_VELOCITY);
	  array_1d<double,3>& current_velocity      = (itNode)->FastGetSolutionStepValue(VELOCITY);
	  array_1d<double,3>& current_residual      = (itNode)->FastGetSolutionStepValue(FORCE_RESIDUAL);
	  array_1d<double,3>& current_displacement  = (itNode)->FastGetSolutionStepValue(DISPLACEMENT);

	  noalias(middle_velocity) = current_velocity;
	  current_residual.clear();
	  current_displacement.clear();
	}

      KRATOS_CATCH("")
    }


    void SchemeCustomInitialization(ModelPart& rModelPart)
    {
      KRATOS_TRY

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

      #pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesContainerType::iterator itNode = NodeBegin + i;

	  // Current step information "N+1" (before step update).
	  const double& nodal_mass                 = (itNode)->FastGetSolutionStepValue(NODAL_MASS);
	  array_1d<double,3>& current_residual     = (itNode)->FastGetSolutionStepValue(FORCE_RESIDUAL);

	  array_1d<double,3>& current_velocity     = (itNode)->FastGetSolutionStepValue(VELOCITY);
	  array_1d<double,3>& current_displacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT);
	  array_1d<double,3>& middle_velocity      = (itNode)->FastGetSolutionStepValue(MIDDLE_VELOCITY);

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
	    middle_velocity[0]      = current_velocity[0];
	  }
	  else if((itNode)->IsFixed(DISPLACEMENT_X)){
	    current_acceleration[0] = 0.0;
	    middle_velocity[0]    = 0.0;
	  }

	  if ((itNode)->IsFixed(VELOCITY_Y)){
	    current_acceleration[1] = 0.0;
	    middle_velocity[1]      = current_velocity[1];
	  }
	  else if ((itNode)->IsFixed(DISPLACEMENT_Y)){
	    current_acceleration[1] = 0.0;
	    middle_velocity[1]    = 0.0;
	  }

	  // For 3D cases
	  if ((itNode)->HasDofFor(DISPLACEMENT_Z))
	    {
	      if ((itNode)->IsFixed(VELOCITY_Z)){
	  	current_acceleration[2] = 0.0;
	  	middle_velocity[2]      = current_velocity[2];
	      }
	      else if ((itNode)->IsFixed(DISPLACEMENT_Z)){
		current_acceleration[2] = 0.0;
		middle_velocity[2]    = 0.0;
	      }
	    }

	  //Solution of the explicit equation:
	  noalias(middle_velocity)   = ( mTime.Middle - mTime.Previous ) * current_acceleration;
	  noalias(current_velocity)  = middle_velocity + ( mTime.Previous - mTime.PreviousMiddle ) * current_acceleration;
	  current_displacement.clear();

	}

      mTime.Previous = mTime.Current;
      mTime.PreviousMiddle = mTime.Middle;

      KRATOS_CATCH("")
    }



    /**
     * It adds the dynamic RHS contribution of the condition: b - D*v
     * @param rCurrentElement: The element to compute
     * @param rRHS_Contribution: The dynamic contribution for the RHS
     * @param rD: The damping matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    void AddDynamicsToRHS(Element::Pointer rCurrentElement,
			  LocalSystemVectorType& rRHS_Contribution,
			  LocalSystemMatrixType& rD,
			  ProcessInfo& rCurrentProcessInfo)
    {
      int thread = OpenMPUtils::ThisThread();

      // Adding damping contribution
      if (rD.size1() != 0)
	{
	  this->GetFirstDerivativesVector(rCurrentElement, mVector[thread]);

	  noalias(rRHS_Contribution) -= prod(rD, mVector[thread]);
	}
    }

    /**
     * It adds the dynamic RHS contribution of the condition: b - D*v
     * @param rCurrentCondition: The condition to compute
     * @param rRHS_Contribution: The dynamic contribution for the RHS
     * @param rD: The damping matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    void AddDynamicsToRHS(Condition::Pointer rCurrentCondition,
			  LocalSystemVectorType& rRHS_Contribution,
			  LocalSystemMatrixType& rD,
			  ProcessInfo& rCurrentProcessInfo)
    {
      int thread = OpenMPUtils::ThisThread();

      // Adding damping contribution
      if (rD.size1() != 0)
	{
	  this->GetFirstDerivativesVector(rCurrentCondition, mVector[thread]);

	  noalias(rRHS_Contribution) -= prod(rD, mVector[thread]);
	}
    }


    /**
     *  Obtain explicit first derivatives (velocity)
     */

    void GetFirstDerivativesVector(Element::Pointer rCurrentElement, Vector& rValues ) //V at time n-1/2 old
    {

      const unsigned int number_of_nodes = rCurrentElement->GetGeometry().size();
      const unsigned int dimension       = rCurrentElement->GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension;

      if ( rValues.size() != element_size )
	rValues.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  unsigned int index = i * dimension;


	  rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue(MIDDLE_VELOCITY);

	  rValues[index]     = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[0];
	  rValues[index + 1] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[1];

	  if ( dimension == 3 )
	    rValues[index + 2] = rCurrentElement->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[2];
	}
    }


    /**
     *  Obtain explicit first derivatives (velocity)
     */

    void GetFirstDerivativesVector(Condition::Pointer rCurrentCondition, Vector& rValues ) //V at time n-1/2 old
    {

      const unsigned int number_of_nodes = rCurrentCondition->GetGeometry().size();
      const unsigned int dimension       = rCurrentCondition->GetGeometry().WorkingSpaceDimension();
      unsigned int       condition_size  = number_of_nodes * dimension;

      if ( rValues.size() != condition_size )
	rValues.resize( condition_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  unsigned int index = i * dimension;

	  rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue(MIDDLE_VELOCITY);

	  rValues[index]     = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[0];
	  rValues[index + 1] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[1];

	  if ( dimension == 3 )
	    rValues[index + 2] = rCurrentCondition->GetGeometry()[i].FastGetSolutionStepValue( MIDDLE_VELOCITY )[2];
	}
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
  }; // Class ExplicitCentralDifferencesScheme
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  //namespace Kratos.

#endif // KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME_H_INCLUDED  defined

