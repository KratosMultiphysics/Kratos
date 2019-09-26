//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:              MSantasusana $
//   Last modified by:    $Co-Author:            JMCarbonell $
//   Date:                $Date:                  April 2014 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_EXPLICIT_TWO_STEP_V_P_STRATEGY)
#define KRATOS_EXPLICIT_TWO_STEP_V_P_STRATEGY

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"

//default builder and solver
#include "custom_strategies/builders_and_solvers/explicit_two_step_v_p_builder_and_solver.hpp"

namespace Kratos
{

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ExplicitTwoStepVPStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
  /** Counted pointer of ClassName */

  KRATOS_CLASS_POINTER_DEFINITION(ExplicitTwoStepVPStrategy);

  typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

  typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

  typedef typename BaseType::TDataType TDataType;

  typedef TSparseSpace SparseSpaceType;

  typedef typename BaseType::TSchemeType TSchemeType;

  typedef typename BaseType::DofsArrayType DofsArrayType;

  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

  typedef typename BaseType::TSystemVectorType TSystemVectorType;

  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

  typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

  typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

  typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

  typedef ModelPart::NodesContainerType NodesArrayType;

  typedef ModelPart::ElementsContainerType ElementsArrayType;

  typedef ModelPart::ConditionsContainerType ConditionsArrayType;

  /** Constructors.
     */
  ExplicitTwoStepVPStrategy(
      ModelPart &model_part,
      bool MoveMeshFlag = true)
      : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
  {
  }

  ExplicitTwoStepVPStrategy(
      ModelPart &model_part,
      typename TSchemeType::Pointer pScheme,
      typename TLinearSolver::Pointer pNewLinearSolver,
      bool CalculateReactions = false,
      bool ReformDofSetAtEachStep = true,
      bool MoveMeshFlag = true)
      : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
  {
    KRATOS_TRY

    // std::cout<<"ExplicitTwoStepVPStrategy "<<std::endl;

    // typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;

    //set flags to default values
    // mCalculateReactionsFlag = CalculateReactions;
    // mReformDofSetAtEachStep = ReformDofSetAtEachStep;

    mCalculateReactionsFlag = false;
    mReformDofSetAtEachStep = true;

    //saving the scheme
    mpScheme = pScheme;

    //saving the linear solver
    mpLinearSolver = pNewLinearSolver; //Not used in explicit strategies

    //setting up the default builder and solver
    mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
        new ExplicitTwoStepVPBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(mpLinearSolver));

    //set flags to start correcty the calculations
    mSolutionStepIsInitialized = false;
    mInitializeWasPerformed = false;

    //set EchoLevel to the default value (only time is displayed)
    SetEchoLevel(1);

    //set RebuildLevel to the default value
    BaseType::SetRebuildLevel(0);

    //set it true for explicit :: taking the default geometry lumping factors
    BaseType::GetModelPart().GetProcessInfo()[COMPUTE_LUMPED_MASS_MATRIX] = true;

    // BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(mpLinearSolver, PRESSURE));

    KRATOS_CATCH("")
  }

  /** Destructor.
     */
  virtual ~ExplicitTwoStepVPStrategy()
  {
  }

  /** Destructor.
     */

  //Set and Get Scheme ... containing Builder, Update and other

  void SetScheme(typename TSchemeType::Pointer pScheme)
  {
    mpScheme = pScheme;
  };

  typename TSchemeType::Pointer GetScheme()
  {
    return mpScheme;
  };

  //Set and Get the BuilderAndSolver

  void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
  {
    mpBuilderAndSolver = pNewBuilderAndSolver;
  };

  typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
  {
    return mpBuilderAndSolver;
  };

  //Ser and Get Flags

  void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
  {
    mInitializeWasPerformed = InitializePerformedFlag;
  }

  bool GetInitializePerformedFlag()
  {
    return mInitializeWasPerformed;
  }

  void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
  {
    mCalculateReactionsFlag = CalculateReactionsFlag;
  }

  bool GetCalculateReactionsFlag()
  {
    return mCalculateReactionsFlag;
  }

  void SetReformDofSetAtEachStepFlag(bool flag)
  {

    mReformDofSetAtEachStep = flag;
  }

  bool GetReformDofSetAtEachStepFlag()
  {
    return mReformDofSetAtEachStep;
  }

  //level of echo for the solving strategy
  // 0 -> mute... no echo at all
  // 1 -> printing time and basic informations
  // 2 -> printing linear solver data
  // 3 -> Print of debug informations:
  //    Echo of stiffness matrix, Dx, b...

  void SetEchoLevel(int Level) override
  {
    BaseType::mEchoLevel = Level;
  }

  //**********************************************************************
  //**********************************************************************

  void Initialize() override
  {
    KRATOS_TRY

    //pointers needed in the solution
    typename TSchemeType::Pointer pScheme = GetScheme();
    typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
    ModelPart &r_model_part = BaseType::GetModelPart();

    TSystemMatrixType mA = TSystemMatrixType(); //dummy initialization. Not used in builder and solver.

    // std::cout<<" INITIALIZE  "<<std::endl;

    //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
    if (pScheme->SchemeIsInitialized() == false)
      pScheme->Initialize(BaseType::GetModelPart());

    //Initialize The Elements - OPERATIONS TO BE DONE ONCE
    if (pScheme->ElementsAreInitialized() == false)
      pScheme->InitializeElements(BaseType::GetModelPart());

    //Initialize The Conditions- OPERATIONS TO BE DONE ONCE
    if (pScheme->ConditionsAreInitialized() == false)
      pScheme->InitializeConditions(BaseType::GetModelPart());

    pBuilderAndSolver->BuildLHS(pScheme, r_model_part, mA); //calculate Mass Matrix
    TSystemVectorType mb = TSystemVectorType();
    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb); //calculate Residual for scheme initialization

    mInitializeWasPerformed = true;

    // this->UpdatePressure(r_model_part);

    //std::cout<<" Rebuild Level "<<BaseType::mRebuildLevel<<std::endl;

    KRATOS_CATCH("")
  }

  //**********************************************************************
  //**********************************************************************

  void InitializeSolutionStep() override
  {
    KRATOS_TRY

    // std::cout<<" InitializeSolutionStep() in explicit two step v p"<<std::endl;

    typename TSchemeType::Pointer pScheme = GetScheme();
    typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
    ModelPart &r_model_part = BaseType::GetModelPart();

    TSystemMatrixType mA = TSystemMatrixType(); //dummy initialization. Not used in builder and solver.
    TSystemVectorType mDx = TSystemVectorType();
    TSystemVectorType mb = TSystemVectorType();

    //initial operations ... things that are constant over the Solution Step
    pScheme->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

    // if(BaseType::mRebuildLevel > 0)
    // {
    //   pBuilderAndSolver->BuildLHS(pScheme, r_model_part, mA); //calculate Mass Matrix
    // }
    pBuilderAndSolver->BuildLHS(pScheme, r_model_part, mA); //calculate Mass Matrix

    mSolutionStepIsInitialized = true;

    KRATOS_CATCH("")
  }

  //**********************************************************************
  //**********************************************************************
  /*
      SOLUTION OF THE PROBLEM OF INTEREST
    */
  //**********************************************************************

  double Solve() override
  {
    KRATOS_TRY

    // std::cout<<" Solve() in explicit two step v p"<<std::endl;

    //pointers needed in the solution
    typename TSchemeType::Pointer pScheme = GetScheme();
    typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
    ModelPart &r_model_part = BaseType::GetModelPart();

    DofsArrayType rDofSet; //dummy initialization. Not used in builder and solver
    TSystemMatrixType mA = TSystemMatrixType();
    TSystemVectorType mDx = TSystemVectorType();
    TSystemVectorType mBmomentum = TSystemVectorType();
    TSystemVectorType mBcontinuity = TSystemVectorType();

    //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
    //if the operations needed were already performed this does nothing
    if (mInitializeWasPerformed == false)
      Initialize();

    //prints informations about the current time
    if (this->GetEchoLevel() == 2 && r_model_part.GetCommunicator().MyPID() == 0)
    {
      std::cout << " " << std::endl;
      std::cout << "CurrentTime = " << r_model_part.GetProcessInfo()[TIME] << std::endl;
    }

    ProcessInfo rCurrentProcessInfo = r_model_part.GetProcessInfo();
    //double currentTime   = rCurrentProcessInfo[TIME];
    //int step   = rCurrentProcessInfo[STEP];
    //double timeStep     = rCurrentProcessInfo[DELTA_TIME];

    //initialize solution step
    // if(mSolutionStepIsInitialized == false)
    // InitializeSolutionStep();

    // it predicts the timeStep, it clear the RHS of momentum and continuity
    pScheme->InitializeSolutionStep(r_model_part, mA, mDx, mBmomentum);

    ////////////////////// starting momentum step solution /////////////////////

    pBuilderAndSolver->Build(pScheme, r_model_part, mA, mBmomentum);

    // if(step==2)
    // InitializeDensity(r_model_part);

    pBuilderAndSolver->BuildRHS(pScheme, r_model_part, mBcontinuity);

    pScheme->Update(r_model_part, rDofSet, mA, mDx, mBmomentum); // Explicitly integrates the equation of motion.

    //Finalisation of the solution step,
    //operations to be done after achieving convergence, for example the
    //Final Residual Vector (mb) has to be saved in there
    //to avoid error accumulation
    pScheme->FinalizeSolutionStep(r_model_part, mA, mDx, mBmomentum);

    //move the mesh if needed
    BaseType::MoveMesh();

    //Cleaning memory after the solution
    pScheme->Clean();

    //reset flags for next step
    mSolutionStepIsInitialized = false;

    ////////////////////// momentum step solution finished /////////////////////

    // pBuilderAndSolver->BuildRHS(pScheme, r_model_part, mBcontinuity);

    pBuilderAndSolver->BuildLHS(pScheme, r_model_part, mA); //calculate Mass Matrix

    // pBuilderAndSolver->Build(pScheme, r_model_part, mA, mb);

    // this->ComputeAndUpdatePressureFromDensity(r_model_part,rCurrentProcessInfo);
    // if(timeStep<currentTime)

    this->UpdateDensity(r_model_part, rCurrentProcessInfo);  // it computes the new density and saves the previous ones in PRESSURE,1
    this->UpdatePressure(r_model_part, rCurrentProcessInfo); // it computes the new pressures

    pScheme->FinalizeSolutionStep(r_model_part, mA, mDx, mBcontinuity);

    pScheme->Clean();

    return 0.00;

    KRATOS_CATCH("")
  }

  //   //**********************************************************************
  //   //**********************************************************************
  //   /*
  //                     SOLUTION OF THE PROBLEM OF INTEREST
  //    */
  //   //**********************************************************************

  // double Solve()
  // {
  //   KRATOS_TRY

  //     // std::cout<<" Solve() in explicit two step v p"<<std::endl;

  //   //pointers needed in the solution
  //   typename TSchemeType::Pointer pScheme = GetScheme();
  //   typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
  //   ModelPart& r_model_part                                   = BaseType::GetModelPart();

  //   DofsArrayType rDofSet; //dummy initialization. Not used in builder and solver
  //   TSystemMatrixType mA  = TSystemMatrixType();
  //   TSystemVectorType mDx = TSystemVectorType();
  //   TSystemVectorType mb  = TSystemVectorType();

  //   //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
  //   //if the operations needed were already performed this does nothing
  //   if(mInitializeWasPerformed == false)
  //     Initialize();

  //   //prints informations about the current time
  //   if (this->GetEchoLevel() == 2 && r_model_part.GetCommunicator().MyPID() == 0 )
  //     {
  // 	std::cout << " " << std::endl;
  // 	std::cout << "CurrentTime = " << r_model_part.GetProcessInfo()[TIME] << std::endl;
  //     }

  //   //initialize solution step
  //   // if(mSolutionStepIsInitialized == false)
  //   // InitializeSolutionStep();

  //   // it predicts the timeStep, it clear the RHS of momentum and continuity
  //   pScheme->InitializeSolutionStep(r_model_part, mA, mDx, mb);

  //   ////////////////////// starting momentum step solution /////////////////////

  //   pBuilderAndSolver->BuildLHS(pScheme, r_model_part, mA); //calculate Mass Matrix

  //   pBuilderAndSolver->BuildRHS(pScheme, r_model_part, mb);

  //   pScheme->Update(r_model_part, rDofSet, mA, mDx, mb); // Explicitly integrates the equation of motion.

  //   //Finalisation of the solution step,
  //   //operations to be done after achieving convergence, for example the
  //   //Final Residual Vector (mb) has to be saved in there
  //   //to avoid error accumulation
  //   pScheme->FinalizeSolutionStep(r_model_part, mA, mDx, mb);

  //   //move the mesh if needed
  //   BaseType::MoveMesh();

  //   //Cleaning memory after the solution
  //   pScheme->Clean();

  //   //reset flags for next step
  //   mSolutionStepIsInitialized = false;

  //   ////////////////////// momentum step solution finished /////////////////////

  //   pBuilderAndSolver->Build(pScheme, r_model_part, mA, mb);

  //   this->UpdatePressure(r_model_part); // Explicitly integrates the equation of motion.

  //   pScheme->FinalizeSolutionStep(r_model_part, mA, mDx, mb);

  //   pScheme->Clean();

  //   return 0.00;

  //   KRATOS_CATCH( "" )

  //     }

  void InitializeDensity(ModelPart &rModelPart)
  {
    KRATOS_TRY

    std::cout << "InitializeDensity() InitializeDensity() InitializeDensity() InitializeDensity()" << std::endl;
    const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

    OpenMPUtils::PartitionVector NodePartition;
    OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

    const int nnodes = static_cast<int>(rModelPart.Nodes().size());
    NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
    for (int i = 0; i < nnodes; i++)
    {
      NodesArrayType::iterator itNode = NodeBegin + i;

      // Current step information "N+1" (before step update).
      // double& previous_density  = (itNode)->FastGetSolutionStepValue(PRESSURE,1);
      double &previous_density = (itNode)->FastGetSolutionStepValue(PRESSURE, 0);
      previous_density = 1000;
    }

    KRATOS_CATCH("")
  }

  void ComputeAndUpdatePressureFromDensity(ModelPart &rModelPart, ProcessInfo rCurrentProcessInfo)
  {
    KRATOS_TRY

    const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

    OpenMPUtils::PartitionVector NodePartition;
    OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

    const int nnodes = static_cast<int>(rModelPart.Nodes().size());
    NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

    // double currentTime   = rCurrentProcessInfo[TIME];  //the first step is (time = initial_time + delta time )
    // double timeStep     = rCurrentProcessInfo[DELTA_TIME];

#pragma omp parallel for firstprivate(NodeBegin)
    for (int i = 0; i < nnodes; i++)
    {
      NodesArrayType::iterator itNode = NodeBegin + i;

      // Current step information "N+1" (before step update).
      const double &bulk_term = (itNode)->FastGetSolutionStepValue(NODAL_MASS);
      double &density_rhs = (itNode)->FastGetSolutionStepValue(NODAL_ERROR);
      double current_density = (itNode)->FastGetSolutionStepValue(PRESSURE, 0);
      double &current_pressure = (itNode)->FastGetSolutionStepValue(PRESSURE, 0);
      // double& previous_density  = (itNode)->FastGetSolutionStepValue(PRESSURE,1);

      if (bulk_term != 0)
      {
        // Solution of the explicit equation:
        if ((itNode)->IsFixed(PRESSURE) == false && (itNode)->IsNot(ISOLATED))
        {
          // (itNode)->FastGetSolutionStepValue(PRESSURE,1)=current_density;
          current_density = density_rhs / bulk_term;
          double bulkModulus = 2100000000;
          double initialDensity = 1000;
          double gammaExponent = 7.0;
          double densityRatio = current_density / initialDensity;
          double powerDensityRatio = pow(densityRatio, gammaExponent);
          current_pressure = -bulkModulus * (powerDensityRatio - 1.0);
          // std::cout<<currentTime<<" s) -> current_density "<<current_density<<" density_rhs= "<<density_rhs<<"   bulk_term "<<bulk_term<<std::endl;
        }
        else if ((itNode)->Is(ISOLATED))
        {
          std::cout << "ISOLATED NODE " << (itNode)->X() << " " << (itNode)->Y() << std::endl;
          current_pressure = 0;
        }
      }
    }

    KRATOS_CATCH("")
  }

  void UpdateDensity(ModelPart &rModelPart, ProcessInfo rCurrentProcessInfo)
  {
    KRATOS_TRY

    const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

    OpenMPUtils::PartitionVector NodePartition;
    OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

    const int nnodes = static_cast<int>(rModelPart.Nodes().size());
    NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
    for (int i = 0; i < nnodes; i++)
    {
      NodesArrayType::iterator itNode = NodeBegin + i;

      // Current step information "N+1" (before step update).
      const double &bulk_term = (itNode)->FastGetSolutionStepValue(NODAL_MASS);
      double &density_rhs = (itNode)->FastGetSolutionStepValue(NODAL_ERROR);
      double &current_density = (itNode)->FastGetSolutionStepValue(DENSITY, 0);
      // double& previous_density  = (itNode)->FastGetSolutionStepValue(PRESSURE,1);

      if (bulk_term != 0)
      {
        // Solution of the explicit equation:
        // (itNode)->FastGetSolutionStepValue(PRESSURE,1)=current_density;
        current_density = density_rhs / bulk_term;
        // std::cout<<currentTime<<" s) -> current_density "<<current_density<<" density_rhs= "<<density_rhs<<"   bulk_term "<<bulk_term<<std::endl;
      }
    }

    KRATOS_CATCH("")
  }

  void UpdatePressure(ModelPart &rModelPart, ProcessInfo rCurrentProcessInfo)
  {
    KRATOS_TRY

    const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

    OpenMPUtils::PartitionVector NodePartition;
    OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

    const int nnodes = static_cast<int>(rModelPart.Nodes().size());
    NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

    // double currentTime   = rCurrentProcessInfo[TIME];  //the first step is (time = initial_time + delta time )
    // double timeStep     = rCurrentProcessInfo[DELTA_TIME];

#pragma omp parallel for firstprivate(NodeBegin)
    for (int i = 0; i < nnodes; i++)
    {
      NodesArrayType::iterator itNode = NodeBegin + i;

      // Current step information "N+1" (before step update).
      double &current_pressure = (itNode)->FastGetSolutionStepValue(PRESSURE, 0);
      double &current_density = (itNode)->FastGetSolutionStepValue(DENSITY, 0);
      // double current_density=current_pressure;
      // double previous_density  = (itNode)->FastGetSolutionStepValue(PRESSURE,1);

      // Solution of the explicit equation:
      if ((itNode)->IsFixed(PRESSURE) == false && (itNode)->IsNot(ISOLATED))
      {
        double bulkModulus = 2100000000;
        double initialDensity = 1000;
        double gammaExponent = 7.0;
        double densityRatio = current_density / initialDensity;
        double powerDensityRatio = pow(densityRatio, gammaExponent);

        current_pressure = -bulkModulus * (powerDensityRatio - 1.0);

        // std::cout<<currentTime<<" s) -> previous_density "<<previous_density<<" current_density "<<current_density<<" densityRatio="<<densityRatio<<"   powerDensityRatio="<<powerDensityRatio<<"    current_pressure  "<<current_pressure<<std::endl;
      }

      if ((itNode)->Is(ISOLATED))
      {
        // std::cout<<"ISOLATED NODE "<<(itNode)->X()<<" "<<(itNode)->Y()<<std::endl;
        current_pressure = 0;
      }
    }

    KRATOS_CATCH("")
  }

  //   void UpdatePressure(ModelPart& rModelPart)
  //   {
  //     KRATOS_TRY

  //       const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

  //     OpenMPUtils::PartitionVector NodePartition;
  //     OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(),NumThreads,NodePartition);

  //     const int nnodes = static_cast<int>(rModelPart.Nodes().size());
  //     NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

  // #pragma omp parallel for firstprivate(NodeBegin)
  //     for(int i = 0;  i < nnodes; i++)
  //       {
  // 	NodesArrayType::iterator itNode = NodeBegin + i;

  // 	// Current step information "N+1" (before step update).
  // 	const double& bulk_term  = (itNode)->FastGetSolutionStepValue(NODAL_MASS);
  // 	double& pressure_rhs  = (itNode)->FastGetSolutionStepValue(NODAL_ERROR);
  // 	double& current_pressure  = (itNode)->FastGetSolutionStepValue(PRESSURE,0);
  // 	double& previous_pressure  = (itNode)->FastGetSolutionStepValue(PRESSURE,1);

  // 	if(bulk_term!=0){
  // 	  // Solution of the explicit equation:
  // 	  if((itNode)->IsFixed(PRESSURE) == false && (itNode)->IsNot(ISOLATED)){
  // 	    previous_pressure=current_pressure;
  // 	    current_pressure =  pressure_rhs/bulk_term;
  // 	  }
  // 	}
  // 	if((itNode)->Is(ISOLATED)){
  // 	  // std::cout<<"ISOLATED NODE "<<(itNode)->X()<<" "<<(itNode)->Y()<<std::endl;
  // 	    current_pressure  =0;
  // 	    previous_pressure  =0;
  // 	}
  //       }

  //     KRATOS_CATCH("")
  //       }

  //**********************************************************************
  //**********************************************************************

  void Clear() override
  {
    KRATOS_TRY
    std::cout << "Explicit strategy Clear function used" << std::endl;

    //setting to zero the internal flag to ensure that the dof sets are recalculated
    //GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
    //GetBuilderAndSolver()->Clear();

    GetScheme()->Clear();

    KRATOS_CATCH("")
  }

  /*@} */
  /**@name Operators
     */
  /*@{ */

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

private:
  /**@name Protected static Member Variables */
  /*@{ */

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

  /*@} */

protected:
  /**@name Static Member Variables */
  /*@{ */

  /*@} */
  /**@name Member Variables */
  /*@{ */

  typename TSchemeType::Pointer mpScheme;

  typename TLinearSolver::Pointer mpLinearSolver;

  typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

  TSystemVectorPointerType mpDx;
  TSystemVectorPointerType mpb;
  TSystemMatrixPointerType mpA;

  /**
       Flag telling if it is needed to reform the DofSet at each
       solution step or if it is possible to form it just once
       - true  => reforme at each time step
       - false => form just one (more efficient)

       Default = false
    */
  bool mReformDofSetAtEachStep;

  /**
       Flag telling if it is needed or not to compute the reactions

       default = true
    */
  bool mCalculateReactionsFlag;

  bool mSolutionStepIsInitialized;

  bool mInitializeWasPerformed;

  bool mComputeTime;

  /*@} */
  /**@name Private Operators*/
  /*@{ */

  //**********************************************************************
  //**********************************************************************

  void CalculateReactions()
  {
  }

  //**********************************************************************
  //**********************************************************************

  /**
     * function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */

  int Check() override
  {
    KRATOS_TRY

    BaseType::Check();

    GetScheme()->Check(BaseType::GetModelPart());

    return 0;

    KRATOS_CATCH("")
  }

  //***************************************************************************
  //***************************************************************************

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

  /** Copy constructor.
     */
  ExplicitTwoStepVPStrategy(const ExplicitTwoStepVPStrategy &Other){};

  /*@} */

}; /* Class ExplicitTwoStepVPStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_TWO_STEP_V_P_STRATEGY  defined */
