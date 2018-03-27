//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               March 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_NEWTON_RAPHSON_STRATEGY_H_INCLUDED)
#define KRATOS_NEWTON_RAPHSON_STRATEGY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_strategies/strategies/solution_strategy.hpp"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
//default convergence criterion
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

#include "solid_mechanics_application_variables.h"

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

/**   
 * @class NewtonRaphsonStrategy
 * @brief This is the base Newton Raphson strategy
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is surpassed) using a Newton Raphson algorithm
 */

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver // = LinearSolver<TSparseSpace,TDenseSpace>
          >
class NewtonRaphsonStrategy : public SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
 public:
  ///@name Type Definitions
  ///@{

  // Counted pointer of ClassName
  KRATOS_CLASS_POINTER_DEFINITION(NewtonRaphsonStrategy);

  typedef SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>            BaseType;

  typedef typename BaseType::LocalFlagType                                 LocalFlagType;
  
  typedef ConvergenceCriteria<TSparseSpace, TDenseSpace>        ConvergenceCriterionType;

  typedef typename BaseType::BuilderAndSolverType                   BuilderAndSolverType;

  typedef typename BaseType::SchemeType                                       SchemeType;

  typedef TLinearSolver                                                 LinearSolverType;

  typedef TSparseSpace                                                   SparseSpaceType;

  typedef typename BaseType::DofsArrayType                                 DofsArrayType;

  typedef typename BaseType::SystemMatrixType                           SystemMatrixType;

  typedef typename BaseType::SystemVectorType                           SystemVectorType;

  typedef typename BaseType::SystemMatrixPointerType             SystemMatrixPointerType;

  typedef typename BaseType::SystemVectorPointerType             SystemVectorPointerType;

  ///@}
  ///@name Life Cycle

  ///@{

  /**
   * Default constructor 
   * @param rModelPart The model part of the problem
   * @param pScheme The integration scheme
   * @param pBuilderAndSolver The builder and solver employed
   * @param pConvergenceCriteria The convergence criteria employed
   * @param rOptions The solution options
   * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
   */  
  NewtonRaphsonStrategy(ModelPart& rModelPart,
                        typename SchemeType::Pointer pScheme,
                        typename BuilderAndSolverType::Pointer pBuilderAndSolver,
                        typename ConvergenceCriterionType::Pointer pConvergenceCriterion,
                        Flags& rOptions,
                        unsigned int MaxIterations = 30)
      : SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, rOptions)
  {
    KRATOS_TRY
 
    // Saving the convergence criteria to be used
    mpConvergenceCriteria = pConvergenceCriterion;

    // Saving the scheme
    mpScheme = pScheme;

    // Setting up the default builder and solver
    mpBuilderAndSolver = pBuilderAndSolver;
       
    // Maximum iterations allowed
    mMaxIterationNumber = MaxIterations;

    // Set Options
    this->mOptions = rOptions;

    // Tells to the builder and solver if the reactions have to be Calculated or not
    if(this->mOptions.Is(LocalFlagType::COMPUTE_REACTIONS))
      mpBuilderAndSolver->SetCalculateReactionsFlag(true);
    else
      mpBuilderAndSolver->SetCalculateReactionsFlag(false);
    
    // Tells to the Builder And Solver if the system matrix and vectors need to be reshaped at each step or not
    if(this->mOptions.Is(LocalFlagType::REFORM_DOFS))
      mpBuilderAndSolver->SetReshapeMatrixFlag(true);
    else
      mpBuilderAndSolver->SetReshapeMatrixFlag(false);

    mpA  = TSparseSpace::CreateEmptyMatrixPointer();
    mpDx = TSparseSpace::CreateEmptyVectorPointer();
    mpb  = TSparseSpace::CreateEmptyVectorPointer();

    KRATOS_CATCH("")
  }
  

  /**
   * Default constructor 
   * @param rModelPart The model part of the problem
   * @param pScheme The integration scheme
   * @param pLinearSolver The linear solver employed
   * @param pConvergenceCriteria The convergence criteria employed
   * @param rOptions The solution options
   * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
   */
  NewtonRaphsonStrategy(ModelPart& rModelPart,
                        typename SchemeType::Pointer pScheme,
                        typename LinearSolverType::Pointer pLinearSolver,
                        typename ConvergenceCriterionType::Pointer pConvergenceCriterion,
                        Flags& rOptions,
                        unsigned int MaxIterations = 30)
      : NewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, typename BuilderAndSolverType::Pointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pLinearSolver)), pConvergenceCriterion, rOptions, MaxIterations)
  {}
    

  /** 
   * @brief Destructor.
   * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
   */
  ~NewtonRaphsonStrategy() override
  {
    Clear();
  }

  ///@}
  ///@name Operators

  ///@{

  ///@}
  ///@name Operations
  ///@{
  
  /**
   * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
   */
  void InitializeSolutionStep() override
  {
    KRATOS_TRY
                    
    //initialize system elements, conditions..
    if(this->IsNot(LocalFlagType::INITIALIZED))
      this->Initialize();

    //set implex
    if(this->mOptions.Is(LocalFlagType::IMPLEX) && this->GetModelPart().GetProcessInfo().Has(IMPLEX))
      this->GetModelPart().GetProcessInfo()[IMPLEX] = 1;
    
    //prints informations about the current time
    //KRATOS_INFO << "  [STEP:" << this->GetModelPart().GetProcessInfo()[STEP] << "  TIME: "<< this->GetModelPart().GetProcessInfo()[TIME]<< "]\n" << LoggerMessage::Category::STATUS;
    
    //set up the system
    if(this->mOptions.Is(LocalFlagType::REFORM_DOFS)) 
      this->SetSystemDofs();   

    //setting up the vectors involved to the correct size
    //double begin_time = OpenMPUtils::GetCurrentTime();
    mpBuilderAndSolver->ResizeAndInitializeVectors(mpScheme, mpA, mpDx, mpb,
                                                   this->GetModelPart().Elements(),
                                                   this->GetModelPart().Conditions(),
                                                   this->GetModelPart().GetProcessInfo());
    //double end_time = OpenMPUtils::GetCurrentTime();

    //KRATOS_INFO("system_resize_time") << ": system_resize_time : " << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;
    
    //initial operations ... things that are constant over the Solution Step
    mpBuilderAndSolver->InitializeSolutionStep(this->GetModelPart(), (*mpA), (*mpDx), (*mpb));
    
    //initial operations ... things that are constant over the Solution Step
    mpScheme->InitializeSolutionStep(this->GetModelPart(), (*mpA), (*mpDx), (*mpb));

    this->Predict();
        
    KRATOS_CATCH("")
  }

  /**
   * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
   */
  void FinalizeSolutionStep() override
  {
    KRATOS_TRY

    //Finalization of the solution step, operations to be done after achieving convergence

    //set implex calculation
    if(this->mOptions.Is(LocalFlagType::IMPLEX) && this->GetModelPart().GetProcessInfo().Has(IMPLEX))
      this->GetModelPart().GetProcessInfo()[IMPLEX] = 0;
    
    //calculate reactions if required
    if(this->mOptions.Is(LocalFlagType::COMPUTE_REACTIONS) || this->mOptions.Is(LocalFlagType::IMPLEX))
      mpBuilderAndSolver->CalculateReactions(mpScheme, this->GetModelPart(), (*mpA), (*mpDx), (*mpb));
       
    //finalize scheme anb builder and solver
    mpScheme->FinalizeSolutionStep(this->GetModelPart(), (*mpA), (*mpDx), (*mpb));
    mpBuilderAndSolver->FinalizeSolutionStep(this->GetModelPart(), (*mpA), (*mpDx), (*mpb));

    //cleaning memory after the solution
    mpScheme->Clean();

    if(this->mOptions.Is(LocalFlagType::REFORM_DOFS)){
      //deallocate the systemvectors
      SparseSpaceType::Clear(mpA);
      SparseSpaceType::Clear(mpDx);
      SparseSpaceType::Clear(mpb);

      this->Clear();
    }

    //this->Finalize();
    
    KRATOS_CATCH("")
  }

  
  /**
   * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
   */   
  bool SolveSolutionStep() override
  {
    KRATOS_TRY

    //initializing the parameters of the Newton-Raphson cicle
    unsigned int iteration_number = 1;
    
    this->Set(LocalFlagType::CONVERGED, this->SolveIteration(iteration_number));
   
    //Iteration Cicle... performed only for NonLinearProblems
    while( this->IsNot(LocalFlagType::CONVERGED) && iteration_number++ < mMaxIterationNumber)
    {
      this->Set(LocalFlagType::CONVERGED, this->SolveIteration(iteration_number));
    }

    //plots a warning if the maximum number of iterations is exceeded
    if(iteration_number >= mMaxIterationNumber)
    {
      KRATOS_WARNING("Max Iterations Exceeded") << " **** Maximum iterations Exceeded [" << iteration_number << "] ****\n"; 
    }
    
    return (this->Is(LocalFlagType::CONVERGED));

    KRATOS_CATCH("")

  }

  
  /**
   * @brief Clears the internal storage
   */
  void Clear() override
  {
    KRATOS_TRY

    // if the preconditioner is saved between solves, it should be cleared here.
    mpBuilderAndSolver->GetLinearSystemSolver()->Clear();

    if(mpA != nullptr)
      SparseSpaceType::Clear(mpA);
    if(mpDx != nullptr)
      SparseSpaceType::Clear(mpDx);
    if(mpb != nullptr)
      SparseSpaceType::Clear(mpb);

    //setting to zero the internal flag to ensure that the dof sets are recalculated
    mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);
    mpBuilderAndSolver->Clear();
    mpScheme->Clear();

    KRATOS_CATCH("")
  }

  /**
   * @brief Function to perform expensive checks.
   * @details It is designed to be called ONCE to verify that the input is correct.
   */
  int Check() override
  {
    KRATOS_TRY

    //check the model part    
    BaseType::Check();

    //check the scheme
    mpScheme->Check(this->GetModelPart());

    //check the builder and solver
    mpBuilderAndSolver->Check(this->GetModelPart());

    //check the convergence criterion
    mpConvergenceCriteria->Check(this->GetModelPart());

    return 0;

    KRATOS_CATCH("")
  }
  

  ///@}
  ///@name Access

  ///@{

  /**
   * @brief Set method for the time scheme
   * @param pScheme The pointer to the time scheme considered
   */
  void SetScheme(typename SchemeType::Pointer pScheme)
  {
    mpScheme = pScheme;
  };

  /**
   * @brief Get method for the time scheme
   * @return mpScheme: The pointer to the time scheme considered
   */
  typename SchemeType::Pointer GetScheme()
  {
    return mpScheme;
  };

  /**
   * @brief Set method for the builder and solver
   * @param pNewBuilderAndSolver The pointer to the builder and solver considered
   */
  void SetBuilderAndSolver(typename BuilderAndSolverType::Pointer pBuilderAndSolver)
  {
    mpBuilderAndSolver = pBuilderAndSolver;
  };

  /**
   * @brief Get method for the builder and solver
   * @return mpBuilderAndSolver: The pointer to the builder and solver considered
   */
  typename BuilderAndSolverType::Pointer GetBuilderAndSolver()
  {
    return mpBuilderAndSolver;
  };


  /**
   * @brief This method sets the flag mMaxIterationNumber
   * @param MaxIterationNumber This is the maximum number of on linear iterations
   */
  void SetMaxIterationNumber(unsigned int MaxIterationNumber)
  {
    mMaxIterationNumber = MaxIterationNumber;
  }

  /**
   * @brief This method gets the flag mMaxIterationNumber
   * @return mMaxIterationNumber: This is the maximum number of on linear iterations
   */
  unsigned int GetMaxIterationNumber()
  {
    return mMaxIterationNumber;
  }
  
  /**
   * @brief This method returns the LHS matrix
   * @return The LHS matrix
   */
  SystemMatrixType &GetSystemMatrix()
  {
    SystemMatrixType &mA = *mpA;

    return mA;
  }

  /**
   * @brief This method directly sets the input as the LHS
   * @param A The LHS matrix
   */
  void GetDirectSystemMatrix(SystemMatrixType& A)
  {
    A = *mpA;
  }
    

  ///@}
  ///@name Inquiry
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
  
  typename SchemeType::Pointer mpScheme; /// The pointer to the time scheme employed
  typename BuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver employed
  typename ConvergenceCriterionType::Pointer mpConvergenceCriteria; /// The pointer to the convergence criteria employed

  SystemVectorPointerType mpDx; /// The incremement in the solution
  SystemVectorPointerType mpb; /// The RHS vector of the system of equations
  SystemMatrixPointerType mpA; /// The LHS matrix of the system of equations

  unsigned int mMaxIterationNumber; /// The maximum number of iterations, 30 by default

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{
    /**
   * @brief Initialization of member variables and prior operations
   */    
  void Initialize() override
  {
    KRATOS_TRY

    //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
    if(mpScheme->SchemeIsInitialized() == false)
      mpScheme->Initialize(this->GetModelPart());

    //Initialize The Elements - OPERATIONS TO BE DONE ONCE
    if(mpScheme->ElementsAreInitialized() == false)
      mpScheme->InitializeElements(this->GetModelPart());
      
    //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
    if(mpScheme->ConditionsAreInitialized() == false)
      mpScheme->InitializeConditions(this->GetModelPart());
    
    //initialisation of the convergence criteria
    if(mpConvergenceCriteria->mConvergenceCriteriaIsInitialized == false)
      mpConvergenceCriteria->Initialize(this->GetModelPart());
    
    //set up the system
    this->SetSystemDofs();

    this->Set(LocalFlagType::INITIALIZED,true);
    
    KRATOS_CATCH("")
  }

  
  /**
   * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
   values of the solution step of interest are assumed equal to the old values
  */
  void Predict() override
  {
    KRATOS_TRY

    mpScheme->Predict(this->GetModelPart(), mpBuilderAndSolver->GetDofSet(), (*mpA), (*mpDx), (*mpb));

    //move the mesh if needed
    if(this->mOptions.Is(LocalFlagType::MOVE_MESH))
      BaseType::MoveMesh();
    
    KRATOS_CATCH("")
  }


  /**
   * @brief Here the database is updated
   * @param A The LHS matrix of the system of equations
   * @param Dx The incremement in the solution
   * @param b The RHS vector of the system of equations
   * @param MoveMesh The flag that allows to move the mesh
   */
  void Update() override
  {
    KRATOS_TRY
           
    mpScheme->Update(this->GetModelPart(), mpBuilderAndSolver->GetDofSet(), (*mpA), (*mpDx), (*mpb));

    // Move the mesh if needed
    if(this->mOptions.Is(LocalFlagType::MOVE_MESH))
      BaseType::MoveMesh();

    KRATOS_CATCH("")
  }
 


  /**
   * @brief Performs all the required operations to reform dofs
   */  
  void SetSystemDofs()
  {
    KRATOS_TRY


    KRATOS_INFO(" Reform Dofs ") << " Flag = " <<this->mOptions.Is(LocalFlagType::REFORM_DOFS) << std::endl;
                                                                                        
    //set up the system, operation performed just once unless it is required to reform the dof set at each iteration
        
    //setting up the list of the DOFs to be solved
    //double begin_time = OpenMPUtils::GetCurrentTime();
    mpBuilderAndSolver->SetUpDofSet(mpScheme, this->GetModelPart());
    //double end_time = OpenMPUtils::GetCurrentTime();
    //KRATOS_INFO("setup_dofs_time") << "setup_dofs_time : " << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;
      
    //shaping correctly the system
    //begin_time = OpenMPUtils::GetCurrentTime();
    mpBuilderAndSolver->SetUpSystem(this->GetModelPart());
    //end_time = OpenMPUtils::GetCurrentTime();
    //KRATOS_INFO("setup_system_time") << ": setup_system_time : " << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;

    KRATOS_CATCH("")
  }


  /**
   * @brief Solves the iteration. This function returns true if a solution has been found, false otherwise.
   */   
  bool SolveIteration(unsigned int& rIterationNumber)
  {
    KRATOS_TRY
        
    bool is_converged = false;
   
    // Warning info
    if(!SparseSpaceType::Size((*mpDx)))
      KRATOS_WARNING("DOFS") << "solution has zero size, no free DOFs" << std::endl;

    // Initialize Iteration
    mpScheme->InitializeNonLinIteration(this->GetModelPart(), (*mpA), (*mpDx), (*mpb));
    
    //setting the number of iteration
    this->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = rIterationNumber;
    
    is_converged = mpConvergenceCriteria->PreCriteria(this->GetModelPart(), mpBuilderAndSolver->GetDofSet(), (*mpA), (*mpDx), (*mpb));

    //function to perform the building and the solving phase.
    if(this->mOptions.IsNot(LocalFlagType::CONSTANT_SYSTEM_LHS)){

      TSparseSpace::SetToZero((*mpA));
      TSparseSpace::SetToZero((*mpDx));
      TSparseSpace::SetToZero((*mpb));
      
      mpBuilderAndSolver->BuildAndSolve(mpScheme, this->GetModelPart(), (*mpA), (*mpDx), (*mpb));      
    }
    else{
      
      TSparseSpace::SetToZero((*mpDx));
      TSparseSpace::SetToZero((*mpb));
      
      mpBuilderAndSolver->BuildRHSAndSolve(mpScheme, this->GetModelPart(), (*mpA), (*mpDx), (*mpb));
    }

    // Debugging info
    KRATOS_TRACE("LHS") << "LHS = " << (*mpA) <<  "\n" << LoggerMessage::Category::CHECKING;
    KRATOS_TRACE("Dx")  << "Solution = " << (*mpDx) <<  "\n" << LoggerMessage::Category::CHECKING;
    KRATOS_TRACE("RHS") << "RHS  = " << (*mpb) << "\n" << LoggerMessage::Category::CHECKING;

    // Updating the results
    this->Update();

    // Finalize Iteration
    mpScheme->FinalizeNonLinIteration(this->GetModelPart(), (*mpA), (*mpDx), (*mpb));

    if(is_converged == true)
    {

      //initialisation of the convergence criteria (after first calculation only)
      if( rIterationNumber == 1 )
        mpConvergenceCriteria->InitializeSolutionStep(this->GetModelPart(), mpBuilderAndSolver->GetDofSet(), (*mpA), (*mpDx), (*mpb));

      
      if(mpConvergenceCriteria->GetActualizeRHSflag() == true)
      {
        TSparseSpace::SetToZero((*mpb));

        mpBuilderAndSolver->BuildRHS(mpScheme, this->GetModelPart(), (*mpb));
      }

      is_converged = mpConvergenceCriteria->PostCriteria(this->GetModelPart(), mpBuilderAndSolver->GetDofSet(), (*mpA), (*mpDx), (*mpb));
    }   
    
    return is_converged;

    KRATOS_CATCH("")
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
  ///@name Private Inquiry
  ///@{
  ///@}
  ///@name Un accessible methods
  ///@{

  /// Copy constructor.
  NewtonRaphsonStrategy(const NewtonRaphsonStrategy &Other){};

  ///@}

}; /// Class NewtonRaphsonStrategy

///@}

///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

///@} addtogroup block
  
}  // namespace Kratos.
#endif // KRATOS_NEWTON_RAPHSON_STRATEGY_H_INCLUDED defined
