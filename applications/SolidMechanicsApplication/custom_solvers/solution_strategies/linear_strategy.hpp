//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               March 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LINEAR_STRATEGY_H_INCLUDED)
#define KRATOS_LINEAR_STRATEGY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_strategies/solution_strategy.hpp"

//default builder and solver
#include "custom_solvers/solution_builders_and_solvers/block_builder_and_solver.hpp"

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
 * @class LinearStrategy
 * @brief This is the base linear strategy jacobi / gauss-seidel linear strategies
 */

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver // = LinearSolver<TSparseSpace,TDenseSpace>
          >
class LinearStrategy : public SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
 public:
  ///@name Type Definitions
  ///@{

  // Counted pointer of ClassName
  KRATOS_CLASS_POINTER_DEFINITION(LinearStrategy);

  typedef SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>            BaseType;

  typedef typename BaseType::LocalFlagType                                 LocalFlagType;

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
   * @param rOptions The solution options
   */
  LinearStrategy(ModelPart& rModelPart,
                 typename SchemeType::Pointer pScheme,
                 typename BuilderAndSolverType::Pointer pBuilderAndSolver,
                 Flags& rOptions)
      : SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, rOptions)
  {
    KRATOS_TRY

    // Saving the scheme
    mpScheme = pScheme;

    // Setting up the default builder and solver
    mpBuilderAndSolver = pBuilderAndSolver;

    // Set Options
    this->mOptions = rOptions;

    // Tells to the builder and solver if the reactions have to be Calculated or not
    mpBuilderAndSolver->GetOptions().Set(LocalFlagType::COMPUTE_REACTIONS, this->mOptions.Is(LocalFlagType::COMPUTE_REACTIONS));

    // Tells to the Builder And Solver if the system matrix and vectors need to be reshaped at each step or not
    mpBuilderAndSolver->GetOptions().Set(LocalFlagType::REFORM_DOFS, this->mOptions.Is(LocalFlagType::REFORM_DOFS));

    mpBuilderAndSolver->SetEchoLevel(this->mEchoLevel);

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
   * @param rOptions The solution options
   */
  LinearStrategy(ModelPart& rModelPart,
                 typename SchemeType::Pointer pScheme,
                 typename LinearSolverType::Pointer pLinearSolver,
                 Flags& rOptions)
      : LinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, Kratos::make_shared<BlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(pLinearSolver), rOptions)
  {}


  /**
   * @brief Destructor.
   * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
   */
  ~LinearStrategy() override
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

    //prints informations about the current time
    //KRATOS_INFO("") << "  [STEP:" << this->GetModelPart().GetProcessInfo()[STEP] << "  TIME: "<< this->GetModelPart().GetProcessInfo()[TIME]<< "]\n" << LoggerMessage::Category::STATUS;

    //set up the system
    if(this->mOptions.Is(LocalFlagType::REFORM_DOFS))
      this->SetSystemDofs();

    //setting up the vectors involved to the correct size
    double begin_time = OpenMPUtils::GetCurrentTime();
    mpBuilderAndSolver->SetUpSystemMatrices(mpScheme, this->GetModelPart(), mpA, mpDx, mpb);
    double end_time = OpenMPUtils::GetCurrentTime();

    if (this->mEchoLevel >= 2)
      KRATOS_INFO("system_resize_time") << ": system_resize_time : " << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;

    //initial operations ... things that are constant over the Solution Step
    mpBuilderAndSolver->InitializeSolutionStep(mpScheme, this->GetModelPart(), mpA, mpDx, mpb);

    //initial operations ... things that are constant over the Solution Step
    mpScheme->InitializeSolutionStep(this->GetModelPart());

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

    //calculate reactions if required
    if(this->mOptions.Is(LocalFlagType::COMPUTE_REACTIONS))
      mpBuilderAndSolver->CalculateReactions(mpScheme, this->GetModelPart(), (*mpA), (*mpDx), (*mpb));

    //finalize scheme anb builder and solver
    mpScheme->FinalizeSolutionStep(this->GetModelPart());
    mpBuilderAndSolver->FinalizeSolutionStep(mpScheme, this->GetModelPart(), mpA, mpDx, mpb);

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

    return this->SolveIteration();

    KRATOS_CATCH("")
  }


  /**
   * @brief Solves the current iteration. This function returns true if a solution has been found, false otherwise.
   */
  bool SolveIteration() override
  {
    KRATOS_TRY

    // Warning info
    if(!SparseSpaceType::Size((*mpDx)))
      KRATOS_WARNING("DOFS") << "solution has zero size, no free DOFs" << std::endl;

    // Initialize Iteration
    mpScheme->InitializeNonLinearIteration(this->GetModelPart());

    //function to perform the building and the solving phase.
    if(this->mOptions.IsNot(LocalFlagType::CONSTANT_SYSTEM_MATRIX)){

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

    // EchoInfo
    mpBuilderAndSolver->EchoInfo(this->GetModelPart(), (*mpA), (*mpDx), (*mpb));

    // Updating the results
    this->Update();

    // Finalize Iteration
    mpScheme->FinalizeNonLinearIteration(this->GetModelPart());

    return true;

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

    return 0;

    KRATOS_CATCH("")
  }


  ///@}
  ///@name Access
  ///@{

  /**
   * @brief This sets the level of echo for the solving strategy
   * @param Level of echo for the solving strategy
   * @details
   * {
   * 0 -> Mute... no echo at all
   * 1 -> Printing time and basic informations
   * 2 -> Printing linear solver data
   * 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
   * }
   */
  void SetEchoLevel(const int Level) override
  {
    BaseType::SetEchoLevel(Level);
    mpBuilderAndSolver->SetEchoLevel(Level);
  }


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

  SystemVectorPointerType mpDx; /// The incremement in the solution
  SystemVectorPointerType mpb; /// The RHS vector of the system of equations
  SystemMatrixPointerType mpA; /// The LHS matrix of the system of equations

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
    if( mpScheme->IsNot(LocalFlagType::INITIALIZED) )
      mpScheme->Initialize(this->GetModelPart());

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

    mpScheme->Predict(this->GetModelPart(), mpBuilderAndSolver->GetDofSet(), (*mpDx));

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

    mpScheme->Update(this->GetModelPart(), mpBuilderAndSolver->GetDofSet(), (*mpDx));

    KRATOS_CATCH("")
  }


  /**
   * @brief Performs all the required operations to reform dofs
   */
  void SetSystemDofs()
  {
    KRATOS_TRY

   if (this->mEchoLevel >= 2)
      KRATOS_INFO(" Reform Dofs ") << " Flag = " <<this->mOptions.Is(LocalFlagType::REFORM_DOFS) << std::endl;

    //set up the system, operation performed just once unless it is required to reform the dof set at each iteration

    //setting up the list of the DOFs to be solved
    double begin_time = OpenMPUtils::GetCurrentTime();
    mpBuilderAndSolver->SetUpDofSet(mpScheme, this->GetModelPart());
    double end_time = OpenMPUtils::GetCurrentTime();
    if (this->mEchoLevel >= 2)
      KRATOS_INFO("setup_dofs_time") << "setup_dofs_time : " << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;

    //shaping correctly the system
    begin_time = OpenMPUtils::GetCurrentTime();
    mpBuilderAndSolver->SetUpSystem();
    end_time = OpenMPUtils::GetCurrentTime();
    if (this->mEchoLevel >= 2)
      KRATOS_INFO("setup_system_time") << ": setup_system_time : " << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;    //set up the system, operation performed just once unless it is required to reform the dof set at each iteration

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
  LinearStrategy(const LinearStrategy &Other){};

  ///@}

}; /// Class LinearStrategy

///@}

///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

///@} addtogroup block

}  // namespace Kratos.
#endif // KRATOS_LINEAR_STRATEGY_H_INCLUDED defined
