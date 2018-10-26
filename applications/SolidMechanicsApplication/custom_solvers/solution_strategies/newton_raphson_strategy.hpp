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
#include "custom_solvers/solution_strategies/linear_strategy.hpp"

//default convergence criterion
#include "custom_solvers/convergence_criteria/convergence_criterion.hpp"

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
class NewtonRaphsonStrategy : public LinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
 public:
  ///@name Type Definitions
  ///@{

  // Counted pointer of ClassName
  KRATOS_CLASS_POINTER_DEFINITION(NewtonRaphsonStrategy);

  typedef LinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>              BaseType;

  typedef typename BaseType::LocalFlagType                                 LocalFlagType;

  typedef ConvergenceCriterion<TSparseSpace, TDenseSpace>       ConvergenceCriterionType;

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
      : LinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pBuilderAndSolver, rOptions)
  {
    KRATOS_TRY

    // Saving the convergence criteria to be used
    mpConvergenceCriteria = pConvergenceCriterion;

    // Maximum iterations allowed
    mMaxIterationNumber = MaxIterations;

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
      : NewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, Kratos::make_shared<BlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(pLinearSolver), pConvergenceCriterion, rOptions, MaxIterations)
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

    //set implex
    if(this->mOptions.Is(LocalFlagType::IMPLEX))
      this->GetModelPart().GetProcessInfo().SetValue(IMPLEX, true);

    BaseType::InitializeSolutionStep();

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
    if(this->mOptions.Is(LocalFlagType::IMPLEX)){
      this->GetModelPart().GetProcessInfo().SetValue(IMPLEX, false);
      if(this->mOptions.IsNot(LocalFlagType::COMPUTE_REACTIONS))
        this->mpBuilderAndSolver->BuildRHS(this->mpScheme, this->GetModelPart(), (*this->mpb));
    }

    BaseType::FinalizeSolutionStep();

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

    //setting the iteration number
    this->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

    this->Set(LocalFlagType::CONVERGED, this->SolveIteration());

    //iteration cycle... performed only for NonLinearProblems
    while( this->IsNot(LocalFlagType::CONVERGED) && iteration_number++ < mMaxIterationNumber)
    {
      //setting the iteration number
      this->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

      this->Set(LocalFlagType::CONVERGED, this->SolveIteration());
    }

    //plots a warning if the maximum number of iterations is exceeded
    if(iteration_number >= mMaxIterationNumber)
    {
      if( this->GetEchoLevel() >= 0 )
        KRATOS_INFO("  [Iterative loop interrupted] ") << "[" << iteration_number << " iterations performed] \n";
    }

    return (this->Is(LocalFlagType::CONVERGED));

    KRATOS_CATCH("")

  }


  /**
   * @brief Solves the iteration. This function returns true if a solution has been found, false otherwise.
   */
  bool SolveIteration() override
  {
    KRATOS_TRY

    bool is_converged = false;

    // Warning info
    if(!SparseSpaceType::Size((*this->mpDx)))
      KRATOS_WARNING("DOFS") << "solution has zero size, no free DOFs" << std::endl;

    // Initialize Iteration
    this->mpScheme->InitializeNonLinearIteration(this->GetModelPart());

    is_converged = mpConvergenceCriteria->PreCriteria(this->GetModelPart(), this->mpBuilderAndSolver->GetDofSet(), (*this->mpA), (*this->mpDx), (*this->mpb));

    // Function to perform the building and the solving phase.
    if(this->mOptions.IsNot(LocalFlagType::CONSTANT_SYSTEM_MATRIX)){

      TSparseSpace::SetToZero((*this->mpA));
      TSparseSpace::SetToZero((*this->mpDx));
      TSparseSpace::SetToZero((*this->mpb));

      this->mpBuilderAndSolver->BuildAndSolve(this->mpScheme, this->GetModelPart(), (*this->mpA), (*this->mpDx), (*this->mpb));
    }
    else{

      TSparseSpace::SetToZero((*this->mpDx));
      TSparseSpace::SetToZero((*this->mpb));

      this->mpBuilderAndSolver->BuildRHSAndSolve(this->mpScheme, this->GetModelPart(), (*this->mpA), (*this->mpDx), (*this->mpb));
    }

    // EchoInfo
    this->mpBuilderAndSolver->EchoInfo(this->GetModelPart(), (*this->mpA), (*this->mpDx), (*this->mpb));

    // Updating the results
    this->Update();

    // Finalize Iteration
    this->mpScheme->FinalizeNonLinearIteration(this->GetModelPart());

    if(is_converged == true)
    {
      //initialisation of the convergence criteria (after first calculation only)
      if( this->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] == 1 ){
        mpConvergenceCriteria->InitializeSolutionStep(this->GetModelPart(), this->mpBuilderAndSolver->GetDofSet(), (*this->mpA), (*this->mpDx), (*this->mpb));
      }

      if(mpConvergenceCriteria->Is(CriterionLocalFlags::UPDATE_RHS))
      {
        TSparseSpace::SetToZero((*this->mpb));

        this->mpBuilderAndSolver->BuildRHS(this->mpScheme, this->GetModelPart(), (*this->mpb));
      }

      is_converged = mpConvergenceCriteria->PostCriteria(this->GetModelPart(), this->mpBuilderAndSolver->GetDofSet(), (*this->mpA), (*this->mpDx), (*this->mpb));
    }

    return is_converged;

    KRATOS_CATCH("")
  }


  /**
   * @brief Clears the internal storage
   */
  void Clear() override
  {
    KRATOS_TRY

    BaseType::Clear();

    KRATOS_CATCH("")
  }

  /**
   * @brief Function to perform expensive checks.
   * @details It is designed to be called ONCE to verify that the input is correct.
   */
  int Check() override
  {
    KRATOS_TRY

    //check linear strategy
    BaseType::Check();

    //check the convergence criterion
    mpConvergenceCriteria->Check(this->GetModelPart());

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
  }

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
  unsigned int GetMaxIterationNumber() override
  {
    return mMaxIterationNumber;
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

  typename ConvergenceCriterionType::Pointer mpConvergenceCriteria; /// The pointer to the convergence criteria employed

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

    BaseType::Initialize();

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
