//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               March 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SEGREGATED_STRATEGY_H_INCLUDED)
#define KRATOS_SEGREGATED_STRATEGY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_strategies/solution_strategy.hpp"

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
 * @class SegregatedStrategy
 * @brief This is the base class for a segregated strategy for the same model part
 */

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver // = LinearSolver<TSparseSpace,TDenseSpace>
          >
class SegregatedStrategy : public SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
 public:
  ///@name Type Definitions
  ///@{

  // Counted pointer of ClassName
  KRATOS_CLASS_POINTER_DEFINITION(SegregatedStrategy);

  typedef SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>            BaseType;
  typedef typename BaseType::Pointer                                     BasePointerType;
  typedef typename BaseType::LocalFlagType                                 LocalFlagType;

  typedef typename std::vector<BasePointerType>                  StrategiesContainerType;
  typedef typename StrategiesContainerType::iterator     StrategiesContainerIteratorType;

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
   * @param rOptions The solution options
   */
  SegregatedStrategy(ModelPart& rModelPart,
                     Flags& rOptions)
      : SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, rOptions)
  {
    KRATOS_TRY

    KRATOS_CATCH("")
  }

  /**
   * Constructor
   * @param rModelPart The model part of the problem
   * @param rOptions The solution options
   * @param rStrategies Vector with the solution strategies
   */
  SegregatedStrategy(ModelPart& rModelPart,
                     Flags& rOptions,
                     StrategiesContainerType& rStrategies)
      : SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, rOptions)
  {
    KRATOS_TRY

    mStrategies = rStrategies;

    KRATOS_CATCH("")
  }

  /**
   * @brief Destructor.
   */
  ~SegregatedStrategy() override
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
    if(this->mOptions.Is(LocalFlagType::IMPLEX) && this->GetModelPart().GetProcessInfo().Has(IMPLEX))
      this->GetModelPart().GetProcessInfo()[IMPLEX] = 1;

    //prints informations about the current time
    //KRATOS_INFO("") << "  [STEP:" << this->GetModelPart().GetProcessInfo()[STEP] << "  TIME: "<< this->GetModelPart().GetProcessInfo()[TIME]<< "]\n" << LoggerMessage::Category::STATUS;

    int counter = 0;
    for(StrategiesContainerIteratorType it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      this->GetModelPart().GetProcessInfo().SetValue(SEGREGATED_STEP, counter);
      (*it)->InitializeSolutionStep();
      ++counter;
    }

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

    int counter = 0;
    for(StrategiesContainerIteratorType it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      this->GetModelPart().GetProcessInfo().SetValue(SEGREGATED_STEP, counter);
      (*it)->FinalizeSolutionStep();
      ++counter;
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

    //setting the iteration number
    this->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

    //maximum number of iterations
    unsigned int MaxIterationNumber = 1;
    for(StrategiesContainerIteratorType it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      if(MaxIterationNumber < (*it)->GetMaxIterationNumber())
        MaxIterationNumber = (*it)->GetMaxIterationNumber();
    }

    //setting the iteration number
    this->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

    this->Set(LocalFlagType::CONVERGED, this->SolveIteration());

    //iteration cycle... performed only for NonLinearProblems
    while( this->IsNot(LocalFlagType::CONVERGED) && ++iteration_number < MaxIterationNumber)
    {
      //setting the iteration number
      this->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

      this->Set(LocalFlagType::CONVERGED, this->SolveIteration());
    }

    //plots a warning if the maximum number of iterations is exceeded
    if(iteration_number >= MaxIterationNumber)
    {
      if( this->GetEchoLevel() >= 0 )
        KRATOS_INFO("  [The iterative loop interrupted] ") << "[" << iteration_number << " iterations performed] \n";
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

    //convergence vector
    std::vector<bool> convergences(mStrategies.size());
    std::fill(convergences.begin(), convergences.end(), false);

    int counter = 0;
    for(StrategiesContainerIteratorType it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      this->GetModelPart().GetProcessInfo().SetValue(SEGREGATED_STEP, counter);
      convergences[counter] = (*it)->SolveIteration();
      ++counter;
    }

    bool convergence = std::all_of(convergences.begin(), convergences.end(), [](bool const n){return n == true;});
    if(convergence == true){
      if( this->GetEchoLevel() >= 0 )
        KRATOS_INFO("  [Convergence Achieved] ") << "[" << this->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] << " iterations performed]\n";
    }

    return (convergence);

    KRATOS_CATCH("")
  }


  /**
   * @brief Clears the internal storage
   */
  void Clear() override
  {
    KRATOS_TRY

    int counter = 0;
    for(StrategiesContainerIteratorType it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      this->GetModelPart().GetProcessInfo().SetValue(SEGREGATED_STEP, counter);
      (*it)->Clear();
      ++counter;
    }

    KRATOS_CATCH("")
  }

  /**
   * @brief Function to perform expensive checks.
   * @details It is designed to be called ONCE to verify that the input is correct.
   */
  int Check() override
  {
    KRATOS_TRY

    int counter = 0;
    for(StrategiesContainerIteratorType it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      this->GetModelPart().GetProcessInfo().SetValue(SEGREGATED_STEP, counter);
      (*it)->Check();
      ++counter;
    }

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
    for(StrategiesContainerIteratorType it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      (*it)->SetEchoLevel(Level);
    }

  }

  /**
   * @brief Set method for the time scheme
   * @param pScheme The pointer to the time scheme considered
   */
  void AddStrategy(const BasePointerType pStrategy)
  {
    mStrategies.push_back(pStrategy);
  };

  /**
   * @brief Get method for the time scheme
   * @return mpScheme: The pointer to the time scheme considered
   */
  BasePointerType GetStrategy(unsigned int Label)
  {
    return mStrategies[Label];
  };


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

  StrategiesContainerType  mStrategies;

  ///@}
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
  ///@{
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
  SegregatedStrategy(const SegregatedStrategy &Other){};

  ///@}

}; /// Class SegregatedStrategy

///@}

///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

///@} addtogroup block

}  // namespace Kratos.
#endif // KRATOS_SEGREGATED_STRATEGY_H_INCLUDED defined
