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

  typedef typename std::vector<BasePointerType>                  StrategiesContainerType;
  
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
  SegregatedStrategy(ModelPart& rModelPart,
                     Flags& rOptions)
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
   * @param pConvergenceCriteria The convergence criteria employed
   * @param rOptions The solution options
   * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
   */
  SegregatedStrategy(ModelPart& rModelPart,
                        typename SchemeType::Pointer pScheme,
                        typename LinearSolverType::Pointer pLinearSolver,
                        typename ConvergenceCriterionType::Pointer pConvergenceCriterion,
                        Flags& rOptions,
                        unsigned int MaxIterations = 30)
      : SegregatedStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, typename BuilderAndSolverType::Pointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pLinearSolver)), pConvergenceCriterion, rOptions, MaxIterations)
  {}
    

  /** 
   * @brief Destructor.
   * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
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
                    
    //initialize system elements, conditions..
    if(this->IsNot(LocalFlagType::INITIALIZED))
      this->Initialize();

    //set implex
    if(this->mOptions.Is(LocalFlagType::IMPLEX) && this->GetModelPart().GetProcessInfo().Has(IMPLEX))
      this->GetModelPart().GetProcessInfo()[IMPLEX] = 1;
    
    //prints informations about the current time
    //KRATOS_INFO("") << "  [STEP:" << this->GetModelPart().GetProcessInfo()[STEP] << "  TIME: "<< this->GetModelPart().GetProcessInfo()[TIME]<< "]\n" << LoggerMessage::Category::STATUS;

    unsigned int counter = 0;
    for(StrategiesContainerType::iterator it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      this->GetModelPart().GetProcessInfo().SetValue(FRACTIONAL_STEP, ++counter);
      it->InitializeSolutionStep();
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
    
    unsigned int counter = 0;
    for(StrategiesContainerType::iterator it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      this->GetModelPart().GetProcessInfo().SetValue(FRACTIONAL_STEP, ++counter);
      it->FinalizeSolutionStep();
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

    //Iteration Cicle... performed only for NonLinearProblems
    while( this->IsNot(LocalFlagType::CONVERGED) && iteration_number++ < mMaxIterationNumber)
    {
      //setting the iteration number
      this->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

      bool is_converged = false;
      unsigned int counter = 0;
      for(StrategiesContainerType::iterator it= mStrategies.begin(); it!= mStrategies.end(); ++it)
      {
        this->GetModelPart().GetProcessInfo().SetValue(FRACTIONAL_STEP, ++counter);
        is_converged = it->SolveIteration();
      }

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

    unsigned int counter = 0;
    for(StrategiesContainerType::iterator it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      this->GetModelPart().GetProcessInfo().SetValue(FRACTIONAL_STEP, ++counter);
      it->Clear();
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

    unsigned int counter = 0;
    for(StrategiesContainerType::iterator it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      this->GetModelPart().GetProcessInfo().SetValue(FRACTIONAL_STEP, ++counter);
      it->Check();
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
    for(StrategiesContainerType::iterator it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      it->SetEchoLevel(Level);
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

  /**
   * @brief Initialization of member variables and prior operations
   */    
  void Initialize() override
  {
    KRATOS_TRY

    unsigned int counter = 0;
    for(StrategiesContainerType::iterator it= mStrategies.begin(); it!= mStrategies.end(); ++it)
    {
      this->GetModelPart().GetProcessInfo().SetValue(FRACTIONAL_STEP, ++counter);
      it->Initialize();
    }

    this->Set(LocalFlagType::INITIALIZED,true);
    
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
