//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               March 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SOLUTION_STRATEGY_H_INCLUDED)
#define KRATOS_SOLUTION_STRATEGY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/solution_scheme.hpp"
#include "custom_solvers/solution_builders_and_solvers/solution_builder_and_solver.hpp"

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

/** @brief Solution strategy base class
 *  @details This is the base class for the strategies
 */

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class SolutionStrategy : public Flags
{
 public:

  ///@name Type Definitions
  ///@{
  typedef SolverLocalFlags                                                           LocalFlagType;
  typedef ModelPart::DofsArrayType                                                   DofsArrayType;

  typedef typename TSparseSpace::MatrixType                                       SystemMatrixType;
  typedef typename TSparseSpace::VectorType                                       SystemVectorType;
  typedef typename TSparseSpace::MatrixPointerType                         SystemMatrixPointerType;
  typedef typename TSparseSpace::VectorPointerType                         SystemVectorPointerType;

  typedef SolutionScheme<TSparseSpace, TDenseSpace>                                     SchemeType;
  typedef SolutionBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>  BuilderAndSolverType;


  /// Pointer definition of SolutionStrategy
  KRATOS_CLASS_POINTER_DEFINITION(SolutionStrategy);


  ///@}
  ///@name Life Cycle
  ///@{


  /// Constructor.
  SolutionStrategy(ModelPart& rModelPart) : Flags(), mrModelPart(rModelPart) { mEchoLevel = 0; }

  /// Constructor.
  SolutionStrategy(ModelPart& rModelPart, Flags& rOptions) : Flags(), mOptions(rOptions), mrModelPart(rModelPart) {mEchoLevel = 0; }

  /// Destructor.
  ~SolutionStrategy() override {}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  /**
   * @brief The problem of interest is solved.
   * @details
   * {
   * This function calls sequentially: InitializeSolutionStep(), SolveSolutionStep() and FinalizeSolutionStep().
   * All those functions can otherwise be called separately.
   * }
   */
  virtual bool Solve()
  {
    KRATOS_TRY

    this->InitializeSolutionStep();

    bool converged = this->SolveSolutionStep();

    // implementation of the adaptive time reduction
    if( this->IsNot(LocalFlagType::ADAPTIVE_SOLUTION) )
      converged = true;
    
    if(converged)
      this->FinalizeSolutionStep();

    return converged;

    KRATOS_CATCH("")
  }

  /**
   * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
   * @details this function must be called only once per step.
   */
  virtual void InitializeSolutionStep() {}

  /**
   * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
   * @details this function must be called only once per step.
   */
  virtual void FinalizeSolutionStep() {}

  /**
   * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
   */
  virtual bool SolveSolutionStep() {return true;}

  /**
   * @brief Solves the current iteration. This function returns true if a solution has been found, false otherwise.
   */
  virtual bool SolveIteration() {return true;}

  /**
   * @brief Clears the internal storage
   */
  virtual void Clear() {}

  /**
   * @brief Function to perform expensive checks.
   * @details It is designed to be called ONCE to verify that the input is correct.
   */
  virtual int Check()
  {
    KRATOS_TRY

    for (ModelPart::ElementsContainerType::iterator it_elem = GetModelPart().ElementsBegin();
         it_elem != GetModelPart().ElementsEnd(); it_elem++)
    {
      it_elem->Check(GetModelPart().GetProcessInfo());
    }

    for (ModelPart::ConditionsContainerType::iterator it_cond = GetModelPart().ConditionsBegin();
         it_cond != GetModelPart().ConditionsEnd(); it_cond++)
    {
      it_cond->Check(GetModelPart().GetProcessInfo());
    }

    return 0;

    KRATOS_CATCH("")
  }

  ///@}
  ///@name Access
  ///@{

  /**
   * @brief This sets the level of echo for the solution strategy
   * @param Level of echo for the solution strategy
   * @details
   * {
   * 0 -> Mute... no echo at all
   * 1 -> Printing time and basic informations
   * 2 -> Printing linear solver data
   * 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
   * }
   */
  virtual void SetEchoLevel(const int Level)
  {
    mEchoLevel = Level;
  }

  /**
   * @brief This returns the level of echo for the solution strategy
   * @details
   * {
   * 0 -> Mute... no echo at all
   * 1 -> Printing time and basic informations
   * 2 -> Printing linear solver data
   * 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
   * }
   * @return Level of echo for the solution strategy
   */
  virtual int GetEchoLevel()
  {
    return mEchoLevel;
  }


  /**
   * @brief Sets strategy options
   */
  void SetOptions(Flags& rOptions)
  {
    mOptions = rOptions;
  }

  /**
   * @brief Get strategy options
   @return mOptions: options member variable
   */
  Flags& GetOptions()
  {
    return mOptions;
  }


  /**
   * @brief This method gets the flag mMaxIterationNumber
   * @return mMaxIterationNumber: This is the maximum number of on linear iterations
   */
  virtual unsigned int GetMaxIterationNumber()
  {
    return 0;
  }

  /**
   * @brief Operations to get the pointer to the model
   * @return mrModelPart: The model part member variable
   */
  inline ModelPart& GetModelPart()
  {
    return mrModelPart;
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

  // Flags to set options
  Flags mOptions;

  // Level of echo for the solution strategy
  int mEchoLevel;

  ///@}
  ///@name Protected Operators
  ///@{

  /**
   * @brief Initialization of member variables and prior operations
   */
  virtual void Initialize(){};

  /**
   * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
   * values of the solution step of interest are assumed equal to the old values
   */
  virtual void Predict(){};

  /**
   * @brief Operation to update the solution ... if it is not called a trivial updater is used in which the
   * values of the solution step of interest are assumed equal to the old values
   */
  virtual void Update(){};

  /**
   * @brief Finalization of member variables and prior operations
   */
  virtual void Finalize(){};

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

 private:

  ///@}
  ///@name Static Member Variables
  ///@{
  ///@}
  ///@name Member Variables
  ///@{

  ModelPart& mrModelPart;

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
  SolutionStrategy(SolutionStrategy const& Other) {};

  ///@}

}; /// Class SolutionStrategy

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.
#endif // KRATOS_SOLUTION_STRATEGY_H_INCLUDED defined
