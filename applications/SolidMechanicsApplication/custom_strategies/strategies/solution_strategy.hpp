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

#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

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

/** @brief Solving strategy local flags class definition
 *  @details This is the base class for strategy local flags
 */
class SolverLocalFlags
{
 public:
  /// Flags for the Strategy control

  /// external type:
  KRATOS_DEFINE_LOCAL_FLAG( INITIALIZED );
  KRATOS_DEFINE_LOCAL_FLAG( CONVERGED );

  /// internal type:
  KRATOS_DEFINE_LOCAL_FLAG( INITIALIZED_ELEMENTS );
  KRATOS_DEFINE_LOCAL_FLAG( INITIALIZED_CONDITIONS );
  
  /// Flags for the Strategy options
  KRATOS_DEFINE_LOCAL_FLAG( MOVE_MESH );
  KRATOS_DEFINE_LOCAL_FLAG( REFORM_DOFS );
  KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_SYSTEM_LHS );
  KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_REACTIONS );
  KRATOS_DEFINE_LOCAL_FLAG( IMPLEX );
};

/** @brief Solving strategy base class
 *  @details This is the base class for strategies
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
  typedef SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>  SolutionStrategyType;

  typedef SolverLocalFlags                                                   LocalFlagType;

  typedef typename TSparseSpace::MatrixType                               SystemMatrixType;
  typedef typename TSparseSpace::VectorType                               SystemVectorType;
  typedef typename TSparseSpace::MatrixPointerType                 SystemMatrixPointerType;
  typedef typename TSparseSpace::VectorPointerType                 SystemVectorPointerType;

  // ---
  
  typedef Scheme<TSparseSpace, TDenseSpace>                                     SchemeType;
  typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>  BuilderAndSolverType;

  // ---
  
  typedef typename ModelPart::DofsArrayType                                  DofsArrayType;

  /// Pointer definition of Constitutive3DLaw 
  KRATOS_CLASS_POINTER_DEFINITION(SolutionStrategy);

  
  ///@}
  ///@name Life Cycle
  ///@{


  /// Default constructor.
  SolutionStrategy(ModelPart& rModelPart) : Flags(), mrModelPart(rModelPart) {}

  /// Default constructor.
  SolutionStrategy(ModelPart& rModelPart, Flags& rOptions) : Flags(), mOptions(rOptions), mrModelPart(rModelPart) {}
  
  /// Destructor.
  virtual ~SolutionStrategy() {}
    
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
    this->SolveSolutionStep();
    this->FinalizeSolutionStep();
    
    return true;

    KRATOS_CATCH("")
  }
  
  /**
   * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
   * @details A member variable should be used as a flag to make sure this function is called only once per step.
   */
  virtual void InitializeSolutionStep() {}
  
  /**
   * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
   * @details A member variable should be used as a flag to make sure this function is called only once per step.
   */
  virtual void FinalizeSolutionStep() {}

  /**
   * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
   */
  virtual bool SolveSolutionStep() {return true;}

  /**
   * @brief This function is designed to move the mesh
   * @note Be careful it just consider displacements, derive this method to adapt to your own strategies (ALE, FSI, etc...)
   */
  virtual void MoveMesh()
  {
    KRATOS_TRY

    if (GetModelPart().NodesBegin()->SolutionStepsDataHas(DISPLACEMENT_X) == false)
    {
      KRATOS_ERROR << "It is impossible to move the mesh since the DISPLACEMENT variable is not in the Model Part. Either set flag MOVE_MESH to false or add DISPLACEMENT to the list of variables" << std::endl;
    }

    const int nnodes = mrModelPart.NumberOfNodes();
    ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

#pragma omp parallel for
    for(int i = 0; i<nnodes; i++)
    {
      ModelPart::NodesContainerType::iterator it_node = it_begin + i;

      noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
      noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT);
    }

    KRATOS_CATCH("")
  }
  
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

    // Check DISPLACEMENT variables for MOVE_MESH
    if(mOptions.Is(LocalFlagType::MOVE_MESH))
    {
      for (ModelPart::NodesContainerType::iterator itNode = mrModelPart.NodesBegin();
           itNode != mrModelPart.NodesEnd(); itNode++)
      {
        if (itNode->SolutionStepsDataHas(DISPLACEMENT) == false)
        {
          KRATOS_ERROR << "ERROR:: Problem on node with Id " << itNode->Id() << "\nIt is impossible to move the mesh since the DISPLACEMENT variable is not in the rModelPart. Either set flag MOVE_MESH to false or add DISPLACEMENT to the list of variables" << std::endl;
        }
      }
    }
  
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
  virtual void SetEchoLevel(const int Level)
  {
    mEchoLevel = Level;
  }

  /**
   * @brief This returns the level of echo for the solving strategy
   * @details
   * {
   * 0 -> Mute... no echo at all
   * 1 -> Printing time and basic informations
   * 2 -> Printing linear solver data
   * 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
   * }
   * @return Level of echo for the solving strategy
   */
  virtual int GetEchoLevel()
  {
    return mEchoLevel;
  }

  
  /**
   * @brief Sets strategy options
   */
  inline void SetOptions(Flags& rOptions)
  {
    mOptions = rOptions;
  }

  /**
   * @brief Get strategy options
   @return mOptions: options member variable
   */
  inline Flags& GetOptions()
  {
    return mOptions;
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

  // Flags to set options
  Flags mOptions;
  
  // Level of echo for the solving strategy
  int mEchoLevel;
  
  ///@}
  ///@name Protected member Variables
  ///@{

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
  
/**
 * Flags for the Strategy control
 */
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, INITIALIZED,               0 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, CONVERGED,                 1 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, INITIALIZED_ELEMENTS,      2 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, INITIALIZED_CONDITIONS,    3 );

/**
 * Flags for the Strategy options
 */
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, MOVE_MESH,                 0 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, REFORM_DOFS,               1 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, CONSTANT_SYSTEM_LHS,       2 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, COMPUTE_REACTIONS,         3 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, IMPLEX,                    4 );


///@}
///@name Input and output
///@{

///@}

///@} addtogroup block
  
}  // namespace Kratos.
#endif // KRATOS_SOLUTION_STRATEGY_H_INCLUDED defined
