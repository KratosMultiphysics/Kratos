//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:           MSantasusana $
//   Last modified by:    $Co-Author:         JMCarbonell $
//   Date:                $Date:               April 2014 $
//   Revision:            $Revision:           March 2018 $
//
//

#if !defined(KRATOS_EXPLICIT_STRATEGY_H_INCLUDED)
#define KRATOS_EXPLICIT_STRATEGY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_strategies/solution_strategy.hpp"
#include "solid_mechanics_application_variables.h"

//default builder and solver
#include "custom_solvers/solution_builders_and_solvers/explicit_builder_and_solver.hpp"


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

template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ExplicitSolutionStrategy : public SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
 public:

  ///@name Type Definitions
  ///@{

  // Counted pointer of ClassName

  KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolutionStrategy);

  typedef SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>            BaseType;

  typedef typename BaseType::LocalFlagType                                 LocalFlagType;

  typedef ConvergenceCriterion<TSparseSpace, TDenseSpace>       ConvergenceCriterionType;

  typedef typename BaseType::BuilderAndSolverType                   BuilderAndSolverType;

  typedef typename BaseType::SchemeType                                       SchemeType;

  typedef TSparseSpace                                                   SparseSpaceType;

  typedef typename BaseType::DofsArrayType                                 DofsArrayType;

  typedef typename BaseType::SystemMatrixType                           SystemMatrixType;

  typedef typename BaseType::SystemVectorType                           SystemVectorType;

  typedef typename BaseType::SystemMatrixPointerType             SystemMatrixPointerType;

  typedef typename BaseType::SystemVectorPointerType             SystemVectorPointerType;

  ///@}
  ///@name Life Cycle

  ///@{


  /// Constructor.
  ExplicitSolutionStrategy(ModelPart& rModelPart,
                   typename SchemeType::Pointer pScheme,
                   Flags& rOptions)
      : SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, rOptions)
  {
    KRATOS_TRY


    //saving the scheme
    mpScheme = pScheme;

    //create explicit builder
    mpBuilderAndSolver = Kratos::make_shared<ExplicitBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >();

    //set lumped mass matrix by default
    if( this->GetModelPart().GetProcessInfo().Has(COMPUTE_LUMPED_MASS_MATRIX) )
      this->GetModelPart().GetProcessInfo()[COMPUTE_LUMPED_MASS_MATRIX] = true;
    KRATOS_CATCH( "" )
  }

  /// Destructor.
  ~ExplicitSolutionStrategy() override {}

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

    //initial operations ... things that are constant over the Solution Step
    mpBuilderAndSolver->InitializeSolutionStep(mpScheme, this->GetModelPart(), mpA, mpDx, mpb);

    //initial operations ... things that are constant over the Solution Step
    mpScheme->InitializeSolutionStep(this->GetModelPart());


    KRATOS_CATCH( "" )
  }


  /**
   * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
   */
  void FinalizeSolutionStep() override
  {
    KRATOS_TRY

    //Finalization of the solution step

    //calculate reactions if required
    // if(mOptions.Is(LocalFlagType::COMPUTE_REACTIONS))
    //   mpBuilderAndSolver->CalculateReactions(pScheme, this->GetModelPart(), (*mpA), (*mpDx), (*mpb));

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

    //compute nodal mass and inertia
    if(this->mOptions.IsNot(LocalFlagType::CONSTANT_SYSTEM_MATRIX))
      mpBuilderAndSolver->BuildLHS(mpScheme, this->GetModelPart(), (*mpA));


    mpBuilderAndSolver->BuildRHS(mpScheme, this->GetModelPart(), (*mpb));

    //update explicitly integrates the equation of motion
    this->Update();

    return true;

    KRATOS_CATCH( "" )
  }

  //**********************************************************************
  //**********************************************************************

  void Clear() override
  {
    KRATOS_TRY

    mpBuilderAndSolver->Clear();
    mpScheme->Clear();

    KRATOS_CATCH( "" )
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
  SystemMatrixPointerType mpA; /// The LHS matrix of the system of equations (dummy)


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
    if (mpScheme->IsNot(LocalFlagType::INITIALIZED))
      mpScheme->Initialize(this->GetModelPart());

    // //Initialize The Elements - OPERATIONS TO BE DONE ONCE
    // if (mpScheme->ElementsAreInitialized() == false)
    //   mpScheme->InitializeElements(this->GetModelPart());

    // //Initialize The Conditions- OPERATIONS TO BE DONE ONCE
    // if (mpScheme->ConditionsAreInitialized() == false)
    //   mpScheme->InitializeConditions(this->GetModelPart());

    //compute nodal mass and inertia
    if(this->mOptions.Is(LocalFlagType::CONSTANT_SYSTEM_MATRIX))
      mpBuilderAndSolver->BuildLHS(mpScheme, this->GetModelPart(), (*mpA));

    this->Set(LocalFlagType::INITIALIZED,true);


    KRATOS_CATCH( "" )
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
   * function to perform expensive checks.
   * It is designed to be called ONCE to verify that the input is correct.
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

    KRATOS_CATCH( "" )
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
  ExplicitSolutionStrategy(const ExplicitSolutionStrategy& Other){};

  ///@}

}; /// Class ExplicitSolutionStrategy

///@}

///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

///@} addtogroup block

}  // namespace Kratos.
#endif // KRATOS_EXPLICIT_STRATEGY_H_INCLUDED  defined
