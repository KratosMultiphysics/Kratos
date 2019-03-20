//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               March 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SOLUTION_BUILDER_AND_SOLVER_H_INCLUDED)
#define KRATOS_SOLUTION_BUILDER_AND_SOLVER_H_INCLUDED

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


/** Short class definition.
Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains this information.
 */


/** @brief Solution Buider and Solver base class
 *  @details This is the base class for the building and solving the solution system
 */
template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class SolutionBuilderAndSolver : public Flags
{
public:

  ///@name Type Definitions
  ///@{

  /// Pointer definition of SolutionBuilderAndSolver
  KRATOS_CLASS_POINTER_DEFINITION(SolutionBuilderAndSolver);

  typedef SolverLocalFlags                                                   LocalFlagType;

  typedef ModelPart::DofsArrayType                                            DofsArrayType;
  typedef typename TSparseSpace::MatrixType                                SystemMatrixType;
  typedef typename TSparseSpace::VectorType                                SystemVectorType;

  typedef typename TSparseSpace::MatrixPointerType                  SystemMatrixPointerType;
  typedef typename TSparseSpace::VectorPointerType                  SystemVectorPointerType;

  typedef typename TDenseSpace::MatrixType                            LocalSystemMatrixType;
  typedef typename TDenseSpace::VectorType                            LocalSystemVectorType;

  typedef SolutionScheme<TSparseSpace, TDenseSpace>                              SchemeType;
  typedef typename SchemeType::Pointer                                    SchemePointerType;

  typedef typename TLinearSolver::Pointer                           LinearSolverPointerType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default Constructor.
  SolutionBuilderAndSolver() : Flags()
  {
    this->Set(LocalFlagType::DOFS_INITIALIZED, false);

    mpLinearSystemSolver = nullptr;
    mEchoLevel = 0;
  }

  /// Constructor.
  SolutionBuilderAndSolver(LinearSolverPointerType pLinearSystemSolver) : Flags()
  {
    this->Set(LocalFlagType::DOFS_INITIALIZED, false);

    mpLinearSystemSolver = pLinearSystemSolver;
    mEchoLevel = 0;
  }

  /// Destructor.
  ~SolutionBuilderAndSolver() override {}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{


  /**
   * @brief Function to perform the building of the LHS,
   * @details Depending on the implementation choosen the size of the matrix could be
   * @details equal to the total number of Dofs or to the number of non constrained dofs
   */
  virtual void BuildLHS(SchemePointerType pScheme,
                        ModelPart& rModelPart,
                        SystemMatrixType& rA)
  {
  }

  /**
   * @brief Function to perform the build of the RHS.
   * @details The vector could be sized as the total number of dofs or as the number of non constrained ones
   */
  virtual void BuildRHS(SchemePointerType pScheme,
                        ModelPart& rModelPart,
                        SystemVectorType& rb)
  {
  }

  /**
   * @brief Function to perform the building of the LHS and RHS
   * @details Equivalent (but generally faster) then performing BuildLHS and BuildRHS
   */
  virtual void Build(SchemePointerType pScheme,
                     ModelPart& rModelPart,
                     SystemMatrixType& rA,
                     SystemVectorType& rb)
  {
  }

  /**
   * @brief This is a call to the linear system solver
   */
  virtual void SystemSolve(SystemMatrixType& rA,
                           SystemVectorType& rDx,
                           SystemVectorType& rb
                           )
  {
  }

  /**
   * @brief Function to perform the building and solving phase at the same time.
   * @details It is ideally the fastest and safer function to use when it is possible to solve just after building
   */
  virtual void BuildAndSolve(SchemePointerType pScheme,
                             ModelPart& rModelPart,
                             SystemMatrixType& rA,
                             SystemVectorType& rDx,
                             SystemVectorType& rb)
  {
  }

  /**
   * @brief Function to perform the building of the RHS and solving phase at the same time.
   * @details  It corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
   */
  virtual void BuildRHSAndSolve(SchemePointerType pScheme,
                                ModelPart& rModelPart,
                                SystemMatrixType& rA,
                                SystemVectorType& rDx,
                                SystemVectorType& rb)
  {
  }

  /**
   * @brief applies the dirichlet conditions.
   * @details This operation may be very heavy or completely
   * @details unexpensive depending on the implementation choosen and on how the System Matrix
   * @details is built. For explanation of how it works for a particular implementation the user
   * @details should refer to the particular Builder And Solver choosen
   */
  virtual void ApplyDirichletConditions(SchemePointerType pScheme,
                                        ModelPart& rModelPart,
                                        SystemMatrixType& rA,
                                        SystemVectorType& rDx,
                                        SystemVectorType& rb)
  {
  }


  // /**
  //  * @brief Performs all the required operations to reform dofs
  //  */
  // virtual void SetSystemDofs(SchemePointerType pScheme,
  //                            ModelPart& rModelPart)
  // {
  //   KRATOS_TRY

  //   if (this->mEchoLevel >= 2)
  //     KRATOS_INFO(" Reform Dofs ") << " Flag = " <<this->mOptions.Is(LocalFlagType::REFORM_DOFS) << std::endl;

  //   //set up the system, operation performed just once unless it is required to reform the dof set at each iteration

  //   //setting up the list of the DOFs to be solved
  //   double begin_time = OpenMPUtils::GetCurrentTime();
  //   this->SetUpDofSet(rModelPart);
  //   double end_time = OpenMPUtils::GetCurrentTime();
  //   if (this->mEchoLevel >= 2)
  //     KRATOS_INFO("setup_dofs_time") << "setup_dofs_time : " << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;

  //   //shaping correctly the system
  //   begin_time = OpenMPUtils::GetCurrentTime();
  //   this->SetUpSystem();
  //   end_time = OpenMPUtils::GetCurrentTime();
  //   if (this->mEchoLevel >= 2)
  //     KRATOS_INFO("setup_system_time") << ": setup_system_time : " << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;

  //   KRATOS_CATCH("")
  // }


  /**
   * @brief Builds the list of the DofSets involved in the problem by "asking" to each element and condition its Dofs.
   * @details The list of dofs is stores insde the SolutionBuilderAndSolver as it is closely connected to the way the matrix and RHS are built
  */
  virtual void SetUpDofSet(SchemePointerType pScheme,
                           ModelPart& rModelPart)
  {
  }

  /**
   * @brief organises the dofset in order to speed up the building phase
  */
  virtual void SetUpSystem()
  {
  }


  /**
   * @brief Resizes and Initializes the system vectors and matrices after SetUpDofSet and SetUpSytem has been called
  */
  virtual void SetUpSystemMatrices(SchemePointerType pScheme,
                                   ModelPart& rModelPart,
                                   SystemMatrixPointerType& pA,
                                   SystemVectorPointerType& pDx,
                                   SystemVectorPointerType& pb)
  {
  }


  /**
   * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
   * @details this function must be called only once per step.
   */
  virtual void InitializeSolutionStep(SchemePointerType pScheme,
                                      ModelPart& rModelPart,
                                      SystemMatrixPointerType& pA,
                                      SystemVectorPointerType& pDx,
                                      SystemVectorPointerType& pb)
  {
  }

  /**
   * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
   * @details this function must be called only once per step.
   */
  virtual void FinalizeSolutionStep(SchemePointerType pScheme,
                                    ModelPart& rModelPart,
                                    SystemMatrixPointerType& pA,
                                    SystemVectorPointerType& pDx,
                                    SystemVectorPointerType& pb)
  {
  }

  /**
   * @brief Calculates system reactions
   * @details A flag controls if reactions must be calculated
   * @details An internal variable to store the reactions vector is needed
   */
  virtual void CalculateReactions(SchemePointerType pScheme,
                                  ModelPart& rModelPart,
                                  SystemMatrixType& rA,
                                  SystemVectorType& rDx,
                                  SystemVectorType& rb)
  {
  }

  /**
   * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
  */
  virtual void Clear()
  {
    this->mDofSet.clear();

    if(this->mpReactionsVector != nullptr)
      TSparseSpace::Clear(this->mpReactionsVector);

    if(this->mpLinearSystemSolver != nullptr)
      this->mpLinearSystemSolver->Clear();

    if (this->mEchoLevel > 0)
    {
      KRATOS_INFO("Builder and Solver Cleared")  << "Clear Function called" << std::endl;
    }
  }


  /**
   * This function is designed to be called once to perform all the checks needed
   * on the input provided. Checks can be "expensive" as the function is designed
   * to catch user's errors.
   * @param rModelPart
   * @return 0 all ok
   */
  virtual int Check(ModelPart& rModelPart)
  {
    KRATOS_TRY

    return 0;

    KRATOS_CATCH("")
  }


  ///@}
  ///@name Access
  ///@{

  /**
   * @brief Get size of the system
   */
  unsigned int GetEquationSystemSize()
  {
    return mEquationSystemSize;
  }


  /**
   * @brief Get linear solver
   */
  LinearSolverPointerType GetLinearSystemSolver()
  {
    return mpLinearSystemSolver;
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
   * @brief This sets the level of echo for the builder and solver
   * @param Level of echo for the builder and solver
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

  virtual int GetEchoLevel()
  {
    return mEchoLevel;
  }


  /**
   * @brief This method returns the components of the system of equations depending of the echo level
   */
  virtual void EchoInfo(ModelPart& rModelPart,
                        SystemMatrixType& rA,
                        SystemVectorType& rDx,
                        SystemVectorType& rb)
  {
    KRATOS_TRY

    if (this->mEchoLevel == 2) //if it is needed to print the debug info
    {
      KRATOS_INFO("Dx")  << "Solution = " << rDx << std::endl;
      KRATOS_INFO("RHS") << "Vector = " << rb << std::endl;
    }
    else if (this->mEchoLevel == 3) //if it is needed to print the debug info
    {
      KRATOS_INFO("LHS") << "Matrix = " << rA << std::endl;
      KRATOS_INFO("Dx")  << "Solution = " << rDx << std::endl;
      KRATOS_INFO("RHS") << "Vector = " << rb << std::endl;

    }
    else if (this->mEchoLevel == 4) //print to matrix market file
    {
      unsigned int iteration_number = 0;
      if( rModelPart.GetProcessInfo().Has(NL_ITERATION_NUMBER) )
        iteration_number = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];

      std::stringstream matrix_market_name;
      matrix_market_name << "A_" << rModelPart.GetProcessInfo()[TIME] << "_" << iteration_number << ".mm";
      TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rA, false);

      std::stringstream matrix_market_vectname;
      matrix_market_vectname << "b_" << rModelPart.GetProcessInfo()[TIME] << "_" << iteration_number << ".mm.rhs";
      TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), rb);
    }

    KRATOS_CATCH("")

  }


  /**
   * @brief Allows to get the list of system Dofs
   */
  virtual DofsArrayType& GetDofSet()
  {
    return mDofSet;
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

  // Pointer to the LinearSolver.
  LinearSolverPointerType mpLinearSystemSolver;

  // Container for de system dofs set
  DofsArrayType mDofSet;

  // Flags to set options
  Flags mOptions;

  // Number of degrees of freedom of the problem to be solved
  unsigned int mEquationSystemSize;

  // Level of echo for the builder and solver
  int mEchoLevel;

  // Reactions vector
  SystemVectorPointerType mpReactionsVector;

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

private:

  ///@}
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
  ///@}

}; /// Class SolutionBuilderAndSolver

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_SOLUTION_BUILDER_AND_SOLVER_H_INCLUDED defined

