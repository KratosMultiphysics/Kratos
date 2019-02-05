//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:   michael.andre@tum.de $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:           September 2016 $
//   Revision:            $Revision:           March 2018 $
//
//

#if !defined(KRATOS_EIGENSOLVER_STRATEGY_H_INCLUDED)
#define  KRATOS_EIGENSOLVER_STRATEGY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_strategies/solution_strategy.hpp"

// Application includes
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

/// Strategy for solving generalized eigenvalue problems.
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class EigensolverStrategy
    : public SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
 public:
  ///@name Type Definitions
  ///@{

  KRATOS_CLASS_POINTER_DEFINITION(EigensolverStrategy);

  typedef SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>          BaseType;

  typedef typename BaseType::LocalFlagType                               LocalFlagType;

  typedef typename BaseType::BuilderAndSolverType                 BuilderAndSolverType;

  typedef typename BaseType::SchemeType                                     SchemeType;

  typedef typename TDenseSpace::VectorType                             DenseVectorType;

  typedef typename TDenseSpace::MatrixType                             DenseMatrixType;

  typedef TSparseSpace SparseSpaceType;

  typedef typename TSparseSpace::VectorPointerType             SparseVectorPointerType;

  typedef typename TSparseSpace::MatrixPointerType             SparseMatrixPointerType;

  typedef typename TSparseSpace::MatrixType                           SparseMatrixType;

  typedef typename TSparseSpace::VectorType                           SparseVectorType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Constructor.
  EigensolverStrategy(ModelPart& rModelPart,
                      typename SchemeType::Pointer pScheme,
                      typename BuilderAndSolverType::Pointer pBuilderAndSolver,
                      Flags & rOptions,
                      bool ComputeModalContribution = false)
      : SolutionStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, rOptions)
  {
    KRATOS_TRY

    mpScheme = pScheme;

    mpBuilderAndSolver = pBuilderAndSolver;

    // set if modal contribution is computed
    mComputeModalContribution = ComputeModalContribution;

    mpMassMatrix = Kratos::make_shared<SparseMatrixType>(0,0);
    mpStiffnessMatrix = Kratos::make_shared<SparseMatrixType>(0,0);

    mpBuilderAndSolver->SetEchoLevel(this->mEchoLevel);

    KRATOS_CATCH("")
  }

  /// Deleted copy constructor.
  EigensolverStrategy(const EigensolverStrategy& Other) = delete;

  /// Destructor.
  ~EigensolverStrategy() override
  {
    // Clear() controls order of deallocation to avoid invalid memory access in some special cases.
    // warning: BaseType::GetModelPart() may be invalid here.
    this->Clear();
  }

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  /**
   * Performs all the required operations that should be done (for each step)
   * before solving the solution step.
   * A member variable should be used as a flag to make sure this function is called only once per step.
   */
  void InitializeSolutionStep() override
  {
    KRATOS_TRY

    if(this->IsNot(LocalFlagType::INITIALIZED))
      this->Initialize();

    SparseVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer(); //dummy
    SparseVectorPointerType pb  = SparseSpaceType::CreateEmptyVectorPointer(); //dummy

    //set up the system
    this->SetSystemDofs();

    //set system sizes
    this->SetSystemMatrices(pDx,pb);

    //initial operations ... things that are constant over the Solution Step
    mpBuilderAndSolver->InitializeSolutionStep(mpScheme, BaseType::GetModelPart(), mpStiffnessMatrix, pDx, pb);

    //initial operations ... things that are constant over the Solution Step
    mpScheme->InitializeSolutionStep(BaseType::GetModelPart());


    KRATOS_CATCH("")
  }

  /**
   * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
   */
  bool SolveSolutionStep() override
  {
    KRATOS_TRY

    // Initialize dummy rhs vector
    SparseVectorPointerType pb  = SparseSpaceType::CreateEmptyVectorPointer(); //dummy
    SparseSpaceType::Resize((*pb),SparseSpaceType::Size1((*mpMassMatrix)));
    SparseSpaceType::Set((*pb),0.0);

    // Generate lhs matrix. the factor 1 is chosen to preserve
    // SPD property
    this->GetModelPart().GetProcessInfo()[BUILD_LEVEL] = 1;
    TSparseSpace::SetToZero((*mpMassMatrix));
    mpBuilderAndSolver->Build(mpScheme,this->GetModelPart(),(*mpMassMatrix),(*pb));
    this->ApplyDirichletConditions((*mpMassMatrix), 1.0);

    // Generate rhs matrix. the factor -1 is chosen to make
    // Eigenvalues corresponding to fixed dofs negative
    this->GetModelPart().GetProcessInfo()[BUILD_LEVEL] = 2;
    TSparseSpace::SetToZero((*mpStiffnessMatrix));
    mpBuilderAndSolver->Build(mpScheme,this->GetModelPart(),(*mpStiffnessMatrix),(*pb));
    ApplyDirichletConditions((*mpStiffnessMatrix), -1.0);

    // Eigenvector matrix and eigenvalue vector are initialized by the solver
    DenseVectorType Eigenvalues;
    DenseMatrixType Eigenvectors;

    // Solve for eigenvalues and eigenvectors
    mpBuilderAndSolver->GetLinearSystemSolver()->Solve((*mpStiffnessMatrix),
                                                       (*mpMassMatrix),
                                                       Eigenvalues,
                                                       Eigenvectors);

    this->AssignVariables(Eigenvalues,Eigenvectors);

    if( mComputeModalContribution == true )
      this->ComputeModalContribution((*mpMassMatrix),Eigenvalues,Eigenvectors);

    return true;

    KRATOS_CATCH("")
  }

  /**
   * Clears the internal storage
   */
  void Clear() override
  {
    KRATOS_TRY

    // if the preconditioner is saved between solves, it should be cleared here
    mpBuilderAndSolver->GetLinearSystemSolver()->Clear();

    if(this->pGetMassMatrix() != nullptr)
    {
      this->pGetMassMatrix() = nullptr;
    }

    if(this->pGetStiffnessMatrix() != nullptr)
    {
      this->pGetStiffnessMatrix() = nullptr;
    }

    mpBuilderAndSolver->Clear();

    mpScheme->Clear();

    KRATOS_CATCH("")
  }

  /**
   * Function to perform expensive checks.
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

    //check internal variable build level
    KRATOS_CHECK_VARIABLE_KEY(BUILD_LEVEL);
    KRATOS_CHECK_VARIABLE_KEY(EIGENVALUE_VECTOR);
    KRATOS_CHECK_VARIABLE_KEY(EIGENVECTOR_MATRIX);

    return 0;

    KRATOS_CATCH("")
  }

  ///@}
  ///@name Access
  ///@{

  /* @brief This sets the level of echo for the solving strategy
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


  void SetScheme(typename SchemeType::Pointer pScheme)
  {
    mpScheme = pScheme;
  };

  typename SchemeType::Pointer& pGetScheme()
  {
    return mpScheme;
  };

  void SetBuilderAndSolver(typename BuilderAndSolverType::Pointer pBuilderAndSolver)
  {
    mpBuilderAndSolver = pBuilderAndSolver;
  };

  typename BuilderAndSolverType::Pointer& pGetBuilderAndSolver()
  {
    return mpBuilderAndSolver;
  };

  SparseMatrixType& GetMassMatrix()
  {
    return *mpMassMatrix;
  }

  SparseMatrixType& GetStiffnessMatrix()
  {
    return *mpStiffnessMatrix;
  }

  SparseMatrixPointerType& pGetMassMatrix()
  {
    return mpMassMatrix;
  }

  SparseMatrixPointerType& pGetStiffnessMatrix()
  {
    return mpStiffnessMatrix;
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

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{

  /**
   * Initialization to be performed once before using the strategy.
   */
  void Initialize() override
  {
    KRATOS_TRY

    //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
    if(mpScheme->IsNot(LocalFlagType::INITIALIZED))
      mpScheme->Initialize(this->GetModelPart());

    //set up the system
    this->SetSystemDofs();

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

  typename SchemeType::Pointer                     mpScheme;

  typename BuilderAndSolverType::Pointer mpBuilderAndSolver;

  SparseMatrixPointerType                      mpMassMatrix;

  SparseMatrixPointerType                 mpStiffnessMatrix;

  bool                            mComputeModalContribution;

  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  /**
   * @brief Performs all the required operations to resize system matrices and vectors
   */
  void SetSystemMatrices(SparseVectorPointerType& pDx, SparseVectorPointerType& pb)
  {
    KRATOS_TRY

    // Mass matrix
    mpBuilderAndSolver->SetUpSystemMatrices(mpScheme, this->GetModelPart(), mpMassMatrix, pDx, pb);

    // Stiffness matrix
    mpBuilderAndSolver->SetUpSystemMatrices(mpScheme, this->GetModelPart(), mpStiffnessMatrix, pDx, pb);

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
      KRATOS_INFO("setup_system_time") << ": setup_system_time : " << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;

    KRATOS_CATCH("")
  }

  /// Apply Dirichlet boundary conditions without modifying dof pattern.
  /**
   *  The dof pattern is preserved to support algebraic multigrid solvers with
   *  component-wise aggregation. Rows and columns of the fixed dofs are replaced
   *  with zeros on the off-diagonal and the diagonal is scaled by factor.
   */
  void ApplyDirichletConditions(SparseMatrixType& rA, double Factor)
  {
    KRATOS_TRY

    const std::size_t SystemSize = rA.size1();
    std::vector<double> ScalingFactors(SystemSize);
    auto& rDofSet = this->pGetBuilderAndSolver()->GetDofSet();
    const int NumDofs = static_cast<int>(rDofSet.size());

    // NOTE: dofs are assumed to be numbered consecutively
#pragma omp parallel for firstprivate(NumDofs)
    for(int k = 0; k<NumDofs; k++)
    {
      auto dof_iterator = std::begin(rDofSet) + k;
      ScalingFactors[k] = (dof_iterator->IsFixed()) ? 0.0 : 1.0;
    }

    double* AValues = std::begin(rA.value_data());
    std::size_t* ARowIndices = std::begin(rA.index1_data());
    std::size_t* AColIndices = std::begin(rA.index2_data());

    // if there is a line of all zeros, put one on the diagonal
    // #pragma omp parallel for firstprivate(SystemSize)
    // for(int k = 0; k < static_cast<int>(SystemSize); ++k)
    // {
    //     std::size_t ColBegin = ARowIndices[k];
    //     std::size_t ColEnd = ARowIndices[k+1];
    //     bool empty = true;
    //     for (auto j = ColBegin; j < ColEnd; ++j)
    //         if(AValues[j] != 0.0)
    //         {
    //             empty = false;
    //             break;
    //         }
    //     if(empty == true)
    //         rA(k,k) = 1.0;
    // }

#pragma omp parallel for
    for (int k = 0; k < static_cast<int>(SystemSize); ++k)
    {
      std::size_t ColBegin = ARowIndices[k];
      std::size_t ColEnd = ARowIndices[k+1];
      if(ScalingFactors[k] == 0.0)
      {
        // row dof is fixed. zero off-diagonal columns and factor diagonal
        for (std::size_t j = ColBegin; j < ColEnd; ++j)
        {
          if(static_cast<int>(AColIndices[j]) != k)
          {
            AValues[j] = 0.0;
          }
          else
          {
            AValues[j] *= Factor;
          }
        }
      }
      else
      {
        // row dof is not fixed. zero columns associated with fixed dofs
        for (std::size_t j = ColBegin; j < ColEnd; ++j)
        {
          AValues[j] *= ScalingFactors[AColIndices[j]];
        }
      }
    }

    KRATOS_CATCH("")
  }

  /// Assign eigenvalues and eigenvectors to kratos variables.
  void AssignVariables(DenseVectorType& rEigenvalues, DenseMatrixType& rEigenvectors)
  {
    KRATOS_TRY

    const std::size_t NumEigenvalues = rEigenvalues.size();

    // store eigenvalues in process info
    this->GetModelPart().GetProcessInfo()[EIGENVALUE_VECTOR] = rEigenvalues;

    for (ModelPart::NodeIterator itNode = this->GetModelPart().NodesBegin(); itNode!= this->GetModelPart().NodesEnd(); itNode++)
    {
      ModelPart::NodeType::DofsContainerType& NodeDofs = itNode->GetDofs();
      const std::size_t NumNodeDofs = NodeDofs.size();
      Matrix& rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);
      if(rNodeEigenvectors.size1() != NumEigenvalues || rNodeEigenvectors.size2() != NumNodeDofs)
      {
        rNodeEigenvectors.resize(NumEigenvalues,NumNodeDofs,false);
      }

      // the jth column index of EIGENVECTOR_MATRIX corresponds to the jth nodal dof. therefore,
      // the dof ordering must not change.
      if(NodeDofs.IsSorted() == false)
      {
        NodeDofs.Sort();
      }

      // fill the EIGENVECTOR_MATRIX
      for (std::size_t i = 0; i < NumEigenvalues; i++)
        for (std::size_t j = 0; j < NumNodeDofs; j++)
        {
          auto itDof = std::begin(NodeDofs) + j;
          rNodeEigenvectors(i,j) = rEigenvectors(i,itDof->EquationId());
        }
    }

    KRATOS_CATCH("")
  }

  /// Compute nodal contributions
  void ComputeModalContribution(SparseMatrixType& rMassMatrix, DenseVectorType& rEigenValues, DenseMatrixType& rEigenVectors)
  {
    KRATOS_TRY

    //Computing modal contribution
    const auto num_eigen_values = rEigenValues.size();
    const auto system_size = rMassMatrix.size1();
    Matrix mass(num_eigen_values,num_eigen_values);
    noalias(mass) = ZeroMatrix(num_eigen_values,num_eigen_values);
    Vector mode_contribution(num_eigen_values);
    noalias(mode_contribution)= ZeroVector(num_eigen_values);
    Vector ratio_mass_mode_contribution(num_eigen_values);
    noalias(ratio_mass_mode_contribution) = ZeroVector(num_eigen_values);
    Matrix eigen_contribution(num_eigen_values, system_size);
    noalias(eigen_contribution)= ZeroMatrix(num_eigen_values, system_size);

    double total_mass= 0.0;
    for (std::size_t i = 0; i < system_size; i++)
    {
      for (std::size_t j = 0; j < system_size; j++)
      {
        total_mass += rMassMatrix(i,j);
      }
    }

    noalias(eigen_contribution) = prod(rEigenVectors,rMassMatrix);
    noalias(mass) = prod(eigen_contribution,trans(rEigenVectors));
    double total_mass_contribution =0.0;

    for (std::size_t i = 0; i < num_eigen_values; i++)
    {
      for (std::size_t j = 0; j < system_size; j++)
      {
        mode_contribution[i] += eigen_contribution(i,j);
      }

      ratio_mass_mode_contribution[i] = (mode_contribution[i]*mode_contribution[i])/(mass(i,i)*total_mass)*100.0;
      total_mass_contribution += ratio_mass_mode_contribution[i];
    }

    if (this->mEchoLevel >= 1){
      KRATOS_INFO(" eigen_ratio") << " ::EIGEN_CONTRIBUTION:: (Mode/Mass) RATIO [" <<ratio_mass_mode_contribution << "]\n" << LoggerMessage::Category::STATUS;
      KRATOS_INFO(" eigen_total") << " ::EIGEN_CONTRIBUTION:: (Mode/Mass) TOTAL [" << total_mass_contribution<< "]\n" << LoggerMessage::Category::STATUS;
    }

    KRATOS_CATCH("")
  }


  ///@}
  ///@name Private  Access
  ///@{

  ///@}
  ///@name Private Inquiry
  ///@{

  ///@}

}; /// Class EigensolverStrategy

///@}

///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_EIGENSOLVER_STRATEGY_H_INCLUDED  defined
