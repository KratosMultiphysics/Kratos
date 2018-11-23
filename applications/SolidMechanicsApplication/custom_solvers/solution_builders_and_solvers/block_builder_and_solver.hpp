//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               March 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_BLOCK_BUILDER_AND_SOLVER_H_INCLUDED)
#define  KRATOS_BLOCK_BUILDER_AND_SOLVER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_builders_and_solvers/solution_builder_and_solver.hpp"
#include "includes/key_hash.h"

#ifdef USE_GOOGLE_HASH
#include "sparsehash/dense_hash_set" //included in external libraries
#else
#include <unordered_set>
#endif

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


/** @brief Solution Buider and Solver based on block matrix
 *  @details This is builder class the system is not reduced, 1 is set in the diagonal
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class BlockBuilderAndSolver : public SolutionBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
 public:

  ///@name Type Definitions
  ///@{

   /// Pointer definition of BlockBuilderAndSolver
  KRATOS_CLASS_POINTER_DEFINITION(BlockBuilderAndSolver);

  typedef SolutionBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>      BaseType;

  typedef typename BaseType::Pointer                                       BasePointerType;

  typedef typename BaseType::LocalFlagType                                   LocalFlagType;
  typedef typename BaseType::DofsArrayType                                   DofsArrayType;

  typedef typename BaseType::SystemMatrixType                             SystemMatrixType;
  typedef typename BaseType::SystemVectorType                             SystemVectorType;
  typedef typename BaseType::SystemMatrixPointerType               SystemMatrixPointerType;
  typedef typename BaseType::SystemVectorPointerType               SystemVectorPointerType;
  typedef typename BaseType::LocalSystemVectorType                   LocalSystemVectorType;
  typedef typename BaseType::LocalSystemMatrixType                   LocalSystemMatrixType;

  typedef typename ModelPart::NodesContainerType                        NodesContainerType;
  typedef typename ModelPart::ElementsContainerType                  ElementsContainerType;
  typedef typename ModelPart::ConditionsContainerType              ConditionsContainerType;

  typedef typename BaseType::SchemePointerType                           SchemePointerType;
  typedef typename BaseType::LinearSolverPointerType               LinearSolverPointerType;

  struct dof_iterator_hash
  {
    size_t operator()(const Node<3>::DofType::Pointer& it) const
    {
      std::size_t seed = 0;
      HashCombine(seed, it->Id());
      HashCombine(seed, (it->GetVariable()).Key());
      return seed;
    }
  };

  struct dof_iterator_equal
  {
    size_t operator()(const Node<3>::DofType::Pointer& it1, const Node<3>::DofType::Pointer& it2) const
    {
      return (((it1->Id() == it2->Id() && (it1->GetVariable()).Key())==(it2->GetVariable()).Key()));
    }
  };


  ///@}
  ///@name Life Cycle
  ///@{

  /// Default Constructor.
  BlockBuilderAndSolver(LinearSolverPointerType pLinearSystemSolver)
      : BaseType(pLinearSystemSolver)
  {
  }

  /// Destructor.
  ~BlockBuilderAndSolver() override
  {
  }

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
  void BuildLHS(SchemePointerType pScheme,
                ModelPart& rModelPart,
                SystemMatrixType& rA) override
  {
    KRATOS_TRY

    SystemVectorType tmp(rA.size1(), 0.0);
    this->Build(pScheme, rModelPart, rA, tmp);

    KRATOS_CATCH("")
  }

  /**
   * @brief Function to perform the build of the RHS.
   * @details The vector could be sized as the total number of dofs or as the number of non constrained ones
   */
  void BuildRHS(SchemePointerType pScheme,
                ModelPart& rModelPart,
                SystemVectorType& rb) override
  {
    KRATOS_TRY

    BuildRHSNoDirichlet(pScheme,rModelPart,rb);

    const int ndofs = static_cast<int>(this->mDofSet.size());

    //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
#pragma omp parallel for firstprivate(ndofs)
    for (int k = 0; k<ndofs; k++)
    {
      typename DofsArrayType::iterator dof_iterator = this->mDofSet.begin() + k;
      const std::size_t i = dof_iterator->EquationId();

      if (dof_iterator->IsFixed())
        rb[i] = 0.0f;
    }

    KRATOS_CATCH("")

  }

  /**
   * @brief Function to perform the building of the LHS and RHS
   * @details Equivalent (but generally faster) then performing BuildLHS and BuildRHS
   */
  void Build(SchemePointerType pScheme,
             ModelPart& rModelPart,
             SystemMatrixType& rA,
             SystemVectorType& rb) override
  {
    KRATOS_TRY

    if (!pScheme)
      KRATOS_ERROR << "No scheme provided!" << std::endl;

    //getting the elements from the model
    const int nelements = static_cast<int>(rModelPart.Elements().size());

    //getting the array of the conditions
    const int nconditions = static_cast<int>(rModelPart.Conditions().size());

    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
    ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
    ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

    //contributions to the system
    LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
    LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

    //vector containing the localization in the system of the different terms
    Element::EquationIdVectorType EquationId;

    // assemble all elements
    double start_build = OpenMPUtils::GetCurrentTime();

#pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
    {
# pragma omp for  schedule(guided, 512) nowait
      for (int k = 0; k < nelements; k++)
      {
        ModelPart::ElementsContainerType::iterator it = el_begin + k;

        //detect if the element is active or not. If the user did not make any choice the element
        //is active by default
        bool element_is_active = true;
        if ((it)->IsDefined(ACTIVE))
          element_is_active = (it)->Is(ACTIVE);

        if (element_is_active)
        {
          //calculate elemental contribution
          pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, rCurrentProcessInfo);

          //assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
          Assemble(rA, rb, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
#else
          Assemble(rA, rb, LHS_Contribution, RHS_Contribution, EquationId);
#endif
          // clean local elemental memory
          pScheme->Clear(*(it.base()));
        }

      }


      //#pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, EquationId ) schedule(dynamic, 1024)
#pragma omp for  schedule(guided, 512)
      for (int k = 0; k < nconditions; k++)
      {
        ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

        //detect if the condition is active or not. If the user did not make any choice the condition
        //is active by default
        bool condition_is_active = true;
        if ((it)->IsDefined(ACTIVE))
          condition_is_active = (it)->Is(ACTIVE);

        if (condition_is_active)
        {
          //calculate elemental contribution
          pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, rCurrentProcessInfo);

          //assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
          Assemble(rA, rb, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
#else
          Assemble(rA, rb, LHS_Contribution, RHS_Contribution, EquationId);
#endif

          // clean local elemental memory
          pScheme->Clear(*(it.base()));
        }
      }
    }

    double stop_build = OpenMPUtils::GetCurrentTime();
    if (this->mEchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0)
      KRATOS_INFO("parallel_build_time") << stop_build - start_build << std::endl;

    if (this->mEchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0){
      KRATOS_INFO("parallel_build") << "finished" << std::endl;
    }

    KRATOS_CATCH("")

  }

  /**
   * @brief This is a call to the linear system solver
   */
  void SystemSolve(SystemMatrixType& rA,
                   SystemVectorType& rDx,
                   SystemVectorType& rb) override
  {
    KRATOS_TRY

    double norm_b;
    if (TSparseSpace::Size(rb) != 0)
      norm_b = TSparseSpace::TwoNorm(rb);
    else
      norm_b = 0.00;

    if (norm_b != 0.00)
    {
      //do solve
      this->mpLinearSystemSolver->Solve(rA, rDx, rb);
    }
    else
      TSparseSpace::SetToZero(rDx);

    //prints informations about the current time
    if (this->mEchoLevel > 1)
    {
      KRATOS_INFO("linear_solver") << *(this->mpLinearSystemSolver) << std::endl;
    }

    KRATOS_CATCH("")
  }

  /**
   * @brief Function to perform the building and solving phase at the same time.
   * @details It is ideally the fastest and safer function to use when it is possible to solve just after building
   */
  void BuildAndSolve(SchemePointerType pScheme,
                     ModelPart& rModelPart,
                     SystemMatrixType& rA,
                     SystemVectorType& rDx,
                     SystemVectorType& rb) override
  {
    KRATOS_TRY

    double begin_time = OpenMPUtils::GetCurrentTime();
    Build(pScheme, rModelPart, rA, rb);
    double end_time = OpenMPUtils::GetCurrentTime();

    if (this->mEchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)
      KRATOS_INFO("system_build_time") << end_time - begin_time << std::endl;

    ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

    if (this->mEchoLevel == 3)
    {
      KRATOS_INFO("LHS before solve") << "Matrix = " << rA << std::endl;
      KRATOS_INFO("Dx before solve")  << "Solution = " << rDx << std::endl;
      KRATOS_INFO("RHS before solve") << "Vector = " << rb << std::endl;
    }

    begin_time = OpenMPUtils::GetCurrentTime();
    SystemSolveWithPhysics(rA, rDx, rb, rModelPart);
    end_time = OpenMPUtils::GetCurrentTime();


    if (this->mEchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)
      KRATOS_INFO("system_solve_time") << end_time - begin_time << std::endl;

    if (this->mEchoLevel == 3)
    {
      KRATOS_INFO("LHS after solve") << "Matrix = " << rA << std::endl;
      KRATOS_INFO("Dx after solve")  << "Solution = " << rDx << std::endl;
      KRATOS_INFO("RHS after solve") << "Vector = " << rb << std::endl;
    }

    KRATOS_CATCH("")
  }

  /**
   * @brief Function to perform the building of the RHS and solving phase at the same time.
   * @details  It corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
   */
  void BuildRHSAndSolve(SchemePointerType pScheme,
                        ModelPart& rModelPart,
                        SystemMatrixType& rA,
                        SystemVectorType& rDx,
                        SystemVectorType& rb) override
  {
    KRATOS_TRY

    BuildRHS(pScheme, rModelPart, rb);
    SystemSolve(rA, rDx, rb);

    KRATOS_CATCH("")
  }


  /**
   * @brief applies the dirichlet conditions.
   * @details This operation may be very heavy or completely
   * @details unexpensive depending on the implementation choosen and on how the System Matrix
   * @details is built. For explanation of how it works for a particular implementation the user
   * @details should refer to the particular Builder And Solver choosen
   */
  void ApplyDirichletConditions(SchemePointerType pScheme,
                                ModelPart& rModelPart,
                                SystemMatrixType& rA,
                                SystemVectorType& rDx,
                                SystemVectorType& rb) override
  {
    std::size_t system_size = rA.size1();
    std::vector<double> scaling_factors (system_size, 0.0f);

    const int ndofs = static_cast<int>(this->mDofSet.size());

    //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
#pragma omp parallel for firstprivate(ndofs)
    for(int k = 0; k<ndofs; k++)
    {
      typename DofsArrayType::iterator dof_iterator = this->mDofSet.begin() + k;
      if(dof_iterator->IsFixed())
        scaling_factors[k] = 0.0f;
      else
        scaling_factors[k] = 1.0f;

    }

    double* Avalues = rA.value_data().begin();
    std::size_t* Arow_indices = rA.index1_data().begin();
    std::size_t* Acol_indices = rA.index2_data().begin();

    //detect if there is a line of all zeros and set the diagonal to a 1 if this happens
#pragma omp parallel for firstprivate(system_size)
    for(int k = 0; k < static_cast<int>(system_size); ++k)
    {
      std::size_t col_begin = Arow_indices[k];
      std::size_t col_end = Arow_indices[k+1];
      bool empty = true;
      for (std::size_t j = col_begin; j < col_end; ++j)
      {
        if(Avalues[j] != 0.0)
        {
          empty = false;
          break;
        }
      }

      if(empty == true)
      {
        rA(k,k) = 1.0;
        rb[k] = 0.0;
      }
    }

#pragma omp parallel for
    for (int k = 0; k < static_cast<int>(system_size); ++k)
    {
      std::size_t col_begin = Arow_indices[k];
      std::size_t col_end = Arow_indices[k+1];
      double k_factor = scaling_factors[k];
      if (k_factor == 0)
      {
        // zero out the whole row, except the diagonal
        for (std::size_t j = col_begin; j < col_end; ++j)
          if (static_cast<int>(Acol_indices[j]) != k )
            Avalues[j] = 0.0;

        // zero out the RHS
        rb[k] = 0.0;
      }
      else
      {
        // zero out the column which is associated with the zero'ed row
        for (std::size_t j = col_begin; j < col_end; ++j)
          if(scaling_factors[ Acol_indices[j] ] == 0 )
            Avalues[j] = 0.0;
      }
    }
  }

  /**
   * @brief Builds the list of the DofSets involved in the problem by "asking" to each element and condition its Dofs.
   * @details The list of dofs is stores insde the BuilderAndSolver as it is closely connected to the way the matrix and RHS are built
  */
  void SetUpDofSet(SchemePointerType pScheme,
                   ModelPart& rModelPart) override
  {
    KRATOS_TRY

    if( this->mEchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)
    {
      KRATOS_INFO("setting_dofs") << "SetUpDofSet starts" << std::endl;
    }

    //Gets the array of elements from the modeler
    ElementsContainerType& rElements = rModelPart.Elements();
    const int nelements = static_cast<int>(rElements.size());

    Element::DofsVectorType ElementalDofList;

    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

    unsigned int nthreads = OpenMPUtils::GetNumThreads();

#ifdef USE_GOOGLE_HASH
    typedef google::dense_hash_set < Node<3>::DofType::Pointer, dof_iterator_hash>  set_type;
#else
    typedef std::unordered_set < Node<3>::DofType::Pointer, dof_iterator_hash>  set_type;
#endif

    std::vector<set_type> dofs_aux_list(nthreads);

    if( this->mEchoLevel > 2)
    {
      KRATOS_INFO("setting_dofs") << "Number of threads:" << nthreads << std::endl;
    }

    for (int i = 0; i < static_cast<int>(nthreads); i++)
    {
#ifdef USE_GOOGLE_HASH
      dofs_aux_list[i].set_empty_key(Node<3>::DofType::Pointer());
#else
      dofs_aux_list[i].reserve(nelements);
#endif
    }

    if( this->mEchoLevel > 2)
    {
      KRATOS_INFO("setting_dofs") << "initialize_elements" << std::endl;
    }

#pragma omp parallel firstprivate(nelements, ElementalDofList)
    {
#pragma omp for schedule(guided, 512) nowait
      for (int i = 0; i < nelements; i++)
      {
        typename ElementsContainerType::iterator it = rElements.begin() + i;
        const unsigned int this_thread_id = OpenMPUtils::ThisThread();

        // gets list of Dof involved on every element
        pScheme->GetElementalDofList(*(it.base()), ElementalDofList, rCurrentProcessInfo);

        dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
      }

      if( this->mEchoLevel > 2)
      {
        KRATOS_INFO("setting_dofs") << "initialize_conditions" << std::endl;
      }

      ConditionsContainerType& rConditions = rModelPart.Conditions();
      const int nconditions = static_cast<int>(rConditions.size());
#pragma omp for  schedule(guided, 512)
      for (int i = 0; i < nconditions; i++)
      {
        typename ConditionsContainerType::iterator it = rConditions.begin() + i;
        const unsigned int this_thread_id = OpenMPUtils::ThisThread();

        // gets list of Dof involved on every condition
        pScheme->GetConditionDofList(*(it.base()), ElementalDofList, rCurrentProcessInfo);
        dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());

      }
    }

    if( this->mEchoLevel > 2)
    {
      KRATOS_INFO("setting_dofs") << "initialize tree reduction" << std::endl;
    }

    // Here we do a reduction in a tree so to have everything on thread 0
    unsigned int old_max = nthreads;
    unsigned int new_max = ceil(0.5*static_cast<double>(old_max));
    while (new_max>=1 && new_max != old_max)
    {
      if( this->mEchoLevel > 2)
      {
        //just for debugging
        KRATOS_INFO("setting_dofs") << "old_max" << old_max << " new_max:" << new_max << std::endl;
        for (int i = 0; i < static_cast<int>(new_max); i++)
        {
          if (i + new_max < old_max)
          {
            KRATOS_INFO("setting_dofs") << i << " - " << i+new_max << std::endl;
          }
        }
      }

#pragma omp parallel for
      for (int i = 0; i < static_cast<int>(new_max); i++)
      {
        if (i + new_max < old_max)
        {
          dofs_aux_list[i].insert(dofs_aux_list[i+new_max].begin(), dofs_aux_list[i+new_max].end());
          dofs_aux_list[i+new_max].clear();
        }
      }

      old_max = new_max;
      new_max = ceil(0.5*static_cast<double>(old_max));

    }

    if( this->mEchoLevel > 2)
    {
      KRATOS_INFO("setting_dofs") << "initializing ordered array filling" << std::endl;
    }

    DofsArrayType Doftemp;
    this->mDofSet = DofsArrayType();

    Doftemp.reserve(dofs_aux_list[0].size());
    for (auto it= dofs_aux_list[0].begin(); it!= dofs_aux_list[0].end(); it++)
    {
      Doftemp.push_back( it->get() );
    }
    Doftemp.Sort();

    this->mDofSet = Doftemp;

    //Throws an exception if there are no Degrees Of Freedom involved in the analysis
    if (this->mDofSet.size() == 0)
    {
      KRATOS_ERROR << "No degrees of freedom!" << std::endl;
    }
    if( this->mEchoLevel > 2)
    {
      KRATOS_INFO("Dofs size") << this->mDofSet.size() << std::endl;
    }

    if( this->mEchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0)
    {
      KRATOS_INFO("setting_dofs") << "Finished setting up the dofs" << std::endl;
    }

    if( this->mEchoLevel > 2)
    {
      KRATOS_INFO("setting_dofs") << "Initializing lock array" << std::endl;
    }

#ifdef _OPENMP
    if (mlock_array.size() != 0)
    {
      for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
      {
        omp_destroy_lock(&mlock_array[i]);
      }
    }
    mlock_array.resize(this->mDofSet.size());

    for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
    {
      omp_init_lock(&mlock_array[i]);
    }
#endif

    if( this->mEchoLevel > 2)
    {
      KRATOS_INFO("setting_dofs") << "End of setupdofset" << std::endl;
    }

    this->Set(LocalFlagType::DOFS_INITIALIZED, true);

    KRATOS_CATCH("")
  }

  /**
   * @brief organises the dofset in order to speed up the building phase
   */
  void SetUpSystem() override
  {

    //int free_id = 0;
    this->mEquationSystemSize = this->mDofSet.size();
    int ndofs = static_cast<int>(this->mDofSet.size());

#pragma omp parallel for firstprivate(ndofs)
    for (int i = 0; i < static_cast<int>(ndofs); i++)
    {
      typename DofsArrayType::iterator dof_iterator = this->mDofSet.begin() + i;
      dof_iterator->SetEquationId(i);
    }

  }

  /**
   * @brief Resizes and Initializes the system vectors and matrices after SetUpDofSet and SetUpSytem has been called
   */
  void SetUpSystemMatrices(SchemePointerType pScheme,
                           ModelPart& rModelPart,
                           SystemMatrixPointerType& pA,
                           SystemVectorPointerType& pDx,
                           SystemVectorPointerType& pb) override
  {
    KRATOS_TRY

    if (pA == nullptr) //if the pointer is not initialized initialize it to an empty matrix
    {
      SystemMatrixPointerType pNewA = Kratos::make_shared<SystemMatrixType>(0, 0);
      pA.swap(pNewA);
    }
    if (pDx == nullptr) //if the pointer is not initialized initialize it to an empty matrix
    {
      SystemVectorPointerType pNewDx = Kratos::make_shared<SystemVectorType>(0);
      pDx.swap(pNewDx);
    }
    if (pb == nullptr) //if the pointer is not initialized initialize it to an empty matrix
    {
      SystemVectorPointerType pNewb = Kratos::make_shared<SystemVectorType>(0);
      pb.swap(pNewb);
    }

    SystemMatrixType& rA = *pA;
    SystemVectorType& rDx = *pDx;
    SystemVectorType& rb = *pb;

    //resizing the system vectors and matrix
    if (rA.size1() == 0 || this->mOptions.Is(LocalFlagType::REFORM_DOFS)) //if the matrix is not initialized
    {
      rA.resize(this->mEquationSystemSize, this->mEquationSystemSize, false);
      ConstructMatrixStructure(pScheme, rA, rModelPart.Elements(), rModelPart.Conditions(), rModelPart.GetProcessInfo());
    }
    else
    {
      if (rA.size1() != this->mEquationSystemSize || rA.size2() != this->mEquationSystemSize)
      {
        KRATOS_WARNING("block builder resize") << "it should not come here -> this is SLOW" << std::endl;
        rA.resize(this->mEquationSystemSize, this->mEquationSystemSize, true);
        ConstructMatrixStructure(pScheme, rA, rModelPart.Elements(), rModelPart.Conditions(), rModelPart.GetProcessInfo());
      }
    }
    if (rDx.size() != this->mEquationSystemSize)
      rDx.resize(this->mEquationSystemSize, false);
    if (rb.size() != this->mEquationSystemSize)
      rb.resize(this->mEquationSystemSize, false);

    KRATOS_CATCH("")

  }


  /**
   * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
   * @details this function must be called only once per step.
   */
  void InitializeSolutionStep(SchemePointerType pScheme,
                              ModelPart& rModelPart,
                              SystemMatrixPointerType& pA,
                              SystemVectorPointerType& pDx,
                              SystemVectorPointerType& pb) override
  {
    KRATOS_TRY

    BaseType::InitializeSolutionStep(pScheme, rModelPart, pA, pDx, pb);

    // // Initialize
    // pScheme->InitializeSolutionStep(rModelPart);

    // // Set up the system dofs and shape
    // if(this->mOptions.Is(LocalFlagType::REFORM_DOFS))
    //   this->SetSystemDofs(pScheme, rModelPart);

    // // Set up system matrices and vectors
    // double begin_time = OpenMPUtils::GetCurrentTime();
    // this->SetUpSystemMatrices(pScheme, pA, pDx, pb);
    // double end_time = OpenMPUtils::GetCurrentTime();

    // if (this->mEchoLevel >= 2)
    //   KRATOS_INFO("system_resize_time") << ": system_resize_time : " << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;

    KRATOS_CATCH("")
  }

  /**
   * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
   * @details this function must be called only once per step.
   */
  void FinalizeSolutionStep(SchemePointerType pScheme,
                            ModelPart& rModelPart,
                            SystemMatrixPointerType& pA,
                            SystemVectorPointerType& pDx,
                            SystemVectorPointerType& pb) override
  {
    KRATOS_TRY

    BaseType::FinalizeSolutionStep(pScheme, rModelPart, pA, pDx, pb);

    KRATOS_CATCH("")
  }

  /**
   * @brief Calculates system reactions
   * @details A flag controls if reactions must be calculated
   * @details An internal variable to store the reactions vector is needed
   */
  void CalculateReactions(SchemePointerType pScheme,
                          ModelPart& rModelPart,
                          SystemMatrixType& rA,
                          SystemVectorType& rDx,
                          SystemVectorType& rb) override
  {
    TSparseSpace::SetToZero(rb);

    //refresh RHS to have the correct reactions
    BuildRHSNoDirichlet(pScheme, rModelPart, rb);

    const int ndofs = static_cast<int>(this->mDofSet.size());

    //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
#pragma omp parallel for firstprivate(ndofs)
    for (int k = 0; k<ndofs; k++)
    {
      typename DofsArrayType::iterator dof_iterator = this->mDofSet.begin() + k;

      const int i = (dof_iterator)->EquationId();
      if ( (dof_iterator)->IsFixed() ) {
         (dof_iterator)->GetSolutionStepReactionValue() = -rb[i];
      } else{
         (dof_iterator)->GetSolutionStepReactionValue() = 0.0;
      }


    }

  }

  /**
   * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
   */
  void Clear() override
  {

    BaseType::Clear();

#ifdef _OPENMP
    for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
      omp_destroy_lock(&mlock_array[i]);
    mlock_array.resize(0);
#endif

  }

  /**
   * This function is designed to be called once to perform all the checks needed
   * on the input provided. Checks can be "expensive" as the function is designed
   * to catch user's errors.
   * @param rModelPart
   * @return 0 all ok
   */
  int Check(ModelPart& r_mUSE_GOOGLE_HASHodel_part) override
  {
    KRATOS_TRY

    return 0;

    KRATOS_CATCH("");
  }

  ///@}
  ///@name Access
  ///@{

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

#ifdef _OPENMP
  std::vector< omp_lock_t > mlock_array;
#endif

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{

  void SystemSolveWithPhysics(SystemMatrixType& rA,
                              SystemVectorType& rDx,
                              SystemVectorType& rb,
                              ModelPart& rModelPart)
  {
    KRATOS_TRY

    double norm_b;
    if (TSparseSpace::Size(rb) != 0)
      norm_b = TSparseSpace::TwoNorm(rb);
    else
      norm_b = 0.00;

    if (norm_b != 0.00)
    {
      //provide physical data as needed
      if(this->mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
        this->mpLinearSystemSolver->ProvideAdditionalData(rA, rDx, rb, this->mDofSet, rModelPart);

      //do solve
      this->mpLinearSystemSolver->Solve(rA, rDx, rb);
    }
    else
    {
      TSparseSpace::SetToZero(rDx);
      KRATOS_WARNING("RHS") << "ATTENTION! setting the RHS to zero!" << std::endl;
    }

    //prints informations about the current time
    if (this->mEchoLevel > 1)
    {
      KRATOS_INFO("LinearSolver") << *(this->mpLinearSystemSolver) << std::endl;
    }

    KRATOS_CATCH("")

  }

  virtual void ConstructMatrixStructure(SchemePointerType pScheme,
                                        SystemMatrixType& rA,
                                        ElementsContainerType& rElements,
                                        ConditionsContainerType& rConditions,
                                        ProcessInfo& rCurrentProcessInfo)
  {
    //filling with zero the matrix (creating the structure)
    double begin_time = OpenMPUtils::GetCurrentTime();

    const std::size_t equation_size = this->mEquationSystemSize;

#ifdef USE_GOOGLE_HASH
    std::vector<google::dense_hash_set<std::size_t> > indices(equation_size);
    const std::size_t empty_key = 2*equation_size + 10;
#else
    std::vector<std::unordered_set<std::size_t> > indices(equation_size);
#endif

#pragma omp parallel for firstprivate(equation_size)
    for (int iii = 0; iii < static_cast<int>(equation_size); iii++)
    {
#ifdef USE_GOOGLE_HASH
      indices[iii].set_empty_key(empty_key);
#else
      indices[iii].reserve(40);
#endif
    }

    Element::EquationIdVectorType ids(3, 0);

    const int nelements = static_cast<int>(rElements.size());
#pragma omp parallel for firstprivate(nelements, ids)
    for(int iii=0; iii<nelements; iii++)
    {
      typename ElementsContainerType::iterator i_element = rElements.begin() + iii;
      pScheme->EquationId( *(i_element.base()) , ids, rCurrentProcessInfo);
      for (std::size_t i = 0; i < ids.size(); i++)
      {
#ifdef _OPENMP
        omp_set_lock(&mlock_array[ids[i]]);
#endif
        auto& row_indices = indices[ids[i]];
        row_indices.insert(ids.begin(), ids.end());

#ifdef _OPENMP
        omp_unset_lock(&mlock_array[ids[i]]);
#endif
      }

    }

    const int nconditions = static_cast<int>(rConditions.size());
#pragma omp parallel for firstprivate(nconditions, ids)
    for (int iii = 0; iii<nconditions; iii++)
    {
      typename ConditionsContainerType::iterator i_condition = rConditions.begin() + iii;
      pScheme->Condition_EquationId( *(i_condition.base()), ids, rCurrentProcessInfo);
      for (std::size_t i = 0; i < ids.size(); i++)
      {
#ifdef _OPENMP
        omp_set_lock(&mlock_array[ids[i]]);
#endif
        auto& row_indices = indices[ids[i]];
        row_indices.insert(ids.begin(), ids.end());
#ifdef _OPENMP
        omp_unset_lock(&mlock_array[ids[i]]);
#endif
      }
    }

    //count the row sizes
    unsigned int nnz = 0;
    for (unsigned int i = 0; i < indices.size(); i++)
      nnz += indices[i].size();

    rA = boost::numeric::ublas::compressed_matrix<double>(indices.size(), indices.size(), nnz);

    double* Avalues = rA.value_data().begin();
    std::size_t* Arow_indices = rA.index1_data().begin();
    std::size_t* Acol_indices = rA.index2_data().begin();

    //filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
    Arow_indices[0] = 0;
    for (int i = 0; i < static_cast<int>(rA.size1()); i++)
      Arow_indices[i+1] = Arow_indices[i] + indices[i].size();


#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rA.size1()); i++)
    {
      const unsigned int row_begin = Arow_indices[i];
      const unsigned int row_end = Arow_indices[i+1];
      unsigned int k = row_begin;
      for (auto it = indices[i].begin(); it != indices[i].end(); it++)
      {
        Acol_indices[k] = *it;
        Avalues[k] = 0.0;
        k++;
      }

      indices[i].clear(); //deallocating the memory

      std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

    }

    rA.set_filled(indices.size()+1, nnz);

    double end_time = OpenMPUtils::GetCurrentTime();
    if (this->mEchoLevel >= 2)
      KRATOS_INFO("BlockBuilderAndSolver") << "construct matrix structure time:" << end_time - begin_time << "\n" << LoggerMessage::Category::STATISTICS;

  }


  void AssembleLHS(SystemMatrixType& rA,
                   LocalSystemMatrixType& rLHS_Contribution,
                   Element::EquationIdVectorType& rEquationId)
  {
    unsigned int local_size = rLHS_Contribution.size1();

    for (unsigned int i_local = 0; i_local < local_size; i_local++)
    {
      unsigned int i_global = rEquationId[i_local];

      for (unsigned int j_local = 0; j_local < local_size; j_local++)
      {
        unsigned int j_global = rEquationId[j_local];

        rA(i_global, j_global) += rLHS_Contribution(i_local, j_local);
      }
    }

  }

  void AssembleRHS(SystemVectorType& rb,
                   LocalSystemVectorType& rRHS_Contribution,
                   Element::EquationIdVectorType& rEquationId)
  {
    unsigned int local_size = rRHS_Contribution.size();

    for (unsigned int i_local = 0; i_local < local_size; i_local++)
    {
      unsigned int i_global = rEquationId[i_local];

      // ASSEMBLING THE SYSTEM VECTOR
      double& b_value = rb[i_global];
      const double& rhs_value = rRHS_Contribution[i_local];

#pragma omp atomic
      b_value += rhs_value;
    }
  }

  void Assemble(SystemMatrixType& rA,
                SystemVectorType& rb,
                const LocalSystemMatrixType& rLHS_Contribution,
                const LocalSystemVectorType& rRHS_Contribution,
                Element::EquationIdVectorType& rEquationId
#ifdef USE_LOCKS_IN_ASSEMBLY
                ,std::vector< omp_lock_t >& lock_array
#endif
                )
  {
    unsigned int local_size = rLHS_Contribution.size1();

    for (unsigned int i_local = 0; i_local < local_size; i_local++)
    {
      unsigned int i_global = rEquationId[i_local];

#ifdef USE_LOCKS_IN_ASSEMBLY
      omp_set_lock(&lock_array[i_global]);
      b[i_global] += rRHS_Contribution(i_local);
#else
      double& r_a = rb[i_global];
      const double& v_a = rRHS_Contribution(i_local);
#pragma omp atomic
      r_a += v_a;
#endif

      AssembleRowContribution(rA, rLHS_Contribution, i_global, i_local, rEquationId);

#ifdef USE_LOCKS_IN_ASSEMBLY
      omp_unset_lock(&lock_array[i_global]);
#endif
      //note that computation of reactions is not performed here!
    }
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

  inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
  {
    std::vector<std::size_t>::iterator i = v.begin();
    std::vector<std::size_t>::iterator endit = v.end();
    while (i != endit && (*i) != candidate)
    {
      ++i;
    }
    if (i == endit)
    {
      v.push_back(candidate);
    }

  }

  void BuildRHSNoDirichlet(SchemePointerType pScheme,
                           ModelPart& rModelPart,
                           SystemVectorType& rb)
  {
    KRATOS_TRY

    //Getting the Elements
    ElementsContainerType& rElements = rModelPart.Elements();

    //getting the array of the conditions
    ConditionsContainerType& ConditionsArray = rModelPart.Conditions();

    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

    //contributions to the system
    LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
    LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

    //vector containing the localization in the system of the different
    //terms
    Element::EquationIdVectorType EquationId;

    // assemble all elements
    //for (typename ElementsContainerType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)

    const int nelements = static_cast<int>(rElements.size());
#pragma omp parallel firstprivate(nelements, RHS_Contribution, EquationId)
    {
#pragma omp for schedule(guided, 512) nowait
      for(int i=0; i<nelements; i++)
      {
        typename ElementsContainerType::iterator it = rElements.begin() + i;
        //detect if the element is active or not. If the user did not make any choice the element
        //is active by default
        bool element_is_active = true;
        if( (it)->IsDefined(ACTIVE) )
          element_is_active = (it)->Is(ACTIVE);

        if(element_is_active)
        {
          //calculate elemental Right Hand Side Contribution
          pScheme->Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, rCurrentProcessInfo);

          //assemble the elemental contribution
          AssembleRHS(rb, RHS_Contribution, EquationId);
        }
      }

      LHS_Contribution.resize(0, 0, false);
      RHS_Contribution.resize(0, false);

      // assemble all conditions
      //for (typename ConditionsContainerType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
      const int nconditions = static_cast<int>(ConditionsArray.size());
      //#pragma omp parallel for firstprivate(nconditions, RHS_Contribution, EquationId) schedule(dynamic, 1024)
#pragma omp for schedule(guided, 512)
      for (int i = 0; i<nconditions; i++)
      {
        auto it = ConditionsArray.begin() + i;
        //detect if the element is active or not. If the user did not make any choice the element
        //is active by default
        bool condition_is_active = true;
        if( (it)->IsDefined(ACTIVE) )
          condition_is_active = (it)->Is(ACTIVE);

        if(condition_is_active)
        {
          //calculate elemental contribution
          pScheme->Condition_Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, rCurrentProcessInfo);

          //assemble the elemental contribution
          AssembleRHS(rb, RHS_Contribution, EquationId);
        }
      }
    }

    KRATOS_CATCH("")
  }


  inline void CreatePartition(unsigned int number_of_threads,
                              const int number_of_rows,
                              vector<unsigned int>& partitions)
  {
    partitions.resize(number_of_threads + 1);
    int partition_size = number_of_rows / number_of_threads;
    partitions[0] = 0;
    partitions[number_of_threads] = number_of_rows;
    for (unsigned int i = 1; i < number_of_threads; i++)
      partitions[i] = partitions[i - 1] + partition_size;
  }

  inline void AssembleRowContribution(SystemMatrixType& rA,
                                      const Matrix& rAlocal,
                                      const unsigned int i,
                                      const unsigned int i_local,
                                      Element::EquationIdVectorType& rEquationId)
  {
    double* values_vector = rA.value_data().begin();
    std::size_t* index1_vector = rA.index1_data().begin();
    std::size_t* index2_vector = rA.index2_data().begin();

    size_t left_limit = index1_vector[i];
    //	size_t right_limit = index1_vector[i+1];

    //find the first entry
    size_t last_pos = ForwardFind(rEquationId[0],left_limit,index2_vector);
    size_t last_found = rEquationId[0];

#ifndef USE_LOCKS_IN_ASSEMBLY
    double& r_a = values_vector[last_pos];
    const double& v_a = rAlocal(i_local,0);
#pragma omp atomic
    r_a +=  v_a;
#else
    values_vector[last_pos] += rAlocal(i_local,0);
#endif

    //now find all of the other entries
    size_t pos = 0;
    for(unsigned int j=1; j<rEquationId.size(); j++)
    {
      unsigned int id_to_find = rEquationId[j];
      if(id_to_find > last_found)
        pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
      else
        pos = BackwardFind(id_to_find,last_pos-1,index2_vector);

#ifndef USE_LOCKS_IN_ASSEMBLY
      double& r = values_vector[pos];
      const double& v = rAlocal(i_local,j);
#pragma omp atomic
      r +=  v;
#else
      values_vector[pos] += rAlocal(i_local,j);
#endif

      last_found = id_to_find;
      last_pos = pos;
    }
  }

  inline unsigned int ForwardFind(const unsigned int id_to_find,
                                  const unsigned int start,
                                  const size_t* index_vector)
  {
    unsigned int pos = start;
    while(id_to_find != index_vector[pos]) pos++;
    return pos;
  }

  inline unsigned int BackwardFind(const unsigned int id_to_find,
                                   const unsigned int start,
                                   const size_t* index_vector)
  {
    unsigned int pos = start;
    while(id_to_find != index_vector[pos]) pos--;
    return pos;
  }

  ///@}
  ///@name Private  Access
  ///@{

  ///@}
  ///@name Serialization
  ///@{

  ///@}
  ///@name Private Inquiry
  ///@{

  ///@}
  ///@name Un accessible methods
  ///@{

  ///@}

}; // Class BlockBuilderAndSolver
///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos

#endif // KRATOS_BLOCK_BUILDER_AND_SOLVER_H_INCLUDED  defined
