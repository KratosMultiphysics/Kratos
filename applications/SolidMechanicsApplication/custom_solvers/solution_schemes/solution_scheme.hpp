//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               March 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SOLUTION_SCHEME_H_INCLUDED)
#define KRATOS_SOLUTION_SCHEME_H_INCLUDED

// System includes
//#include <set>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/model_part.h"
#include "custom_processes/solver_process.hpp"
#include "custom_solvers/time_integration_methods/time_integration_method.hpp"

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

/** @brief Solution scheme base class
 *  @details This is the base class for the schemes
 */

template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class SolutionScheme : public Flags
{
 public:

  ///@name Type Definitions
  ///@{
  typedef SolutionScheme<TSparseSpace,TDenseSpace>              SolutionSchemeType;
  typedef typename SolutionSchemeType::Pointer           SolutionSchemePointerType;
  typedef SolverLocalFlags                                           LocalFlagType;
  typedef typename ModelPart::DofsArrayType                          DofsArrayType;

  typedef typename TSparseSpace::MatrixType                       SystemMatrixType;
  typedef typename TSparseSpace::VectorType                       SystemVectorType;
  typedef typename TDenseSpace::MatrixType                   LocalSystemMatrixType;
  typedef typename TDenseSpace::VectorType                   LocalSystemVectorType;

  typedef ModelPart::NodesContainerType                         NodesContainerType;
  typedef ModelPart::ElementsContainerType                   ElementsContainerType;
  typedef ModelPart::ConditionsContainerType               ConditionsContainerType;

  typedef ModelPart::NodeType                                                NodeType;
  typedef array_1d<double, 3>                                              VectorType;
  typedef Variable<VectorType>                                     VariableVectorType;
  typedef TimeIntegrationMethod<VariableVectorType,VectorType>  IntegrationVectorType;
  typedef typename IntegrationVectorType::Pointer        IntegrationVectorPointerType;
  typedef std::vector<IntegrationVectorPointerType>      IntegrationMethodsVectorType;

  typedef Variable<double>                                         VariableScalarType;
  typedef TimeIntegrationMethod<VariableScalarType, double>     IntegrationScalarType;
  typedef typename IntegrationScalarType::Pointer        IntegrationScalarPointerType;
  typedef std::vector<IntegrationScalarPointerType>      IntegrationMethodsScalarType;

  typedef typename SolverProcess::Pointer                          ProcessPointerType;
  typedef std::vector<ProcessPointerType>                    ProcessPointerVectorType;

  /// Pointer definition of SolutionScheme
  KRATOS_CLASS_POINTER_DEFINITION(SolutionScheme);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default Constructor.
  SolutionScheme() : Flags() {SetDefaultFlags();}

  /// Constructor.
  SolutionScheme(Flags& rOptions) : Flags(), mOptions(rOptions) {SetDefaultFlags();}

  /// Constructor.
  SolutionScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods, Flags& rOptions) : Flags(), mOptions(rOptions), mTimeVectorIntegrationMethods(rTimeVectorIntegrationMethods) {SetDefaultFlags();}

  /// Constructor.
  SolutionScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods) : Flags(), mTimeVectorIntegrationMethods(rTimeVectorIntegrationMethods) {SetDefaultFlags();}

  /// Constructor.
  SolutionScheme(IntegrationMethodsScalarType& rTimeScalarIntegrationMethods, Flags& rOptions) : Flags(), mOptions(rOptions), mTimeScalarIntegrationMethods(rTimeScalarIntegrationMethods) {SetDefaultFlags();}

  /// Constructor.
  SolutionScheme(IntegrationMethodsScalarType& rTimeScalarIntegrationMethods) : Flags(), mTimeScalarIntegrationMethods(rTimeScalarIntegrationMethods) {SetDefaultFlags();}

  /// Constructor.
  SolutionScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods,
                 IntegrationMethodsScalarType& rTimeScalarIntegrationMethods,
                 Flags& rOptions)
      : Flags(), mOptions(rOptions), mTimeVectorIntegrationMethods(rTimeVectorIntegrationMethods), mTimeScalarIntegrationMethods(rTimeScalarIntegrationMethods) {SetDefaultFlags();}

  /// Constructor.
  SolutionScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods,
                 IntegrationMethodsScalarType& rTimeScalarIntegrationMethods)
      : Flags(), mTimeVectorIntegrationMethods(rTimeVectorIntegrationMethods), mTimeScalarIntegrationMethods(rTimeScalarIntegrationMethods) {SetDefaultFlags();}

  /// Copy contructor.
  SolutionScheme(SolutionScheme& rOther) : mOptions(rOther.mOptions)
                                         , mProcesses(rOther.mProcesses)
  {
    std::copy(std::begin(rOther.mTimeVectorIntegrationMethods), std::end(rOther.mTimeVectorIntegrationMethods), std::back_inserter(mTimeVectorIntegrationMethods));
  }

  /// Clone.
  virtual SolutionSchemePointerType Clone()
  {
    return SolutionSchemePointerType( new SolutionScheme(*this) );
  }

  /// Destructor.
  ~SolutionScheme() override {}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  /**
   * @brief SetDefaultSchemeFlags.
   * @details This is intended to be called with the constructor
   */
  void SetDefaultFlags()
  {
    KRATOS_TRY

    if( this->mOptions.IsNotDefined(LocalFlagType::MOVE_MESH) )
      mOptions.Set(LocalFlagType::MOVE_MESH,true); //default : lagrangian mesh update

    if( this->mOptions.IsNotDefined(LocalFlagType::UPDATE_VARIABLES) )
      mOptions.Set(LocalFlagType::UPDATE_VARIABLES,true); //default : derivatives update

    if( this->mOptions.IsNotDefined(LocalFlagType::INCREMENTAL_SOLUTION) )
      mOptions.Set(LocalFlagType::INCREMENTAL_SOLUTION,true); //default : dof is the variable increment

    KRATOS_CATCH("")
  }

  /**
   * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
   * @details This is intended to be called just once when the strategy is initialized
   */
  virtual void Initialize(ModelPart& rModelPart)
  {
    KRATOS_TRY

    for(typename IntegrationMethodsVectorType::iterator it=mTimeVectorIntegrationMethods.begin();
        it!=mTimeVectorIntegrationMethods.end(); ++it)
      (*it)->SetParameters(rModelPart.GetProcessInfo());

    for(typename IntegrationMethodsScalarType::iterator it=mTimeScalarIntegrationMethods.begin();
        it!=mTimeScalarIntegrationMethods.end(); ++it)
      (*it)->SetParameters(rModelPart.GetProcessInfo());

    this->InitializeElements(rModelPart);

    this->InitializeConditions(rModelPart);

    this->Set(LocalFlagType::INITIALIZED, true);

    KRATOS_CATCH("")
  }

  /**
   * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
   * @details this function must be called only once per step.
   */
  virtual void InitializeSolutionStep(ModelPart& rModelPart)
  {
    KRATOS_TRY

    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

    for(typename ProcessPointerVectorType::iterator it=mProcesses.begin(); it!=mProcesses.end(); ++it)
      (*it)->ExecuteInitializeSolutionStep();

#pragma omp parallel for
    for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++)
    {
      auto itElem = rModelPart.ElementsBegin() + i;
      itElem->InitializeSolutionStep(rCurrentProcessInfo);
    }

#pragma omp parallel for
    for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++)
    {
      auto itCond = rModelPart.ConditionsBegin() + i;
      itCond->InitializeSolutionStep(rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
  }


  /**
   * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
   * @details this function must be called only once per step.
   */

  virtual void FinalizeSolutionStep(ModelPart& rModelPart)
  {
    KRATOS_TRY

    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

    for(typename ProcessPointerVectorType::iterator it=mProcesses.begin(); it!=mProcesses.end(); ++it)
      (*it)->ExecuteFinalizeSolutionStep();


#pragma omp parallel for
    for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++)
    {
      auto itElem = rModelPart.ElementsBegin() + i;
      itElem->FinalizeSolutionStep(rCurrentProcessInfo);
    }


#pragma omp parallel for
    for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++)
    {
      auto itCond = rModelPart.ConditionsBegin() + i;
      itCond->FinalizeSolutionStep(rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
  }

  /**
   * @brief Performs all the required operations that should be done (for each iteration) before solving a solution iteration.
   * @details this function must be called only once per iteration
   */
  virtual void InitializeNonLinearIteration(ModelPart& rModelPart)
  {
    KRATOS_TRY

    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

    for(typename ProcessPointerVectorType::iterator it=mProcesses.begin(); it!=mProcesses.end(); ++it)
      (*it)->ExecuteInitialize(); //corresponds to ExecuteInitializeNonLinearIteration()

#pragma omp parallel for
    for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++)
    {
      auto itElem = rModelPart.ElementsBegin() + i;
      itElem->InitializeNonLinearIteration(rCurrentProcessInfo);
    }


#pragma omp parallel for
    for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++)
    {
      auto itCond = rModelPart.ConditionsBegin() + i;
      itCond->InitializeNonLinearIteration(rCurrentProcessInfo);
    }


    KRATOS_CATCH("")
  }


  /**
   * @brief Performs all the required operations that should be done (for each iteration) after solving a solution iteration.
   * @details this function must be called only once per iteration
   */
  virtual void FinalizeNonLinearIteration(ModelPart& rModelPart)
  {
    KRATOS_TRY

    for(typename ProcessPointerVectorType::iterator it=mProcesses.begin(); it!=mProcesses.end(); ++it)
      (*it)->ExecuteFinalize(); //corresponds to ExecuteFinalizeNonLinearIteration()

    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

#pragma omp parallel for
    for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++)
    {
      auto itElem = rModelPart.ElementsBegin() + i;
      itElem->FinalizeNonLinearIteration(rCurrentProcessInfo);
    }


#pragma omp parallel for
    for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++)
    {
      auto itCond = rModelPart.ConditionsBegin() + i;
      itCond->FinalizeNonLinearIteration(rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
  }


  /**
   * @brief Performing the prediction of the solution.
   * @details this function must be called only once per step.
   */
  virtual void Predict(ModelPart& rModelPart,
                       DofsArrayType& rDofSet,
                       SystemVectorType& rDx)
  {
    KRATOS_TRY

    KRATOS_CATCH("")
  }

  /**
   * @brief Performing the update of the solution.
   * @details this function must be called only once per iteration
   */
  virtual void Update(ModelPart& rModelPart,
                      DofsArrayType& rDofSet,
                      SystemVectorType& rDx)
  {
    KRATOS_TRY


    KRATOS_CATCH("")
  }


  /**
   * @brief Performing the update of the solution Dofs (total solution)
   * @details this function must be called only once per iteration
   */
  static inline void SetSolution(ModelPart& rModelPart,
                                 DofsArrayType& rDofSet,
                                 SystemVectorType& rDx)
  {
    KRATOS_TRY

    const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

    // Update of displacement (by DOF)
    OpenMPUtils::PartitionVector DofPartition;
    OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofPartition);

    const int ndof = static_cast<int>(rDofSet.size());
    typename DofsArrayType::iterator DofBegin = rDofSet.begin();

    #pragma omp parallel for firstprivate(DofBegin)
    for(int i = 0;  i < ndof; i++)
    {
      typename DofsArrayType::iterator itDof = DofBegin + i;

      if (itDof->IsFree() )
      {
        itDof->GetSolutionStepValue() = TSparseSpace::GetValue(rDx,itDof->EquationId());
      }
    }

    KRATOS_CATCH("")
  }


  /**
   * @brief Performing the update of the solution Dofs (incremental solution)
   * @details this function must be called only once per iteration
   */
  static inline void AddSolution(ModelPart& rModelPart,
                                 DofsArrayType& rDofSet,
                                 SystemVectorType& rDx)
  {
    KRATOS_TRY

    const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

    // Update of displacement (by DOF)
    OpenMPUtils::PartitionVector DofPartition;
    OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofPartition);

    const int ndof = static_cast<int>(rDofSet.size());
    typename DofsArrayType::iterator DofBegin = rDofSet.begin();

    #pragma omp parallel for firstprivate(DofBegin)
    for(int i = 0;  i < ndof; i++)
    {
      typename DofsArrayType::iterator itDof = DofBegin + i;

      if (itDof->IsFree() )
      {
        itDof->GetSolutionStepValue() += TSparseSpace::GetValue(rDx,itDof->EquationId());
      }
    }

    KRATOS_CATCH("")
  }

  /**
   * @brief Performing the update of the solution Dofs
   * @details this function must be called only once per iteration
   */
  virtual void UpdateDofs(ModelPart& rModelPart,
                          DofsArrayType& rDofSet,
                          SystemVectorType& rDx)
  {
    KRATOS_TRY

    if( mOptions.Is(LocalFlagType::INCREMENTAL_SOLUTION) )
      AddSolution(rModelPart,rDofSet,rDx);  //dof = incremental variable
    else
      SetSolution(rModelPart,rDofSet,rDx);  //dof = total variable

    KRATOS_CATCH("")
  }



  /**
   * @brief Performing the update of the solution variables
   * @details this function must be called only once per iteration
   */
  virtual void UpdateVariables(ModelPart& rModelPart)
  {
    KRATOS_TRY

    if( this->mOptions.Is(LocalFlagType::UPDATE_VARIABLES) ){

      // Updating time derivatives (nodally for efficiency)
      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

      #pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
      {
        NodesContainerType::iterator itNode = NodeBegin + i;

        this->IntegrationMethodUpdate(*itNode);
      }
    }

    KRATOS_CATCH("")
  }

  /**
   * @brief Performing the prediction of the solution variables
   * @details this function must be called only once per iteration
   */
  virtual void PredictVariables(ModelPart& rModelPart)
  {
    KRATOS_TRY

    // Updating time derivatives (nodally for efficiency)
    const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector NodePartition;
    OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

    const int nnodes = static_cast<int>(rModelPart.Nodes().size());
    NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

    #pragma omp parallel for firstprivate(NodeBegin)
    for(int i = 0;  i < nnodes; i++)
    {
      NodesContainerType::iterator itNode = NodeBegin + i;

      this->IntegrationMethodPredict(*itNode);
    }

    KRATOS_CATCH("")
  }


  /**
   * @brief This function is designed to move the mesh
   * @note Be careful it just consider displacements, derive this method to adapt to your own strategies (ALE, FSI, etc...)
   */
  virtual void MoveMesh(ModelPart& rModelPart)
  {
    KRATOS_TRY

    if( this->mOptions.Is(LocalFlagType::MOVE_MESH) ){

      if (rModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT_X) == false)
      {
        KRATOS_ERROR << "It is impossible to move the mesh since the DISPLACEMENT variable is not in the Model Part. Add DISPLACEMENT to the list of variables" << std::endl;
      }

      bool DisplacementIntegration = false;
      for(typename IntegrationMethodsVectorType::iterator it=mTimeVectorIntegrationMethods.begin();
          it!=mTimeVectorIntegrationMethods.end(); ++it)
      {
        if( "DISPLACEMENT" == (*it)->GetVariableName() ){
          DisplacementIntegration = true;
          break;
        }
      }

      if(DisplacementIntegration == true){

        // Update mesh positions : node coordinates
        const int nnodes = rModelPart.NumberOfNodes();
        ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();

#pragma omp parallel for
        for(int i = 0; i<nnodes; i++)
        {
          ModelPart::NodesContainerType::iterator it_node = it_begin + i;

          noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
          noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT);
        }

      }

    }

    KRATOS_CATCH("")
  }

  /**
   * @brief Liberates internal storage for an element
   */
  virtual void Clear(Element::Pointer rCurrentElement)
  {
    rCurrentElement->CleanMemory();
  }

  /**
   * @brief Liberates internal storage for a condition
   */
  virtual void Clear(Condition::Pointer rCurrentCondition)
  {
    rCurrentCondition->CleanMemory();
  }

  /**
   * @brief Liberates internal storage.
   */
  virtual void Clear() {}

  /**
   * @brief This function is designed to be called once to perform all the checks needed
   * @return 0 all ok
   */
  virtual int Check(ModelPart& rModelPart)
  {
    KRATOS_TRY

    for(ModelPart::ElementsContainerType::iterator it=rModelPart.ElementsBegin();
        it!=rModelPart.ElementsEnd(); ++it)
    {
      it->Check(rModelPart.GetProcessInfo());
    }

    for(ModelPart::ConditionsContainerType::iterator it=rModelPart.ConditionsBegin();
        it!=rModelPart.ConditionsEnd(); ++it)
    {
      it->Check(rModelPart.GetProcessInfo());
    }

    for(typename IntegrationMethodsVectorType::iterator it=mTimeVectorIntegrationMethods.begin();
        it!=mTimeVectorIntegrationMethods.end(); ++it)
      (*it)->Check(rModelPart.GetProcessInfo());

    for(typename IntegrationMethodsScalarType::iterator it=mTimeScalarIntegrationMethods.begin();
        it!=mTimeScalarIntegrationMethods.end(); ++it)
      (*it)->Check(rModelPart.GetProcessInfo());

    KRATOS_CATCH("")

    return 0;
  }

  /** these functions is designed to be called in the builder and solver to introduce
      the selected time integration scheme. It "asks" the matrix needed to the element and
      performs the operations needed to introduce the seected time integration scheme.
      this function calculates at the same time the contribution to the LHS and to the RHS
      of the system
  */

  /**
   * This function is designed to be called in the builder and solver to introduce
   * @param pCurrentElement: The element to compute
   * @param rLHS_Contribution: The LHS matrix contribution
   * @param rRHS_Contribution: The RHS vector contribution
   * @param rEquationId: The ID's of the element degrees of freedom
   * @param rCurrentProcessInfo: The current process info instance
   */
  virtual void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                            LocalSystemMatrixType& rLHS_Contribution,
                                            LocalSystemVectorType& rRHS_Contribution,
                                            Element::EquationIdVectorType& rEquationId,
                                            ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentElement->CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);
  }

  /**
   * This function is designed to calculate just the RHS contribution
   * @param pCurrentElement: The element to compute
   * @param rRHS_Contribution: The RHS vector contribution
   * @param rEquationId: The ID's of the element degrees of freedom
   * @param rCurrentProcessInfo: The current process info instance
   */
  virtual void Calculate_RHS_Contribution(Element::Pointer pCurrentElement,
                                          LocalSystemVectorType& rRHS_Contribution,
                                          Element::EquationIdVectorType& rEquationId,
                                          ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentElement->CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
    pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);
  }
  /**
   * This function is designed to calculate just the LHS contribution
   * @param pCurrentElement: The element to compute
   * @param rLHS_Contribution: The LHS matrix contribution
   * @param rEquationId: The ID's of the element degrees of freedom
   * @param rCurrentProcessInfo: The current process info instance
   */
  virtual void Calculate_LHS_Contribution(Element::Pointer pCurrentElement,
                                          LocalSystemMatrixType& rLHS_Contribution,
                                          Element::EquationIdVectorType& rEquationId,
                                          ProcessInfo& rCurrentProcessInfo)
  {
    std::cout<< " it is C_LHS "<<std::endl;
    pCurrentElement->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
    pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);
  }

  virtual void EquationId(Element::Pointer pCurrentElement,
                          Element::EquationIdVectorType& rEquationId,
                          ProcessInfo& rCurrentProcessInfo)
  {
    (pCurrentElement)->EquationIdVector(rEquationId, rCurrentProcessInfo);
  }

  /**
   * Functions totally analogous to the precedent but applied to the "condition" objects
   * @param pCurrentCondition: The condition to compute
   * @param rLHS_Contribution: The LHS matrix contribution
   * @param rRHS_Contribution: The RHS vector contribution
   * @param rEquationId: The ID's of the element degrees of freedom
   * @param rCurrentProcessInfo: The current process info instance
   */
  virtual void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
                                                      LocalSystemMatrixType& rLHS_Contribution,
                                                      LocalSystemVectorType& rRHS_Contribution,
                                                      Element::EquationIdVectorType& rEquationId,
                                                      ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentCondition->CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    pCurrentCondition->EquationIdVector(rEquationId, rCurrentProcessInfo);
  }
  /**
   * Functions totally analogous to the precedent but applied to the "condition" objects
   * @param rCurrentCondition: The condition to compute
   * @param rLHS_Contribution: The LHS matrix contribution
   * @param rRHS_Contribution: The RHS vector contribution
   * @param rEquationId: The ID's of the element degrees of freedom
   * @param rCurrentProcessInfo: The current process info instance
   */
  virtual void Condition_Calculate_RHS_Contribution(Condition::Pointer pCurrentCondition,
                                                    LocalSystemVectorType& rRHS_Contribution,
                                                    Element::EquationIdVectorType& rEquationId,
                                                    ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentCondition->CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
    pCurrentCondition->EquationIdVector(rEquationId, rCurrentProcessInfo);
  }

  /**
   * Functions that calculates the RHS of a "condition" object
   * @param rCurrentCondition: The condition to compute
   * @param rRHS_Contribution: The RHS vector contribution
   * @param rEquationId: The ID's of the condition degrees of freedom
   * @param rCurrentProcessInfo: The current process info instance
   */
  virtual void Condition_Calculate_LHS_Contribution(Condition::Pointer pCurrentCondition,
                                                    LocalSystemMatrixType& rLHS_Contribution,
                                                    Element::EquationIdVectorType& rEquationId,
                                                    ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentCondition->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
    pCurrentCondition->EquationIdVector(rEquationId, rCurrentProcessInfo);
  }

  virtual void Condition_EquationId(Condition::Pointer pCurrentCondition,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo)
  {
    (pCurrentCondition)->EquationIdVector(rEquationId, rCurrentProcessInfo);
  }

  /**
   * Function that returns the list of Degrees of freedom to be assembled in the system for a Given Element
   * @param rCurrentElement: The element to compute
   * @param rElementalDofsList: The element dofs list
   * @param rCurrentProcessInfo: The current process info instance
   */
  virtual void GetElementalDofList(Element::Pointer pCurrentElement,
                                   Element::DofsVectorType& rElementalDofList,
                                   ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentElement->GetDofList(rElementalDofList, rCurrentProcessInfo);
  }

  /**
   * Function that returns the list of Degrees of freedom to be assembled in the system for a Given Element
   * @param rCurrentCondition: The condition to compute
   * @param rConditionDofsList: The condition dofs list
   * @param rCurrentProcessInfo: The current process info instance
   */
  virtual void GetConditionDofList(Condition::Pointer pCurrentCondition,
                                   Element::DofsVectorType& rConditionDofList,
                                   ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentCondition->GetDofList(rConditionDofList, rCurrentProcessInfo);
  }

  ///@}
  ///@name Access
  ///@{


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
   * @brief Set process to execute after move_mesh
   */
  void SetProcess( ProcessPointerType pProcess )
  {
    mProcesses.push_back(pProcess); //NOTE: order set = order of execution
  }

  /**
   * @brief Set list of processes to execute after move_mesh
   */
  void SetProcessVector( ProcessPointerVectorType& rProcessVector )
  {
    mProcesses = rProcessVector;
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

  // Flags to set options
  Flags mOptions;

  // Time integration methods
  IntegrationMethodsVectorType mTimeVectorIntegrationMethods;
  IntegrationMethodsScalarType mTimeScalarIntegrationMethods;

  // Processes called after move mesh is called
  ProcessPointerVectorType mProcesses;

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{

  /**
   * @brief Initialize the elements.
   * @details This is intended to be called just once when the strategy is initialized
   */
  virtual void InitializeElements(ModelPart& rModelPart)
  {
    KRATOS_TRY

    if( this->IsNot(LocalFlagType::ELEMENTS_INITIALIZED) ){

#pragma omp parallel for
      for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++)
      {
        auto itElem = rModelPart.ElementsBegin() + i;
        itElem->Initialize();
      }

    this->Set(LocalFlagType::ELEMENTS_INITIALIZED, true);

    }

    KRATOS_CATCH("")
  }

  /**
   * @brief Initialize the conditions.
   * @details This is intended to be called just once when the strategy is initialized
   */
  virtual void InitializeConditions(ModelPart& rModelPart)
  {
    KRATOS_TRY

    if( this->IsNot(LocalFlagType::ELEMENTS_INITIALIZED) )
      KRATOS_ERROR << "Before initilizing Conditions, initialize Elements FIRST" << std::endl;

    if( this->IsNot(LocalFlagType::CONDITIONS_INITIALIZED) ){

#pragma omp parallel for
      for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++)
      {
        auto itCond = rModelPart.ConditionsBegin() + i;
        itCond->Initialize();
      }

      this->Set(LocalFlagType::CONDITIONS_INITIALIZED, true);
    }

    KRATOS_CATCH("")
  }

  /**
   * @brief Initialize the conditions.
   * @details this function must be called only once per iteration
   */
  virtual void InitializeNonLinearIteration(Condition::Pointer rCurrentCondition,
                                            ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    rCurrentCondition->InitializeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH("")
  }

  /**
   * @brief Initialize the elements.
   * @details This is intended to be called every iteration
   */
  virtual void InitializeNonLinearIteration(Element::Pointer rCurrentElement,
                                            ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    rCurrentElement->InitializeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH("")
  }


  virtual void IntegrationMethodUpdate(NodeType& rNode)
  {
    for(typename IntegrationMethodsVectorType::iterator it=mTimeVectorIntegrationMethods.begin();
        it!=mTimeVectorIntegrationMethods.end(); ++it)
      (*it)->Update(rNode);
    for(typename IntegrationMethodsScalarType::iterator it=mTimeScalarIntegrationMethods.begin();
        it!=mTimeScalarIntegrationMethods.end(); ++it)
      (*it)->Update(rNode);
  }

  virtual void IntegrationMethodPredict(NodeType& rNode)
  {
    for(typename IntegrationMethodsVectorType::iterator it=mTimeVectorIntegrationMethods.begin();
        it!=mTimeVectorIntegrationMethods.end(); ++it)
      (*it)->Predict(rNode);
    for(typename IntegrationMethodsScalarType::iterator it=mTimeScalarIntegrationMethods.begin();
        it!=mTimeScalarIntegrationMethods.end(); ++it)
      (*it)->Predict(rNode);
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
  ///@name Serialization
  ///@{

  ///@}
  ///@name Private Inquiry
  ///@{

  ///@}
  ///@name Un accessible methods
  ///@{

  ///@}
}; // Class SolutionScheme
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SOLUTION_SCHEME_H_INCLUDED defined
