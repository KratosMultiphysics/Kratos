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

/** @brief Solver local flags class definition
 *  @details This is the base class for solver local flags
 */
class SolverLocalFlags
{
 public:
  /// Flags for the solution control:
  KRATOS_DEFINE_LOCAL_FLAG( INITIALIZED );
  KRATOS_DEFINE_LOCAL_FLAG( CONVERGED );
  KRATOS_DEFINE_LOCAL_FLAG( DOFS_INITIALIZED );
  KRATOS_DEFINE_LOCAL_FLAG( ELEMENTS_INITIALIZED );
  KRATOS_DEFINE_LOCAL_FLAG( CONDITIONS_INITIALIZED );
  
  /// Flags for the solution options:
  KRATOS_DEFINE_LOCAL_FLAG( REFORM_DOFS );
  KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_REACTIONS );
  KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_SYSTEM_MATRIX );
  KRATOS_DEFINE_LOCAL_FLAG( RAYLEIGH_DAMPING );
  KRATOS_DEFINE_LOCAL_FLAG( IMPLEX );
  
};


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

  typedef ModelPart::NodeType                                             NodeType;
  typedef array_1d<double, 3>                                           VectorType;
  typedef Variable<VectorType>                                        VariableType; 
  typedef TimeIntegrationMethod<VariableType,VectorType>           IntegrationType;
  typedef typename IntegrationType::Pointer                 IntegrationPointerType;

  /// Pointer definition of SolutionScheme
  KRATOS_CLASS_POINTER_DEFINITION(SolutionScheme);
  
  ///@}
  ///@name Life Cycle
  ///@{

  /// Default Constructor.
  SolutionScheme() : Flags() {}

  /// Constructor.
  SolutionScheme(Flags& rOptions) : Flags(), mOptions(rOptions) {}
  
  /// Constructor.
  SolutionScheme(IntegrationPointerType pIntegrationMethod, Flags& rOptions) : Flags(), mOptions(rOptions), mpIntegrationMethod(pIntegrationMethod) {}
  
  /// Copy contructor.
  SolutionScheme(SolutionScheme& rOther) : mOptions(rOther.mOptions), mpIntegrationMethod(rOther.mpIntegrationMethod) {}

  /// Clone.
  virtual SolutionSchemePointerType Clone()
  {
    return SolutionSchemePointerType( new SolutionScheme(*this) );
  }

  /// Destructor.
  virtual ~SolutionScheme() {}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{


  /**
   * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
   * @details This is intended to be called just once when the strategy is initialized
   */
  virtual void Initialize(ModelPart& rModelPart)
  {
    KRATOS_TRY

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
   * @brief This function is designed to move the mesh
   * @note Be careful it just consider displacements, derive this method to adapt to your own strategies (ALE, FSI, etc...)
   */
  virtual void MoveMesh(ModelPart& rModelPart)
  {
    KRATOS_TRY

    if (rModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT_X) == false)
    {
      KRATOS_ERROR << "It is impossible to move the mesh since the DISPLACEMENT variable is not in the Model Part. Add DISPLACEMENT to the list of variables" << std::endl;
    }

    const int nnodes = rModelPart.NumberOfNodes();
    ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();

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

    
    KRATOS_CATCH("")

    return 0;
  }

  /** this function is designed to be called in the builder and solver to introduce
      the selected time integration scheme. It "asks" the matrix needed to the element and
      performs the operations needed to introduce the seected time integration scheme.

      this function calculates at the same time the contribution to the LHS and to the RHS
      of the system
  */
  virtual void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                            LocalSystemMatrixType& rLHS_Contribution,
                                            LocalSystemVectorType& rRHS_Contribution,
                                            Element::EquationIdVectorType& rEquationId,
                                            ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentElement->CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
  }

  virtual void Calculate_RHS_Contribution(Element::Pointer pCurrentElement,
                                          LocalSystemVectorType& rRHS_Contribution,
                                          Element::EquationIdVectorType& rEquationId,
                                          ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentElement->CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
  }

  virtual void Calculate_LHS_Contribution(Element::Pointer pCurrentElement,
                                          LocalSystemMatrixType& rLHS_Contribution,
                                          Element::EquationIdVectorType& rEquationId,
                                          ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentElement->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
  }

  virtual void EquationId(Element::Pointer pCurrentElement,
                          Element::EquationIdVectorType& rEquationId,
                          ProcessInfo& rCurrentProcessInfo)
  {
    (pCurrentElement)->EquationIdVector(rEquationId, rCurrentProcessInfo);
  }

  /** functions totally analogous to the precedent but applied to
      the "condition" objects
  */
  virtual void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
                                                      LocalSystemMatrixType& rLHS_Contribution,
                                                      LocalSystemVectorType& rRHS_Contribution,
                                                      Element::EquationIdVectorType& rEquationId,
                                                      ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentCondition->CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
  }

  virtual void Condition_Calculate_RHS_Contribution(Condition::Pointer pCurrentCondition,
                                                    LocalSystemVectorType& rRHS_Contribution,
                                                    Element::EquationIdVectorType& rEquationId,
                                                    ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentCondition->CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
  }

  virtual void Condition_Calculate_LHS_Contribution(Condition::Pointer pCurrentCondition,
                                                    LocalSystemMatrixType& rLHS_Contribution,
                                                    Element::EquationIdVectorType& rEquationId,
                                                    ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentCondition->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
  }

  virtual void Condition_EquationId(Condition::Pointer pCurrentCondition,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo)
  {
    (pCurrentCondition)->EquationIdVector(rEquationId, rCurrentProcessInfo);
  }

  /** Function that returns the list of Degrees of freedom to be
      assembled in the system for a Given Element
  */
  virtual void GetElementalDofList(Element::Pointer pCurrentElement,
                                   Element::DofsVectorType& rElementalDofList,
                                   ProcessInfo& rCurrentProcessInfo)
  {
    pCurrentElement->GetDofList(rElementalDofList, rCurrentProcessInfo);
  }

  /** Function that returns the list of Degrees of freedom to be
      assembled in the system for a Given Element
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

  // Time integration method
  IntegrationPointerType    mpIntegrationMethod;
  
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

/**
 * Flags for the Strategy control
 */
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, INITIALIZED,               0 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, CONVERGED,                 1 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, DOFS_INITIALIZED,          2 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, ELEMENTS_INITIALIZED,      3 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, CONDITIONS_INITIALIZED,    4 );

/**
 * Flags for the Strategy options
 */
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, REFORM_DOFS,               0 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, COMPUTE_REACTIONS,         1 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, CONSTANT_SYSTEM_MATRIX,    2 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, RAYLEIGH_DAMPING,          3 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, IMPLEX,                    4 );
///@}
///@name Input and output
///@{

  
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SOLUTION_SCHEME_H_INCLUDED defined


