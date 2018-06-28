//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Klaus B Sautter (based on the work of JMCarbonel)
//
//

#if !defined(KRATOS_MECHANICAL_EXPLICIT_STRATEGY)
#define KRATOS_MECHANICAL_EXPLICIT_STRATEGY

/* System includes */

/* Project includes */
#include "solving_strategies/strategies/solving_strategy.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/variable_utils.h"

namespace Kratos {

/**
 * @class MechanicalExplicitStrategy
 *
 * @brief This strategy is used for the explicit time integration
 *
 * @author Klauss B Sautter (based on the work of JMCarbonel)
 */

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class MechanicalExplicitStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
  /** Counted pointer of ClassName */

  KRATOS_CLASS_POINTER_DEFINITION(MechanicalExplicitStrategy);

  typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

  typedef typename BaseType::TSchemeType TSchemeType;
  typedef typename BaseType::DofsArrayType DofsArrayType;
  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
  typedef typename BaseType::TSystemVectorType TSystemVectorType;
  typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
  typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
  typedef typename BaseType::NodesArrayType NodesArrayType;
  typedef typename BaseType::ElementsArrayType ElementsArrayType;
  typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

  /** Constructors.
   */

  MechanicalExplicitStrategy(ModelPart &rModelPart, typename TSchemeType::Pointer pScheme,
                   bool CalculateReactions = false,
                   bool ReformDofSetAtEachStep = false,
                   bool MoveMeshFlag = true)
      : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
            rModelPart, MoveMeshFlag) {
    KRATOS_TRY

    // set flags to default values
    this->mCalculateReactionsFlag = CalculateReactions;
    this->mReformDofSetAtEachStep = ReformDofSetAtEachStep;

    // saving the scheme
    this->mpScheme = pScheme;

    // set EchoLevel to the default value (only time is displayed)
    BaseType::SetEchoLevel(1);

    // set RebuildLevel to the default value
    BaseType::SetRebuildLevel(0);

    KRATOS_CATCH("")
  }

  /** Destructor.
   */
  virtual ~MechanicalExplicitStrategy() {}

  /** Destructor.
   */

  // Set and Get Scheme ... containing Update and other

  void SetScheme(typename TSchemeType::Pointer pScheme) {
    this->mpScheme = pScheme;
  };

  typename TSchemeType::Pointer GetScheme() { return this->mpScheme; };

  // Set and Get Flags

  void SetInitializePerformedFlag(bool InitializePerformedFlag = true) {
    this->mInitializeWasPerformed = InitializePerformedFlag;
  }

  bool GetInitializePerformedFlag() { return this->mInitializeWasPerformed; }

  void SetCalculateReactionsFlag(bool CalculateReactionsFlag) {
    this->mCalculateReactionsFlag = CalculateReactionsFlag;
  }

  bool GetCalculateReactionsFlag() { return this->mCalculateReactionsFlag; }

  void SetReformDofSetAtEachStepFlag(bool flag) {
    this->mReformDofSetAtEachStep = flag;
  }

  bool GetReformDofSetAtEachStepFlag() { return this->mReformDofSetAtEachStep; }

  //**********************************************************************
  //**********************************************************************

  void Initialize() override {
    KRATOS_TRY
    if (!this->mInitializeWasPerformed){
      // pointers needed in the solution
      typename TSchemeType::Pointer pScheme = GetScheme();
      ModelPart &r_model_part = BaseType::GetModelPart();

      TSystemMatrixType matrix_a_dummy = TSystemMatrixType();

      // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
      if (!pScheme->SchemeIsInitialized())pScheme->Initialize(r_model_part);

      // Initialize The Elements - OPERATIONS TO BE DONE ONCE
      if (!pScheme->ElementsAreInitialized())pScheme->InitializeElements(r_model_part);

      // Initialize The Conditions- OPERATIONS TO BE DONE ONCE
      if (!pScheme->ConditionsAreInitialized())pScheme->InitializeConditions(r_model_part);

      // Set Nodal Mass to zero
      NodesArrayType &r_nodes = r_model_part.Nodes();
      ElementsArrayType &r_elements = r_model_part.Elements();
      ProcessInfo &r_current_process_info = r_model_part.GetProcessInfo();

      VariableUtils().SetNonHistoricalVariable(NODAL_MASS, 0.0, r_nodes);
      const array_1d<double, 3> zero_array(3, 0.0);
      if (r_model_part.NodesBegin()->HasDofFor(ROTATION_Z))
        VariableUtils().SetNonHistoricalVariable(NODAL_INERTIA, zero_array, r_nodes);

      auto it_elem = r_model_part.ElementsBegin();
      // #pragma omp parallel for firstprivate(it_elem)
      for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
        // Getting nodal mass and inertia from element
        Vector dummy_vector;
        // this function needs to be implemented in the respective
        // element to provide inertias and nodal masses
        (it_elem + i)
            ->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR,
                                      NODAL_INERTIA, r_current_process_info);
      }

      this->mInitializeWasPerformed = true;
    }
    KRATOS_CATCH("")
  }

  //**********************************************************************
  //**********************************************************************

  void InitializeSolutionStep() override {
    KRATOS_TRY

    typename TSchemeType::Pointer pScheme = GetScheme();
    ModelPart &r_model_part = BaseType::GetModelPart();

    TSystemMatrixType matrix_a_dummy = TSystemMatrixType();
    TSystemVectorType mDx = TSystemVectorType();
    TSystemVectorType mb = TSystemVectorType();

    // initial operations ... things that are constant over the Solution Step
    pScheme->InitializeSolutionStep(r_model_part, matrix_a_dummy, mDx, mb);

    ProcessInfo &r_current_process_info = r_model_part.GetProcessInfo();
    ElementsArrayType &r_elements = r_model_part.Elements();

    if (BaseType::mRebuildLevel > 0) {
    auto it_elem = r_model_part.ElementsBegin();
// #pragma omp parallel for firstprivate(it_elem)
    for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
        // Getting nodal mass and inertia from element
        Vector dummy_vector;
        // this function needs to be implemented in the respective
        // element to provide inertias and nodal masses
        (it_elem + i)
            ->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR,
                                    NODAL_INERTIA, r_current_process_info);
    }
    }

    KRATOS_CATCH("")
  }

  void CalculateAndAddRHS(typename TSchemeType::Pointer pScheme,
                          ModelPart &rModelPart) {
    KRATOS_TRY

    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    ConditionsArrayType &pConditions = rModelPart.Conditions();
    ElementsArrayType &pElements = rModelPart.Elements();

    typename ConditionsArrayType::ptr_iterator it_begin_condition =
        pConditions.ptr_begin();
#pragma omp parallel for firstprivate(it_begin_condition)
    for (int i = 0; i < static_cast<int>(pConditions.size()); ++i) {
      LocalSystemVectorType RHS_Condition_Contribution =
          LocalSystemVectorType(0);
      Element::EquationIdVectorType equation_id_vector_dummy; // Dummy

      pScheme->Condition_Calculate_RHS_Contribution(
          *(it_begin_condition + i), RHS_Condition_Contribution,
          equation_id_vector_dummy, rCurrentProcessInfo);
    }

    typename ElementsArrayType::ptr_iterator it_begin_element =
        pElements.ptr_begin();
#pragma omp parallel for firstprivate(it_begin_element)
    for (int i = 0; i < static_cast<int>(pElements.size()); ++i) {
      LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
      Element::EquationIdVectorType equation_id_vector_dummy; // Dummy

      pScheme->Calculate_RHS_Contribution(
          *(it_begin_element + i), RHS_Contribution, equation_id_vector_dummy,
          rCurrentProcessInfo);
    }
    KRATOS_CATCH("")
  }
  //**********************************************************************
  //**********************************************************************
  bool SolveSolutionStep() override
  {
    typename TSchemeType::Pointer pScheme = GetScheme();
    ModelPart &r_model_part = BaseType::GetModelPart();
    DofsArrayType dof_set_dummy;
    TSystemMatrixType mA = TSystemMatrixType();
    TSystemVectorType mDx = TSystemVectorType();
    TSystemVectorType mb = TSystemVectorType();

    pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

    this->CalculateAndAddRHS(pScheme, r_model_part);

    pScheme->Update(r_model_part, dof_set_dummy, mA, mDx,
                    mb); // Explicitly integrates the equation of motion.
    return true;
  }
  //**********************************************************************
  //**********************************************************************
  void FinalizeSolutionStep() override
  {
    typename TSchemeType::Pointer pScheme = GetScheme();
    ModelPart &r_model_part = BaseType::GetModelPart();
    TSystemMatrixType mA = TSystemMatrixType();
    TSystemVectorType mDx = TSystemVectorType();
    TSystemVectorType mb = TSystemVectorType();
    // Finalisation of the solution step,
    // operations to be done after achieving convergence, for example the
    // Final Residual Vector (mb) has to be saved in there
    // to avoid error accumulation
    pScheme->FinalizeSolutionStep(r_model_part, mA, mDx, mb);

    // move the mesh if needed
    if (BaseType::MoveMeshFlag() == true)
      BaseType::MoveMesh();

    // Cleaning memory after the solution
    pScheme->Clean();
  }

  //**********************************************************************
  //**********************************************************************
  void Clear() override {
    KRATOS_TRY

    KRATOS_INFO("MechanicalExplicitStrategy")
      << "Clear function used" << std::endl;

    GetScheme()->Clear();
    mInitializeWasPerformed = false;

    KRATOS_CATCH("")
  }
  //**********************************************************************
  //**********************************************************************

  /**
   * function to perform expensive checks.
   * It is designed to be called ONCE to verify that the input is correct.
   */

  int Check() override {
    KRATOS_TRY

    BaseType::Check();

    GetScheme()->Check(BaseType::GetModelPart());

    return 0;

    KRATOS_CATCH("")
  }

  /*@} */
  /**@name Operators
   */
  /*@{ */

  /*@} */
  /**@name Operations */
  /*@{ */

  /*@} */
  /**@name Access */

  /*@{ */

  /*@} */
  /**@name Inquiry */
  /*@{ */

  /*@} */
  /**@name Friends */
  /*@{ */

  /*@} */

private:
  /**@name Protected static Member Variables */
  /*@{ */

  /*@} */
  /**@name Protected member Variables */
  /*@{ */

  /*@} */
  /**@name Protected Operators*/
  /*@{ */

  /*@} */
  /**@name Protected Operations*/
  /*@{ */

  /*@} */
  /**@name Protected  Access */
  /*@{ */

  /*@} */
  /**@name Protected Inquiry */
  /*@{ */

  /*@} */
  /**@name Protected LifeCycle */
  /*@{ */

  /*@} */

protected:
  /**@name Static Member Variables */
  /*@{ */

  /*@} */
  /**@name Member Variables */
  /*@{ */

  typename TSchemeType::Pointer mpScheme;

  typename TLinearSolver::Pointer mpLinearSolver;

  TSystemVectorPointerType mpDx;
  TSystemVectorPointerType mpb;
  TSystemMatrixPointerType mpA;

  /**
  Flag telling if it is needed to reform the DofSet at each
  solution step or if it is possible to form it just once
  - true  => reforme at each time step
  - false => form just one (more efficient)

  Default = false
   */
  bool mReformDofSetAtEachStep;

  /**
  Flag telling if it is needed or not to compute the reactions

  default = true
   */
  bool mCalculateReactionsFlag;

  bool mInitializeWasPerformed = false;

  /*@} */
  /**@name Private Operators*/
  /*@{ */


  //***************************************************************************
  //***************************************************************************

  /*@} */
  /**@name Private Operations*/
  /*@{ */

  /*@} */
  /**@name Private  Access */
  /*@{ */

  /*@} */
  /**@name Private Inquiry */
  /*@{ */

  /*@} */
  /**@name Un accessible methods */
  /*@{ */

  /** Copy constructor.
   */
  MechanicalExplicitStrategy(const MechanicalExplicitStrategy &Other){};

  /*@} */

}; /* Class MechanicalExplicitStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_STRATEGY  defined */
