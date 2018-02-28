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

#if !defined(KRATOS_EXPLICIT_STRATEGY)
#define KRATOS_EXPLICIT_STRATEGY

/* System includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos {

/**
 * @class ExplicitStrategy
 *
 * @brief This strategy is used for the explicit time integration
 *
 * @author Klauss B Sautter (based on the work of JMCarbonel)
 */

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ExplicitStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
  /** Counted pointer of ClassName */

  KRATOS_CLASS_POINTER_DEFINITION(ExplicitStrategy);

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

  ExplicitStrategy(ModelPart &rModelPart, typename TSchemeType::Pointer pScheme,
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

    // saving the linear solver
    // mpLinearSolver = pNewLinearSolver; //Not used in explicit strategies

    // set flags to start correcty the calculations
    this->mSolutionStepIsInitialized = false;
    this->mInitializeWasPerformed = false;

    // set EchoLevel to the deffault value (only time is displayed)
    SetEchoLevel(1);

    // set RebuildLevel to the deffault value
    BaseType::SetRebuildLevel(0);

    KRATOS_CATCH("")
  }

  /** Destructor.
   */
  virtual ~ExplicitStrategy() {}

  /** Destructor.
   */

  // Set and Get Scheme ... containing Update and other

  void SetScheme(typename TSchemeType::Pointer pScheme) {
    this->mpScheme = pScheme;
  };

  typename TSchemeType::Pointer GetScheme() { return this->mpScheme; };

  // Ser and Get Flags

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

  // level of echo for the solving strategy
  // 0 -> mute... no echo at all
  // 1 -> printing time and basic informations
  // 2 -> printing linear solver data
  // 3 -> Print of debug informations:
  //    Echo of stiffness matrix, Dx, b...

  void SetEchoLevel(int Level) { this->BaseType::mEchoLevel = Level; }

  //**********************************************************************
  //**********************************************************************

  void Initialize() {
    KRATOS_TRY
    // pointers needed in the solution
    typename TSchemeType::Pointer pScheme = GetScheme();
    ModelPart &r_model_part = BaseType::GetModelPart();

    TSystemMatrixType matrix_a_dummy = TSystemMatrixType();

    // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
    if (pScheme->SchemeIsInitialized() == false)
      pScheme->Initialize(r_model_part);

    // Initialize The Elements - OPERATIONS TO BE DONE ONCE
    if (pScheme->ElementsAreInitialized() == false)
      pScheme->InitializeElements(r_model_part);

    // Initialize The Conditions- OPERATIONS TO BE DONE ONCE
    if (pScheme->ConditionsAreInitialized() == false)
      pScheme->InitializeConditions(r_model_part);

    // Set Nodal Mass to zero
    NodesArrayType &r_nodes = r_model_part.Nodes();
    ElementsArrayType &r_elements = r_model_part.Elements();
    ProcessInfo &r_current_process_info = r_model_part.GetProcessInfo();

    auto it_node = r_model_part.NodesBegin();
#pragma omp parallel for firstprivate(it_node)
    for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
      (it_node + i)->SetValue(NODAL_MASS, 0.00);
    }

    if (r_model_part.NodesBegin()->HasDofFor(ROTATION_Z)) {
#pragma omp parallel for
      for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
        array_1d<double, 3> &r_nodal_inertia =
            (it_node + i)->GetValue(NODAL_INERTIA);
        noalias(r_nodal_inertia) = ZeroVector(3);
      }
    }

    auto it_elem = r_model_part.ElementsBegin();
#pragma omp parallel for firstprivate(it_elem)
    for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
      // Getting nodal mass and inertia from element
      Vector dummy_vector;
      // this function needs to be implemented in the respective
      // element to provide inertias and nodal masses
      (it_elem + i)
          ->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR,
                                    NODAL_INERTIA, r_current_process_info);
    }

    mInitializeWasPerformed = true;

    // std::cout<<" Rebuild Level "<<BaseType::mRebuildLevel<<std::endl;

    KRATOS_CATCH("")
  }

  //**********************************************************************
  //**********************************************************************

  void InitializeSolutionStep() {
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
#pragma omp parallel for firstprivate(it_elem)
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

    mSolutionStepIsInitialized = true;

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
  /*
                    SOLUTION OF THE PROBLEM OF INTEREST
   */
  //**********************************************************************

  double Solve() override {
    KRATOS_TRY
    DofsArrayType dof_set_dummy;
    TSystemMatrixType mA = TSystemMatrixType();
    TSystemVectorType mDx = TSystemVectorType();
    TSystemVectorType mb = TSystemVectorType();

    // pointers needed in the solution
    typename TSchemeType::Pointer pScheme = GetScheme();
    ModelPart &r_model_part = BaseType::GetModelPart();

    // OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
    // if the operations needed were already performed this does nothing
    if (mInitializeWasPerformed == false)
      Initialize();

    // prints informations about the current time
    if (this->GetEchoLevel() == 2 &&
        r_model_part.GetCommunicator().MyPID() == 0) {
      std::cout << " " << std::endl;
      std::cout << "CurrentTime = " << r_model_part.GetProcessInfo()[TIME]
                << std::endl;
    }

    // initialize solution step
    if (mSolutionStepIsInitialized == false)
      InitializeSolutionStep();

    this->CalculateAndAddRHS(pScheme, r_model_part);

    pScheme->Update(r_model_part, dof_set_dummy, mA, mDx,
                    mb); // Explicitly integrates the equation of motion.
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

    // reset flags for next step
    mSolutionStepIsInitialized = false;
    return 0.00;

    KRATOS_CATCH("")
  }

  //**********************************************************************
  //**********************************************************************

  void Clear() {
    KRATOS_TRY
    std::cout << "Explicit strategy Clear function used" << std::endl;

    GetScheme()->Clear();

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

  bool mSolutionStepIsInitialized;

  bool mInitializeWasPerformed;

  bool mComputeTime;

  /*@} */
  /**@name Private Operators*/
  /*@{ */

  //**********************************************************************
  //**********************************************************************

  void CalculateReactions() {}

  //**********************************************************************
  //**********************************************************************

  /**
   * function to perform expensive checks.
   * It is designed to be called ONCE to verify that the input is correct.
   */

  int Check() {
    KRATOS_TRY

    BaseType::Check();

    GetScheme()->Check(BaseType::GetModelPart());

    return 0;

    KRATOS_CATCH("")
  }

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
  ExplicitStrategy(const ExplicitStrategy &Other){};

  /*@} */

}; /* Class ExplicitStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_STRATEGY  defined */
