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
//  Main authors:    Klaus B Sautter (based on the work of MSantasusana)
//
//

#if !defined(KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME_HPP_INCLUDED)
#define KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME_HPP_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/scheme.h"
#include "utilities/variable_utils.h"

namespace Kratos {

/**@name Kratos Globals */
/*@{ */

/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */

/**@name  Enum's */
/*@{ */

/*@} */
/**@name  Functions */
/*@{ */

/*@} */
/**@name Kratos Classes */
/*@{ */

template <class TSparseSpace,
          class TDenseSpace //= DenseSpace<double>
          >
class ExplicitCentralDifferencesScheme
    : public Scheme<TSparseSpace, TDenseSpace> {

public:
  /**@name Type Definitions */
  /*@{ */
  KRATOS_CLASS_POINTER_DEFINITION(ExplicitCentralDifferencesScheme);

  typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

  typedef typename BaseType::DofsArrayType DofsArrayType;
  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
  typedef typename BaseType::TSystemVectorType TSystemVectorType;
  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
  typedef ModelPart::ElementsContainerType ElementsArrayType;
  typedef ModelPart::NodesContainerType NodesArrayType;
  typedef std::size_t SizeType;
  typedef typename ModelPart::NodeIterator NodeIterator;

  ExplicitCentralDifferencesScheme(const double MaximumDeltaTime,
                                   const double DeltaTimeFraction,
                                   const double DeltaTimePredictionLevel)
      : Scheme<TSparseSpace, TDenseSpace>() {
    mDeltaTime.PredictionLevel = DeltaTimePredictionLevel;

    mDeltaTime.Maximum = MaximumDeltaTime;

    mDeltaTime.Fraction = DeltaTimeFraction;
  }

  /** Destructor.
   */
  virtual ~ExplicitCentralDifferencesScheme() {}

  /*@} */
  /**@name Operators
   */
  /*@{ */

  virtual int Check(ModelPart &rModelPart) override {
    KRATOS_TRY

    BaseType::Check(rModelPart);

    KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2)
        << "Insufficient buffer size for Central Difference Scheme. It has to "
           "be > 2"
        << std::endl;

    return 0;
    KRATOS_CATCH("");
  }

  virtual void Initialize(ModelPart &rModelPart) override {
    KRATOS_TRY

    if ((mDeltaTime.PredictionLevel > 0) && (!BaseType::SchemeIsInitialized())) {
      CalculateDeltaTime(rModelPart);
    }

    ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();

    // Preparing the time values for the first step (where time = initial_time +
    // dt)
    mTime.Current =
        r_current_process_info[TIME] + r_current_process_info[DELTA_TIME];
    mTime.Delta = r_current_process_info[DELTA_TIME];
    mTime.Middle = mTime.Current - 0.5 * mTime.Delta;
    mTime.Previous = mTime.Current - mTime.Delta;
    mTime.PreviousMiddle = mTime.Current - 1.5 * mTime.Delta;

    if (!BaseType::SchemeIsInitialized())InitializeExplicitScheme(rModelPart);
    else SchemeCustomInitialization(rModelPart);

    BaseType::SetSchemeIsInitialized();
    KRATOS_CATCH("")
  }

  //***************************************************************************

  virtual void InitializeSolutionStep(ModelPart &rModelPart,
                                      TSystemMatrixType &rA,
                                      TSystemVectorType &rDx,
                                      TSystemVectorType &rb) override {
    KRATOS_TRY
    BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);
    if (mDeltaTime.PredictionLevel > 1)CalculateDeltaTime(rModelPart);
    InitializeResidual(rModelPart);
    KRATOS_CATCH("")
  }
  //***************************************************************************
  
  void InitializeNonLinIteration(
      ModelPart& rModelPart,
      TSystemMatrixType& A,
      TSystemVectorType& Dx,
      TSystemVectorType& b
      ) override
  {
      KRATOS_TRY;

      ProcessInfo& current_process_info = rModelPart.GetProcessInfo();
      
      #pragma omp parallel for
      for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i) {
          auto it_elem = rModelPart.ElementsBegin() + i;
          it_elem->InitializeNonLinearIteration(current_process_info);
      }
      
      
      #pragma omp parallel for
      for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
          auto it_elem = rModelPart.ConditionsBegin() + i;
          it_elem->InitializeNonLinearIteration(current_process_info);
      }     
      
      KRATOS_CATCH( "" );
  }


  //**************************************************************************

  void InitializeResidual(ModelPart &rModelPart) {
    KRATOS_TRY
    NodesArrayType &r_nodes = rModelPart.Nodes();
    VariableUtils().SetVectorVar(FORCE_RESIDUAL,ZeroVector(3),r_nodes);
    VariableUtils().SetVectorVar(MOMENT_RESIDUAL,ZeroVector(3),r_nodes);
    KRATOS_CATCH("")
  }

  //***************************************************************************

  void CalculateDeltaTime(ModelPart &rModelPart) {
    KRATOS_TRY

    ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
    ElementsArrayType &r_elements = rModelPart.Elements();

    const double safety_factor = 0.5;

    double delta_time = mDeltaTime.Maximum / safety_factor;

    double stable_delta_time = 1000.00;

    auto it_begin = rModelPart.ElementsBegin();

#pragma omp parallel for firstprivate(it_begin)
      for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
        bool check_has_all_variables = true;
        double E(0.00), nu(0.00), roh(0.00), alpha(0.00), beta(0.00);
        // get geometric and material properties
        if ((it_begin + i)->GetProperties().Has(RAYLEIGH_ALPHA)) {
          alpha = (it_begin + i)->GetProperties()[RAYLEIGH_ALPHA];
        }
        if ((it_begin + i)->GetProperties().Has(RAYLEIGH_BETA)) {
          beta = (it_begin + i)->GetProperties()[RAYLEIGH_BETA];
        }
        if ((it_begin + i)->GetProperties().Has(YOUNG_MODULUS)) {
          E = (it_begin + i)->GetProperties()[YOUNG_MODULUS];
        } else
          check_has_all_variables = false;
        if ((it_begin + i)->GetProperties().Has(POISSON_RATIO)) {
          nu = (it_begin + i)->GetProperties()[POISSON_RATIO];
        }
        if ((it_begin + i)->GetProperties().Has(DENSITY)) {
          roh = (it_begin + i)->GetProperties()[DENSITY];
        } else
          check_has_all_variables = false;

        if (check_has_all_variables) {
          const double length = (it_begin + i)->GetGeometry().Length();

          // compute courant criterion
          const double bulk_modulus = E / (3.0 * (1.0 - 2.0 * nu));
          const double wavespeed = std::sqrt(bulk_modulus / roh);
          const double w = 2.0 * wavespeed / length; // frequency

          const double psi = 0.5 * (alpha / w + beta * w); // critical ratio;
          stable_delta_time = (2.0 / w) * (std::sqrt(1.0 + psi * psi) - psi);

          if (stable_delta_time > 0.00) {
            #pragma omp critical
            if (stable_delta_time < delta_time) delta_time = stable_delta_time;
          }
        } else
          KRATOS_ERROR << "not enough parameters for prediction level "
                       << mDeltaTime.PredictionLevel << std::endl;
      }

    stable_delta_time = delta_time * safety_factor;

    if (stable_delta_time < mDeltaTime.Maximum) {

      r_current_process_info[DELTA_TIME] = stable_delta_time;
    }

    if (mDeltaTime.PredictionLevel > 1) {
      std::cout << "  [EXPLICIT PREDICTION LEVEL " << mDeltaTime.PredictionLevel
                << " ] : (computed stable time step = " << stable_delta_time
                << " s)" << std::endl;
      std::cout << "  Using  = " << r_current_process_info[DELTA_TIME]
                << " s as time step DELTA_TIME)" << std::endl;
    }
    KRATOS_CATCH("")
  }

  //***************************************************************************

  void InitializeExplicitScheme(ModelPart &rModelPart) {
    KRATOS_TRY

    NodesArrayType &r_nodes = rModelPart.Nodes();

    SizeType dim(3);

    auto i_begin = rModelPart.NodesBegin();
    //allocating memory outside of omp
    for (SizeType i=0;i<r_nodes.size();++i)
    {
      (i_begin + i)->SetValue(MIDDLE_VELOCITY,ZeroVector(3));
      (i_begin + i)->SetValue(MIDDLE_ANGULAR_VELOCITY,ZeroVector(3));
      (i_begin + i)->SetValue(NODAL_INERTIA,ZeroVector(3));
    }

    const bool has_dof_for_rot_z = i_begin->HasDofFor(ROTATION_Z);
#pragma omp parallel for firstprivate(i_begin)
    for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
      array_1d<double, 3> &r_middle_velocity =
          (i_begin + i)->GetValue(MIDDLE_VELOCITY);
      array_1d<double, 3> &r_current_velocity =
          (i_begin + i)->FastGetSolutionStepValue(VELOCITY);
      array_1d<double, 3> &r_current_residual =
          (i_begin + i)->FastGetSolutionStepValue(FORCE_RESIDUAL);
      // array_1d<double,3>& r_current_displacement  =
      // i_begin->FastGetSolutionStepValue(DISPLACEMENT);

      for (SizeType j = 0; j < dim; j++) {

        r_middle_velocity[j] = r_current_velocity[j];
        r_current_residual[j] = 0.0;
        // r_current_displacement[j] = 0.0; // this might be wrong for presribed
        // displacement
      }

      if (has_dof_for_rot_z) {
        array_1d<double, 3> &r_middle_angular_velocity =
            (i_begin + i)->GetValue(MIDDLE_ANGULAR_VELOCITY);
        array_1d<double, 3> &r_current_angular_velocity =
            (i_begin + i)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
        array_1d<double, 3> &r_current_residual_moment =
            (i_begin + i)->FastGetSolutionStepValue(MOMENT_RESIDUAL);
        // array_1d<double,3>& current_rotation              =
        // i->FastGetSolutionStepValue(ROTATION);

        for (SizeType j = 0; j < dim; j++) {

          r_middle_angular_velocity[j] = r_current_angular_velocity[j];
          r_current_residual_moment[j] = 0.0;
          // current_rotation[j] = 0.0; // this might be wrong for presribed
          // rotations
        }
      }
    }

    KRATOS_CATCH("")
  }

  /**
      Performing the update of the solution.
  */
  //***************************************************************************
  virtual void Update(ModelPart &rModelPart, DofsArrayType &rDofSet,
                      TSystemMatrixType &A, TSystemVectorType &Dx,
                      TSystemVectorType &b) override {
    KRATOS_TRY
    ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
    NodesArrayType &r_nodes = rModelPart.Nodes();
    // Step Update
    mTime.Current = r_current_process_info[TIME]; // the first step is time =
                                                  // initial_time ( 0.0) + delta
                                                  // time
    mTime.Delta = r_current_process_info[DELTA_TIME];

    mTime.Middle = 0.5 * (mTime.Previous + mTime.Current);

    auto i_begin = rModelPart.NodesBegin();
    const bool has_dof_for_rot_z = i_begin->HasDofFor(ROTATION_Z);
#pragma omp parallel for firstprivate(i_begin)
    for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
      // Current step information "N+1" (before step update).
      this->UpdateTranslationalDegreesOfFreedom(i_begin + i);
      if (has_dof_for_rot_z)
        this->UpdateRotationalDegreesOfFreedom(i_begin + i);
    } // for Node parallel

    mTime.Previous = mTime.Current;
    mTime.PreviousMiddle = mTime.Middle;

    KRATOS_CATCH("")
  }

  void UpdateTranslationalDegreesOfFreedom(NodeIterator CurrentNode) {
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double &nodal_mass = (CurrentNode)->GetValue(NODAL_MASS);
    array_1d<double, 3> &r_current_residual =
        (CurrentNode)->FastGetSolutionStepValue(FORCE_RESIDUAL);

    array_1d<double, 3> &r_current_velocity =
        (CurrentNode)->FastGetSolutionStepValue(VELOCITY);
    array_1d<double, 3> &r_current_displacement =
        (CurrentNode)->FastGetSolutionStepValue(DISPLACEMENT);
    array_1d<double, 3> &r_middle_velocity =
        (CurrentNode)->GetValue(MIDDLE_VELOCITY);

    array_1d<double, 3> &r_current_acceleration =
        (CurrentNode)->FastGetSolutionStepValue(ACCELERATION);

    // Solution of the explicit equation:
    if (nodal_mass > numerical_limit)
      r_current_acceleration = r_current_residual / nodal_mass;
    else
      r_current_acceleration = ZeroVector(3);

    SizeType DoF = 3;
    bool fix_displacements[3] = {false, false, false};

    fix_displacements[0] = ((CurrentNode)->pGetDof(DISPLACEMENT_X))->IsFixed();
    fix_displacements[1] = ((CurrentNode)->pGetDof(DISPLACEMENT_Y))->IsFixed();
    fix_displacements[2] = ((CurrentNode)->pGetDof(DISPLACEMENT_Z))->IsFixed();

    for (SizeType j = 0; j < DoF; j++) {

      if (fix_displacements[j]) {
        r_current_acceleration[j] = 0.0;
        r_middle_velocity[j] = 0.0;
      }

      r_current_velocity[j] =
          r_middle_velocity[j] +
          (mTime.Previous - mTime.PreviousMiddle) *
              r_current_acceleration[j]; //+ actual_velocity;
      r_middle_velocity[j] =
          r_current_velocity[j] +
          (mTime.Middle - mTime.Previous) * r_current_acceleration[j];
      r_current_displacement[j] =
          r_current_displacement[j] + mTime.Delta * r_middle_velocity[j];
    } // for DoF
  }

  void UpdateRotationalDegreesOfFreedom(NodeIterator CurrentNode) {
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    SizeType DoF = 3;
    ////// ROTATION DEGRESS OF FREEDOM
    const array_1d<double, 3> &nodal_inertia =
        (CurrentNode)->GetValue(NODAL_INERTIA);
    array_1d<double, 3> &r_current_residual_moment =
        (CurrentNode)->FastGetSolutionStepValue(MOMENT_RESIDUAL);
    array_1d<double, 3> &r_current_angular_velocity =
        (CurrentNode)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
    array_1d<double, 3> &r_current_rotation =
        (CurrentNode)->FastGetSolutionStepValue(ROTATION);
    array_1d<double, 3> &r_middle_angular_velocity =
        (CurrentNode)->GetValue(MIDDLE_ANGULAR_VELOCITY);
    array_1d<double, 3> &r_current_angular_acceleration =
        (CurrentNode)->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

    for (SizeType kk = 0; kk < DoF; ++kk) {
      if (nodal_inertia[kk] > numerical_limit)
        r_current_angular_acceleration[kk] =
            r_current_residual_moment[kk] / nodal_inertia[kk];
      else
        r_current_angular_acceleration[kk] = 0.00;
    }

    bool fix_rotation[3] = {false, false, false};
    fix_rotation[0] = ((CurrentNode)->pGetDof(ROTATION_X))->IsFixed();
    fix_rotation[1] = ((CurrentNode)->pGetDof(ROTATION_Y))->IsFixed();
    fix_rotation[2] = ((CurrentNode)->pGetDof(ROTATION_Z))->IsFixed();

    for (SizeType j = 0; j < DoF; j++) {
      if (fix_rotation[j]) {
        r_current_angular_acceleration[j] = 0.00;
        r_middle_angular_velocity[j] = 0.00;
      }
      r_current_angular_velocity[j] = r_middle_angular_velocity[j] +
                                      (mTime.Previous - mTime.PreviousMiddle) *
                                          r_current_angular_acceleration[j];
      r_middle_angular_velocity[j] =
          r_current_angular_velocity[j] +
          (mTime.Middle - mTime.Previous) * r_current_angular_acceleration[j];
      r_current_rotation[j] =
          r_current_rotation[j] + mTime.Delta * r_middle_angular_velocity[j];
    }
  }

  virtual void SchemeCustomInitialization(ModelPart &rModelPart) {
    KRATOS_TRY
    NodesArrayType &r_nodes = rModelPart.Nodes();
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    auto i_begin = rModelPart.NodesBegin();
    const bool has_dof_for_rot_z = i_begin->HasDofFor(ROTATION_Z);
#pragma omp parallel for firstprivate(i_begin)
    for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
      // Current step information "N+1" (before step update).

      const double &nodal_mass = (i_begin + i)->GetValue(NODAL_MASS);
      array_1d<double, 3> &r_current_residual =
          (i_begin + i)->FastGetSolutionStepValue(FORCE_RESIDUAL);

      array_1d<double, 3> &r_current_velocity =
          (i_begin + i)->FastGetSolutionStepValue(VELOCITY);
      // array_1d<double,3>& r_current_displacement    =
      // i->FastGetSolutionStepValue(DISPLACEMENT);
      array_1d<double, 3> &r_middle_velocity =
          (i_begin + i)->GetValue(MIDDLE_VELOCITY);

      array_1d<double, 3> &r_current_acceleration =
          (i_begin + i)->FastGetSolutionStepValue(ACCELERATION);

      // Solution of the explicit equation:
      if (nodal_mass > numerical_limit)
        r_current_acceleration = r_current_residual / nodal_mass;
      else
        r_current_acceleration = ZeroVector(3);

      SizeType DoF = 3;
      bool fix_displacements[3] = {false, false, false};

      fix_displacements[0] =
          ((i_begin + i)->pGetDof(DISPLACEMENT_X))->IsFixed();
      fix_displacements[1] =
          ((i_begin + i)->pGetDof(DISPLACEMENT_Y))->IsFixed();
      fix_displacements[2] =
          ((i_begin + i)->pGetDof(DISPLACEMENT_Z))->IsFixed();

      for (SizeType j = 0; j < DoF; j++) {
        if (fix_displacements[j]) {
          r_current_acceleration[j] = 0.0;
          r_middle_velocity[j] = 0.0;
        }

        r_middle_velocity[j] =
            0.0 + (mTime.Middle - mTime.Previous) * r_current_acceleration[j];
        r_current_velocity[j] =
            r_middle_velocity[j] +
            (mTime.Previous - mTime.PreviousMiddle) *
                r_current_acceleration[j]; //+ actual_velocity;
        // r_current_displacement[j]  = 0.0;

      } // for DoF
      ////// ROTATION DEGRESS OF FREEDOM
      if (has_dof_for_rot_z) {

        const array_1d<double, 3> &nodal_inertia =
            (i_begin + i)->GetValue(NODAL_INERTIA);
        array_1d<double, 3> &r_current_residual_moment =
            (i_begin + i)->FastGetSolutionStepValue(MOMENT_RESIDUAL);
        array_1d<double, 3> &r_current_angular_velocity =
            (i_begin + i)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
        // array_1d<double,3>& current_rotation                 =
        // i_begin->FastGetSolutionStepValue(ROTATION);
        array_1d<double, 3> &r_middle_angular_velocity =
            (i_begin + i)->GetValue(MIDDLE_ANGULAR_VELOCITY);
        array_1d<double, 3> &r_current_angular_acceleration =
            (i_begin + i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

        for (SizeType kk = 0; kk < DoF; ++kk) {
          if (nodal_inertia[kk] > numerical_limit)
            r_current_angular_acceleration[kk] =
                r_current_residual_moment[kk] / nodal_inertia[kk];
          else
            r_current_angular_acceleration[kk] = 0.00;
        }

        DoF = 3;
        bool fix_rotation[3] = {false, false, false};
        fix_rotation[0] = ((i_begin + i)->pGetDof(ROTATION_X))->IsFixed();
        fix_rotation[1] = ((i_begin + i)->pGetDof(ROTATION_Y))->IsFixed();
        fix_rotation[2] = ((i_begin + i)->pGetDof(ROTATION_Z))->IsFixed();

        for (SizeType j = 0; j < DoF; j++) {
          if (fix_rotation[j]) {
            r_current_angular_acceleration[j] = 0.00;
            r_middle_angular_velocity[j] = 0.00;
          }

          r_middle_angular_velocity[j] = 0.00 +
                                         (mTime.Middle - mTime.Previous) *
                                             r_current_angular_acceleration[j];
          r_current_angular_velocity[j] =
              r_middle_angular_velocity[j] +
              (mTime.Previous - mTime.PreviousMiddle) *
                  r_current_angular_acceleration[j];
          // current_rotation[j]          = 0.00;
        } // trans DoF
      }   // Rot DoF
    }     // for node parallel

    mTime.Previous = mTime.Current;
    mTime.PreviousMiddle = mTime.Middle;
    KRATOS_CATCH("")
  }
  //***************************************************************************
  //***************************************************************************
  void Calculate_RHS_Contribution(Element::Pointer pCurrentElement,
                                  LocalSystemVectorType &RHS_Contribution,
                                  Element::EquationIdVectorType &EquationId,
                                  ProcessInfo &rCurrentProcessInfo) override {
    KRATOS_TRY
    this->TCalculate_RHS_Contribution(pCurrentElement, RHS_Contribution,
                                      rCurrentProcessInfo);
    KRATOS_CATCH("")
  }

  //***************************************************************************
  //***************************************************************************

  void Condition_Calculate_RHS_Contribution(
      Condition::Pointer pCurrentCondition,
      LocalSystemVectorType &RHS_Contribution,
      Element::EquationIdVectorType &EquationId,
      ProcessInfo &rCurrentProcessInfo) override {
    KRATOS_TRY
    this->TCalculate_RHS_Contribution(pCurrentCondition, RHS_Contribution,
                                      rCurrentProcessInfo);
    KRATOS_CATCH("")
  }

  //***************************************************************************
  //***************************************************************************
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

protected:
  /**@name Protected static Member Variables */
  /*@{ */

  struct DeltaTimeParameters {
    double PredictionLevel; // 0, 1, 2

    double Maximum;  // maximum delta time
    double Fraction; // fraction of the delta time
  };

  struct TimeVariables {
    double PreviousMiddle; // n-1/2
    double Previous;       // n
    double Middle;         // n+1/2
    double Current;        // n+1

    double Delta; // time step
  };

  TimeVariables mTime;
  DeltaTimeParameters mDeltaTime;

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

private:
  /**@name Static Member Variables */
  /*@{ */

  /*@} */
  /**@name Member Variables */
  /*@{ */

  /*@} */
  /**@name Private Operators*/
  /*@{ */

  /*@} */
  /**@name Private Operations*/
  /*@{ */
  template <typename TObjectType>
  void TCalculate_RHS_Contribution(TObjectType pCurrentEntity,
                                   LocalSystemVectorType &RHS_Contribution,
                                   ProcessInfo &rCurrentProcessInfo) {
    Matrix dummy_lhs;
    (pCurrentEntity)
        ->CalculateLocalSystem(dummy_lhs, RHS_Contribution,
                               rCurrentProcessInfo);

    (pCurrentEntity)
        ->AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR,
                                  FORCE_RESIDUAL, rCurrentProcessInfo);
    (pCurrentEntity)
        ->AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR,
                                  MOMENT_RESIDUAL, rCurrentProcessInfo);
  }
  /*@} */
  /**@name Private  Access */
  /*@{ */

  /*@} */
  /**@name Private Inquiry */
  /*@{ */

  /*@} */
  /**@name Un accessible methods */
  /*@{ */

  /*@} */

}; /* Class Scheme */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME  defined */
