//
//   Project Name:        KratosPFEMFluidDynamicsApplication $
//   Last modified by:    $Author:                   AFranci $
//   Date:                $Date:                January 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#ifndef KRATOS_TWO_STEP_V_P_DEM_COUPLING_STRATEGY_H
#define KRATOS_TWO_STEP_V_P_DEM_COUPLING_STRATEGY_H

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"
#include "utilities/openmp_utils.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
/* #include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h" */
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

#include "custom_utilities/solver_settings.h"

#include "custom_strategies/strategies/gauss_seidel_linear_strategy.h"

#include "pfem_fluid_dynamics_application_variables.h"
#include "two_step_v_p_DEM_coupling_strategy.h"

#include <stdio.h>
#include <math.h>

namespace Kratos
{

///@addtogroup PFEMFluidDynamicsApplication
///@{

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

template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver>
class TwoStepVPDEMcouplingStrategy : public TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
  ///@name Type Definitions
  ///@{
  KRATOS_CLASS_POINTER_DEFINITION(TwoStepVPDEMcouplingStrategy);

  /// Counted pointer of TwoStepVPDEMcouplingStrategy
  //typedef boost::shared_ptr< TwoStepVPDEMcouplingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> > Pointer;

  typedef TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

  typedef typename BaseType::TDataType TDataType;

  //typedef typename BaseType::DofSetType DofSetType;

  typedef typename BaseType::DofsArrayType DofsArrayType;

  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

  typedef typename BaseType::TSystemVectorType TSystemVectorType;

  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

  typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

  typedef typename TwoStepVPDEMcouplingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

  typedef TwoStepVPSolverSettings<TSparseSpace, TDenseSpace, TLinearSolver> SolverSettingsType;

  using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mVelocityTolerance;
  using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mPressureTolerance;
  using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mMaxPressureIter;
  using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mDomainSize;
  using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mTimeOrder;
  using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mReformDofSet;
  ///@}
  ///@name Life Cycle
  ///@{

  TwoStepVPDEMcouplingStrategy(ModelPart &rModelPart,
                               SolverSettingsType &rSolverConfig) : BaseType(rModelPart)
  {
    TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::InitializeStrategy(rSolverConfig);
  }

  TwoStepVPDEMcouplingStrategy(ModelPart &rModelPart,
                               /*SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& rSolverConfig,*/
                               typename TLinearSolver::Pointer pVelocityLinearSolver,
                               typename TLinearSolver::Pointer pPressureLinearSolver,
                               bool ReformDofSet = true,
                               double VelTol = 0.0001,
                               double PresTol = 0.0001,
                               int MaxPressureIterations = 1, // Only for predictor-corrector
                               unsigned int TimeOrder = 2,
                               unsigned int DomainSize = 2) : BaseType(rModelPart,
                                                                       pVelocityLinearSolver,
                                                                       pPressureLinearSolver,
                                                                       ReformDofSet,
                                                                       VelTol,
                                                                       PresTol,
                                                                       MaxPressureIterations,
                                                                       TimeOrder,
                                                                       DomainSize)
  {
  }

  /// Destructor.
  virtual ~TwoStepVPDEMcouplingStrategy() {}

  void CalculateTemporalVariables() override
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    Vector &BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

    for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
         i != rModelPart.NodesEnd(); ++i)
    {

      array_1d<double, 3> &CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 0);
      array_1d<double, 3> &PreviousVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 1);

      array_1d<double, 3> &CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
      array_1d<double, 3> &PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

      /* if((i)->IsNot(ISOLATED) || (i)->Is(SOLID)){ */
      if ((i)->IsNot(ISOLATED) && ((i)->IsNot(RIGID) || (i)->Is(SOLID)))
      {
        this->UpdateAccelerations(CurrentAcceleration, CurrentVelocity, PreviousAcceleration, PreviousVelocity, BDFcoeffs);
      }
      else if ((i)->Is(RIGID))
      {
        array_1d<double, 3> Zeros(3, 0.0);
        (i)->FastGetSolutionStepValue(ACCELERATION, 0) = Zeros;
        (i)->FastGetSolutionStepValue(ACCELERATION, 1) = Zeros;
      }
      else
      {
        (i)->FastGetSolutionStepValue(PRESSURE, 0) = 0.0;
        (i)->FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
        (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0) = 0.0;
        (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1) = 0.0;
        (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0) = 0.0;
        (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 1) = 0.0;
        if ((i)->SolutionStepsDataHas(VOLUME_ACCELERATION))
        {
          array_1d<double, 3> &VolumeAcceleration = (i)->FastGetSolutionStepValue(VOLUME_ACCELERATION);
          (i)->FastGetSolutionStepValue(ACCELERATION, 0) = VolumeAcceleration;
          (i)->FastGetSolutionStepValue(VELOCITY, 0) += VolumeAcceleration * rCurrentProcessInfo[DELTA_TIME];
        }
      }

      const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      unsigned int timeStep = rCurrentProcessInfo[STEP];
      if (timeStep == 1)
      {
        (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0) = 0;
        (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1) = 0;
        (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0) = 0;
        (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 1) = 0;
      }
      else
      {
        double &CurrentPressure = (i)->FastGetSolutionStepValue(PRESSURE, 0);
        double &PreviousPressure = (i)->FastGetSolutionStepValue(PRESSURE, 1);
        double &CurrentPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0);
        double &CurrentPressureAcceleration = (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0);

        CurrentPressureAcceleration = CurrentPressureVelocity / timeInterval;

        CurrentPressureVelocity = (CurrentPressure - PreviousPressure) / timeInterval;

        CurrentPressureAcceleration += -CurrentPressureVelocity / timeInterval;

        double &previousFluidFraction = (i)->FastGetSolutionStepValue(FLUID_FRACTION_OLD);
        previousFluidFraction = (i)->FastGetSolutionStepValue(FLUID_FRACTION);
      }
    }
  }

  void CalculateDisplacementsAndPorosity() override
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    const double TimeStep = rCurrentProcessInfo[DELTA_TIME];

    for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
         i != rModelPart.NodesEnd(); ++i)
    {

      array_1d<double, 3> &CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 0);
      array_1d<double, 3> &PreviousVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 1);

      array_1d<double, 3> &CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double, 3> &PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);

      const double &currentFluidFraction = (i)->FastGetSolutionStepValue(FLUID_FRACTION);
      const double &previousFluidFraction = (i)->FastGetSolutionStepValue(FLUID_FRACTION_OLD);
      double &currentFluidFractionRate = (i)->FastGetSolutionStepValue(FLUID_FRACTION_RATE);

      /* if( i->IsFixed(DISPLACEMENT_X) == false ) */
      CurrentDisplacement[0] = 0.5 * TimeStep * (CurrentVelocity[0] + PreviousVelocity[0]) + PreviousDisplacement[0];

      /* if( i->IsFixed(DISPLACEMENT_Y) == false ) */
      CurrentDisplacement[1] = 0.5 * TimeStep * (CurrentVelocity[1] + PreviousVelocity[1]) + PreviousDisplacement[1];

      /* if( i->IsFixed(DISPLACEMENT_Z) == false ) */
      CurrentDisplacement[2] = 0.5 * TimeStep * (CurrentVelocity[2] + PreviousVelocity[2]) + PreviousDisplacement[2];

      currentFluidFractionRate = (currentFluidFraction - previousFluidFraction) / TimeStep;
    }
  }

  ///@name Static Member Variables
  ///@{

  ///@}
  ///@name Member Variables
  ///@{
};
} // namespace Kratos
#endif