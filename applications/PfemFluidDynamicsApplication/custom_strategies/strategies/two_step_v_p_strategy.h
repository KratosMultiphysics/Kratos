//
//   Project Name:        KratosPFEMFluidDynamicsApplication $
//   Last modified by:    $Author:                   AFranci $
//   Date:                $Date:                January 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#ifndef KRATOS_TWO_STEP_V_P_STRATEGY_H
#define KRATOS_TWO_STEP_V_P_STRATEGY_H

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
class TwoStepVPStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
  ///@name Type Definitions
  ///@{
  KRATOS_CLASS_POINTER_DEFINITION(TwoStepVPStrategy);

  /// Counted pointer of TwoStepVPStrategy
  //typedef boost::shared_ptr< TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver> > Pointer;

  typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

  typedef typename BaseType::TDataType TDataType;

  //typedef typename BaseType::DofSetType DofSetType;

  typedef typename BaseType::DofsArrayType DofsArrayType;

  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

  typedef typename BaseType::TSystemVectorType TSystemVectorType;

  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

  typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

  typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

  typedef TwoStepVPSolverSettings<TSparseSpace, TDenseSpace, TLinearSolver> SolverSettingsType;

  ///@}
  ///@name Life Cycle
  ///@{

  TwoStepVPStrategy(ModelPart &rModelPart,
                    SolverSettingsType &rSolverConfig) : BaseType(rModelPart)
  {
    InitializeStrategy(rSolverConfig);
  }

  TwoStepVPStrategy(ModelPart &rModelPart,
                    /*SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& rSolverConfig,*/
                    typename TLinearSolver::Pointer pVelocityLinearSolver,
                    typename TLinearSolver::Pointer pPressureLinearSolver,
                    bool ReformDofSet = true,
                    double VelTol = 0.0001,
                    double PresTol = 0.0001,
                    int MaxPressureIterations = 1, // Only for predictor-corrector
                    unsigned int TimeOrder = 2,
                    unsigned int DomainSize = 2) : BaseType(rModelPart), // Move Mesh flag, pass as input?
                                                   mVelocityTolerance(VelTol),
                                                   mPressureTolerance(PresTol),
                                                   mMaxPressureIter(MaxPressureIterations),
                                                   mDomainSize(DomainSize),
                                                   mTimeOrder(TimeOrder),
                                                   mReformDofSet(ReformDofSet)
  {
    KRATOS_TRY;

    BaseType::SetEchoLevel(1);

    // Check that input parameters are reasonable and sufficient.
    this->Check();

    bool CalculateNormDxFlag = true;

    bool ReformDofAtEachIteration = false; // DofSet modifiaction is managed by the fractional step strategy, auxiliary strategies should not modify the DofSet directly.

    // Additional Typedefs
    //typedef typename Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3 > > > VarComponent;
    typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    //initializing fractional velocity solution step
    typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;
    typename SchemeType::Pointer pScheme;

    typename SchemeType::Pointer Temp = typename SchemeType::Pointer(new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>());
    pScheme.swap(Temp);

    //CONSTRUCTION OF VELOCITY
    BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pVelocityLinearSolver));
    /* BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (pVelocityLinearSolver)); */

    this->mpMomentumStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pVelocityLinearSolver, vel_build, ReformDofAtEachIteration, CalculateNormDxFlag));

    this->mpMomentumStrategy->SetEchoLevel(BaseType::GetEchoLevel());

    vel_build->SetCalculateReactionsFlag(false);

    /* BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pPressureLinearSolver, PRESSURE)); */
    BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pPressureLinearSolver));

    this->mpPressureStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pPressureLinearSolver, pressure_build, ReformDofAtEachIteration, CalculateNormDxFlag));

    this->mpPressureStrategy->SetEchoLevel(BaseType::GetEchoLevel());

    pressure_build->SetCalculateReactionsFlag(false);

    KRATOS_CATCH("");
  }

  /// Destructor.
  virtual ~TwoStepVPStrategy() {}

  int Check() override
  {
    KRATOS_TRY;

    // Check elements and conditions in the model part
    int ierr = BaseType::Check();
    if (ierr != 0)
      return ierr;

    if (DELTA_TIME.Key() == 0)
      KRATOS_THROW_ERROR(std::runtime_error, "DELTA_TIME Key is 0. Check that the application was correctly registered.", "");
    if (BDF_COEFFICIENTS.Key() == 0)
      KRATOS_THROW_ERROR(std::runtime_error, "BDF_COEFFICIENTS Key is 0. Check that the application was correctly registered.", "");

    ModelPart &rModelPart = BaseType::GetModelPart();

    if (mTimeOrder == 2 && rModelPart.GetBufferSize() < 3)
      KRATOS_THROW_ERROR(std::invalid_argument, "Buffer size too small for fractional step strategy (BDF2), needed 3, got ", rModelPart.GetBufferSize());
    if (mTimeOrder == 1 && rModelPart.GetBufferSize() < 2)
      KRATOS_THROW_ERROR(std::invalid_argument, "Buffer size too small for fractional step strategy (Backward Euler), needed 2, got ", rModelPart.GetBufferSize());

    const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();

    for (ModelPart::ElementIterator itEl = rModelPart.ElementsBegin(); itEl != rModelPart.ElementsEnd(); ++itEl)
    {
      ierr = itEl->Check(rCurrentProcessInfo);
      if (ierr != 0)
        break;
    }

    /* for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond) */
    /* { */
    /*     ierr = itCond->Check(rCurrentProcessInfo); */
    /*     if (ierr != 0) break; */
    /* } */

    return ierr;

    KRATOS_CATCH("");
  }

  double Solve() override
  {
    // Initialize BDF2 coefficients
    ModelPart &rModelPart = BaseType::GetModelPart();
    this->SetTimeCoefficients(rModelPart.GetProcessInfo());
    double NormDp = 0.0;
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    double currentTime = rCurrentProcessInfo[TIME];
    double timeInterval = rCurrentProcessInfo[DELTA_TIME];
    bool timeIntervalChanged = rCurrentProcessInfo[TIME_INTERVAL_CHANGED];
    unsigned int stepsWithChangedDt = rCurrentProcessInfo[STEPS_WITH_CHANGED_DT];

    unsigned int maxNonLinearIterations = mMaxPressureIter;

    KRATOS_INFO("TwoStepVPStrategy") << "\n                   Solve with two_step_vp strategy at t=" << currentTime << "s" << std::endl;

    if ((timeIntervalChanged == true && currentTime > 10 * timeInterval) || stepsWithChangedDt > 0)
    {
      maxNonLinearIterations *= 2;
    }
    if (currentTime < 10 * timeInterval)
    {
      if (BaseType::GetEchoLevel() > 1)
        std::cout << "within the first 10 time steps, I consider the given iteration number x3" << std::endl;
      maxNonLinearIterations *= 3;
    }
    if (currentTime < 20 * timeInterval && currentTime >= 10 * timeInterval)
    {
      if (BaseType::GetEchoLevel() > 1)
        std::cout << "within the second 10 time steps, I consider the given iteration number x2" << std::endl;
      maxNonLinearIterations *= 2;
    }
    bool momentumConverged = true;
    bool continuityConverged = false;
    bool fixedTimeStep = false;

    bool momentumAlreadyConverged = false;
    bool continuityAlreadyConverged = false;
    /* boost::timer solve_step_time; */

    // Iterative solution for pressure
    /* unsigned int timeStep = rCurrentProcessInfo[STEP]; */
    /* if(timeStep==1){ */
    /* 	unsigned int iter=0; */
    /* 	continuityConverged = this->SolveContinuityIteration(iter,maxNonLinearIterations); */
    /* }else if(timeStep==2){ */
    /* 	unsigned int iter=0; */
    /* 	momentumConverged = this->SolveMomentumIteration(iter,maxNonLinearIterations,fixedTimeStep); */
    /* }else{ */

    this->UnactiveSliverElements();

    for (unsigned int it = 0; it < maxNonLinearIterations; ++it)
    {
      if (BaseType::GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)
        std::cout << "----- > iteration: " << it << std::endl;

      momentumConverged = this->SolveMomentumIteration(it, maxNonLinearIterations, fixedTimeStep);

      this->UpdateTopology(rModelPart, BaseType::GetEchoLevel());

      if ((momentumConverged == true || it == maxNonLinearIterations - 1) && momentumAlreadyConverged == false)
      {
        // std::ofstream myfile;
        // myfile.open ("momentumConvergedIteration.txt",std::ios::app);
        // myfile << currentTime << "\t" << it << "\n";
        // myfile.close();
        momentumAlreadyConverged = true;
      }
      if ((continuityConverged == true || it == maxNonLinearIterations - 1) && continuityAlreadyConverged == false)
      {
        // std::ofstream myfile;
        // myfile.open ("continuityConvergedIteration.txt",std::ios::app);
        // myfile << currentTime << "\t" << it << "\n";
        // myfile.close();
        continuityAlreadyConverged = true;
      }

      if (fixedTimeStep == false)
      {
        continuityConverged = this->SolveContinuityIteration(it, maxNonLinearIterations);
      }
      if (it == maxNonLinearIterations - 1 || ((continuityConverged && momentumConverged) && it > 2))
      {
        //this->ComputeErrorL2Norm();
        //this->ComputeErrorL2NormCasePoiseuille();
        this->UpdateStressStrain();
        // std::ofstream myfile;
        // myfile.open ("maxConvergedIteration.txt",std::ios::app);
        // myfile << currentTime << "\t" << it << "\n";
        // myfile.close();
      }

      if ((continuityConverged && momentumConverged) && it > 2)
      {
        rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, false);
        rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, false);

        KRATOS_INFO("TwoStepVPStrategy") << "V-P strategy converged in " << it + 1 << " iterations." << std::endl;

        break;
      }
      if (fixedTimeStep == true)
      {
        break;
      }
    }

    /* } */

    if (!continuityConverged && !momentumConverged && BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
      std::cout << "Convergence tolerance not reached." << std::endl;

    /* std::cout << "solve_step_time : " << solve_step_time.elapsed() << std::endl; */

    if (mReformDofSet)
      this->Clear();

    return NormDp;
  }

  void FinalizeSolutionStep() override
  {
    /* this->UpdateStressStrain(); */
  }

  void InitializeSolutionStep() override
  {
  }

  void UpdateTopology(ModelPart &rModelPart, unsigned int echoLevel)
  {
    KRATOS_TRY;

    this->CalculateDisplacementsAndPorosity();
    BaseType::MoveMesh();
    /* BoundaryNormalsCalculationUtilities BoundaryComputation; */
    /* BoundaryComputation.CalculateWeightedBoundaryNormals(rModelPart, echoLevel); */

    KRATOS_CATCH("");
  }

  void UnactiveSliverElements()
  {
    KRATOS_TRY;

    ModelPart &rModelPart = BaseType::GetModelPart();
    const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
    MesherUtilities MesherUtils;
    double ModelPartVolume = MesherUtils.ComputeModelPartVolume(rModelPart);
    double CriticalVolume = 0.001 * ModelPartVolume / double(rModelPart.Elements().size());
    double ElementalVolume = 0;

#pragma omp parallel
    {
      ModelPart::ElementIterator ElemBegin;
      ModelPart::ElementIterator ElemEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElemBegin, ElemEnd);
      for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
      {
        unsigned int numNodes = itElem->GetGeometry().size();
        if (numNodes == (dimension + 1))
        {
          if (dimension == 2)
          {
            ElementalVolume = (itElem)->GetGeometry().Area();
          }
          else if (dimension == 3)
          {
            ElementalVolume = (itElem)->GetGeometry().Volume();
          }

          if (ElementalVolume < CriticalVolume)
          {
            // std::cout << "sliver element: it has Volume: " << ElementalVolume << " vs CriticalVolume(meanVol/1000): " << CriticalVolume<< std::endl;
            (itElem)->Set(ACTIVE, false);
          }
          else
          {
            (itElem)->Set(ACTIVE, true);
          }
        }
      }
    }
    KRATOS_CATCH("");
  }

  void CalculatePressureVelocity()
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
    unsigned int timeStep = rCurrentProcessInfo[STEP];

    for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
         i != rModelPart.NodesEnd(); ++i)
    {
      if (timeStep == 1)
      {
        (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0) = 0;
        (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1) = 0;
      }
      else
      {
        double &CurrentPressure = (i)->FastGetSolutionStepValue(PRESSURE, 0);
        double &PreviousPressure = (i)->FastGetSolutionStepValue(PRESSURE, 1);
        double &CurrentPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0);
        CurrentPressureVelocity = (CurrentPressure - PreviousPressure) / timeInterval;
      }
    }
  }

  void CalculatePressureAcceleration()
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
    unsigned int timeStep = rCurrentProcessInfo[STEP];

    for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
         i != rModelPart.NodesEnd(); ++i)
    {
      if (timeStep == 1)
      {
        (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0) = 0;
        (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 1) = 0;
      }
      else
      {
        double &CurrentPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0);
        double &PreviousPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1);
        double &CurrentPressureAcceleration = (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0);
        CurrentPressureAcceleration = (CurrentPressureVelocity - PreviousPressureVelocity) / timeInterval;
      }
    }
  }

  virtual void CalculateTemporalVariables()
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
        UpdateAccelerations(CurrentAcceleration, CurrentVelocity, PreviousAcceleration, PreviousVelocity, BDFcoeffs);
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
      }
    }
  }

  void CalculateAccelerations()
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
        UpdateAccelerations(CurrentAcceleration, CurrentVelocity, PreviousAcceleration, PreviousVelocity, BDFcoeffs);
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
    }
  }

  inline void UpdateAccelerations(array_1d<double, 3> &CurrentAcceleration,
                                  const array_1d<double, 3> &CurrentVelocity,
                                  array_1d<double, 3> &PreviousAcceleration,
                                  const array_1d<double, 3> &PreviousVelocity,
                                  Vector &BDFcoeffs)
  {
    /* noalias(PreviousAcceleration)=CurrentAcceleration; */
    noalias(CurrentAcceleration) = -BDFcoeffs[1] * (CurrentVelocity - PreviousVelocity) - PreviousAcceleration;
    // std::cout<<"rBDFCoeffs[0] is "<<rBDFCoeffs[0]<<std::endl;//3/(2*delta_t)
    // std::cout<<"rBDFCoeffs[1] is "<<rBDFCoeffs[1]<<std::endl;//-2/(delta_t)
    // std::cout<<"rBDFCoeffs[2] is "<<rBDFCoeffs[2]<<std::endl;//1/(2*delta_t)
  }

  virtual void CalculateDisplacementsAndPorosity()
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

      /* if( i->IsFixed(DISPLACEMENT_X) == false ) */
      CurrentDisplacement[0] = 0.5 * TimeStep * (CurrentVelocity[0] + PreviousVelocity[0]) + PreviousDisplacement[0];

      /* if( i->IsFixed(DISPLACEMENT_Y) == false ) */
      CurrentDisplacement[1] = 0.5 * TimeStep * (CurrentVelocity[1] + PreviousVelocity[1]) + PreviousDisplacement[1];

      /* if( i->IsFixed(DISPLACEMENT_Z) == false ) */
      CurrentDisplacement[2] = 0.5 * TimeStep * (CurrentVelocity[2] + PreviousVelocity[2]) + PreviousDisplacement[2];

      // currentFluidFractionRate = (currentFluidFraction - previousFluidFraction)/TimeStep;
    }
  }

  void UpdateStressStrain()
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();

#pragma omp parallel
    {
      ModelPart::ElementIterator ElemBegin;
      ModelPart::ElementIterator ElemEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElemBegin, ElemEnd);
      for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
      {
        /* itElem-> InitializeElementStrainStressState(); */
        itElem->InitializeSolutionStep(rCurrentProcessInfo);
      }
    }

    /* this->CalculateAccelerations(); */
    /* this->CalculatePressureVelocity(); */
    /* this->CalculatePressureAcceleration(); */

    this->CalculateTemporalVariables();
  }

  void Clear() override
  {
    mpMomentumStrategy->Clear();
    mpPressureStrategy->Clear();
  }

  ///@}
  ///@name Access
  ///@{

  void SetEchoLevel(int Level) override
  {
    BaseType::SetEchoLevel(Level);
    int StrategyLevel = Level > 0 ? Level - 1 : 0;
    mpMomentumStrategy->SetEchoLevel(StrategyLevel);
    mpPressureStrategy->SetEchoLevel(StrategyLevel);
  }

  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  std::string Info() const override
  {
    std::stringstream buffer;
    buffer << "TwoStepVPStrategy";
    return buffer.str();
  }

  /// Print information about this object.
  void PrintInfo(std::ostream &rOStream) const override
  {
    rOStream << "TwoStepVPStrategy";
  }

  /// Print object's data.
  void PrintData(std::ostream &rOStream) const override
  {
  }

  ///@}
  ///@name Friends
  ///@{

  ///@}

protected:
  ///@name Protected Life Cycle
  ///@{

  ///@}
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

  /// Calculate the coefficients for time iteration.
  /**
     * @param rCurrentProcessInfo ProcessInfo instance from the fluid ModelPart. Must contain DELTA_TIME and BDF_COEFFICIENTS variables.
     */
  void SetTimeCoefficients(ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    if (mTimeOrder == 2)
    {
      //calculate the BDF coefficients
      double Dt = rCurrentProcessInfo[DELTA_TIME];
      double OldDt = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

      double Rho = OldDt / Dt;
      double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

      Vector &BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
      BDFcoeffs.resize(3, false);

      BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho);        //coefficient for step n+1 (3/2Dt if Dt is constant)
      BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
      BDFcoeffs[2] = TimeCoeff;                                  //coefficient for step n-1 (1/2Dt if Dt is constant)
    }
    else if (mTimeOrder == 1)
    {
      double Dt = rCurrentProcessInfo[DELTA_TIME];
      double TimeCoeff = 1.0 / Dt;

      Vector &BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
      BDFcoeffs.resize(2, false);

      BDFcoeffs[0] = TimeCoeff;  //coefficient for step n+1 (1/Dt)
      BDFcoeffs[1] = -TimeCoeff; //coefficient for step n (-1/Dt)
    }

    KRATOS_CATCH("");
  }

  bool SolveMomentumIteration(unsigned int it, unsigned int maxIt, bool &fixedTimeStep)
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    int Rank = rModelPart.GetCommunicator().MyPID();
    bool ConvergedMomentum = false;
    double NormDv = 0;
    fixedTimeStep = false;
    // build momentum system and solve for fractional step velocity increment
    rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP, 1);

    /* std::cout<<"---- m o m e n t u m   e q u a t i o n s ----"<<std::endl; */
    if (it == 0)
    {
      mpMomentumStrategy->InitializeSolutionStep();
    }
    /* else{ */
    /* 	NormDv = mpMomentumStrategy->Solve(); */
    /* } */
    NormDv = mpMomentumStrategy->Solve();

    if (BaseType::GetEchoLevel() > 1 && Rank == 0)
      std::cout << "-------------- s o l v e d ! ------------------" << std::endl;

    double DvErrorNorm = 0;
    ConvergedMomentum = this->CheckVelocityConvergence(NormDv, DvErrorNorm);

    unsigned int iterationForCheck = 2;
    KRATOS_INFO("TwoStepVPStrategy") << "iteration(" << it << ") Velocity error: " << DvErrorNorm << " velTol: " << mVelocityTolerance << std::endl;

    // Check convergence
    if (it == maxIt - 1)
    {

      KRATOS_INFO("TwoStepVPStrategy") << "iteration(" << it << ") Final Velocity error: " << DvErrorNorm << " velTol: " << mVelocityTolerance << std::endl;

      fixedTimeStep = this->FixTimeStepMomentum(DvErrorNorm);
    }
    else if (it > iterationForCheck)
    {
      fixedTimeStep = this->CheckMomentumConvergence(DvErrorNorm);
    }

    // ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
    // double currentTime = rCurrentProcessInfo[TIME];
    // double tolerance=0.0000000001;
    // if(currentTime>(0.25-tolerance) && currentTime<(0.25+tolerance)){
    // 	std::ofstream myfile;
    //   myfile.open ("velocityConvergenceAt025s.txt",std::ios::app);
    // 	myfile << it << "\t" << DvErrorNorm << "\n";
    //   myfile.close();
    // }
    // else if(currentTime>(0.5-tolerance) && currentTime<(0.5+tolerance)){
    // 	std::ofstream myfile;
    //   myfile.open ("velocityConvergenceAt05s.txt",std::ios::app);
    // 	myfile << it << "\t" << DvErrorNorm << "\n";
    //   myfile.close();
    // }
    // else if(currentTime>(0.75-tolerance) && currentTime<(0.75+tolerance)){
    // 	std::ofstream myfile;
    //   myfile.open ("velocityConvergenceAt075s.txt",std::ios::app);
    // 	myfile << it << "\t" << DvErrorNorm << "\n";
    //   myfile.close();
    // }
    // else if(currentTime>(1.0-tolerance) && currentTime<(1.0+tolerance)){
    // 	std::ofstream myfile;
    //   myfile.open ("velocityConvergenceAt100s.txt",std::ios::app);
    // 	myfile << it << "\t" << DvErrorNorm << "\n";
    //   myfile.close();
    // }

    if (!ConvergedMomentum && BaseType::GetEchoLevel() > 0 && Rank == 0)
      std::cout << "Momentum equations did not reach the convergence tolerance." << std::endl;

    return ConvergedMomentum;
  }

  bool SolveContinuityIteration(unsigned int it, unsigned int maxIt)
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    int Rank = rModelPart.GetCommunicator().MyPID();
    bool ConvergedContinuity = false;
    double NormDp = 0;

    // 2. Pressure solution
    rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP, 5);

    /* std::cout<<"     ---- c o n t i n u i t y   e q u a t i o n ----"<<std::endl; */

    if (it == 0)
    {
      mpPressureStrategy->InitializeSolutionStep();
    }
    /* else{ */
    /* 	NormDp = mpPressureStrategy->Solve(); */
    /* } */
    NormDp = mpPressureStrategy->Solve();

    if (BaseType::GetEchoLevel() > 0 && Rank == 0)
      std::cout << "The norm of pressure is: " << NormDp << std::endl;

    double DpErrorNorm = 0;
    ConvergedContinuity = this->CheckPressureConvergence(NormDp, DpErrorNorm);
    KRATOS_INFO("TwoStepVPStrategy") << "                    iteration(" << it << ") Pressure error: " << DpErrorNorm << " presTol: " << mPressureTolerance << std::endl;

    // Check convergence
    if (it == maxIt - 1)
    {

      KRATOS_INFO("TwoStepVPStrategy") << "       iteration(" << it << ") Final Pressure error: " << DpErrorNorm << " presTol: " << mPressureTolerance << std::endl;

      ConvergedContinuity = this->FixTimeStepContinuity(DpErrorNorm);
    }

    // ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
    // double currentTime = rCurrentProcessInfo[TIME];
    // double tolerance=0.0000000001;
    // if(currentTime>(0.25-tolerance) && currentTime<(0.25+tolerance)){
    // 	std::ofstream myfile;
    //  myfile.open ("pressureConvergenceAt025s.txt",std::ios::app);
    // 	myfile << it << "\t" << DpErrorNorm << "\n";
    //  myfile.close();
    // }
    // else if(currentTime>(0.5-tolerance) && currentTime<(0.5+tolerance)){
    // 	std::ofstream myfile;
    //  myfile.open ("pressureConvergenceAt05s.txt",std::ios::app);
    // 	myfile << it << "\t" << DpErrorNorm << "\n";
    //  myfile.close();
    // }
    // else if(currentTime>(0.75-tolerance) && currentTime<(0.75+tolerance)){
    // 	std::ofstream myfile;
    //  myfile.open ("pressureConvergenceAt075s.txt",std::ios::app);
    // 	myfile << it << "\t" << DpErrorNorm << "\n";
    //  myfile.close();
    // }
    // else if(currentTime>(1.0-tolerance) && currentTime<(1.0+tolerance)){
    // 	std::ofstream myfile;
    //  myfile.open ("pressureConvergenceAt100s.txt",std::ios::app);
    // 	myfile << it << "\t" << DpErrorNorm << "\n";
    //  myfile.close();
    // }

    if (!ConvergedContinuity && BaseType::GetEchoLevel() > 0 && Rank == 0)
      std::cout << "Continuity equation did not reach the convergence tolerance." << std::endl;

    return ConvergedContinuity;
  }

  void ComputeErrorL2Norm()
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    const double currentTime = rCurrentProcessInfo[TIME];
    const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

    long double sumErrorL2Velocity = 0;
    long double sumErrorL2VelocityX = 0;
    long double sumErrorL2VelocityY = 0;
    long double sumErrorL2Pressure = 0;
    long double sumErrorL2TauXX = 0;
    long double sumErrorL2TauYY = 0;
    long double sumErrorL2TauXY = 0;

#pragma omp parallel
    {
      ModelPart::ElementIterator ElemBegin;
      ModelPart::ElementIterator ElemEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElemBegin, ElemEnd);
      for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
      {

        Element::GeometryType &geometry = itElem->GetGeometry();
        long double nodalArea = 0;

        if (dimension == 2)
        {
          nodalArea = geometry.Area() / 3.0;
        }
        else if (dimension == 3)
        {
          nodalArea = geometry.Volume() * 0.25;
        }

        long double bariPosX = 0;
        long double bariPosY = 0;

        long double eleErrorL2Velocity = 0;
        long double eleErrorL2VelocityX = 0;
        long double eleErrorL2VelocityY = 0;
        long double eleErrorL2Pressure = 0;

        //ShapeFunctionDerivativesArrayType DN_DX;
        Matrix NContainer;
        NContainer = geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
        //this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);

        const Vector &N = row(NContainer, 0);
        //  itElem->EvaluateInPoint(elementalPressure,PRESSURE,N);
        const unsigned int NumNodes = geometry.size();

        double elementalPressure = N[0] * geometry(0)->FastGetSolutionStepValue(PRESSURE);
        double elementalVelocityX = N[0] * geometry(0)->FastGetSolutionStepValue(VELOCITY_X);
        double elementalVelocityY = N[0] * geometry(0)->FastGetSolutionStepValue(VELOCITY_Y);
        ;

        for (unsigned int i = 1; i < NumNodes; i++)
        {
          elementalPressure += N[i] * geometry(i)->FastGetSolutionStepValue(PRESSURE);
          elementalVelocityX += N[i] * geometry(i)->FastGetSolutionStepValue(VELOCITY_X);
          elementalVelocityY += N[i] * geometry(i)->FastGetSolutionStepValue(VELOCITY_Y);
        }

        for (unsigned int i = 0; i < geometry.size(); i++)
        {

          // index = i*dimension;
          const long double nodalPosX = geometry(i)->X();
          const long double nodalPosY = geometry(i)->Y();
          // const long double velX      = geometry(i)->FastGetSolutionStepValue(VELOCITY_X);
          // const long double velY      = geometry(i)->FastGetSolutionStepValue(VELOCITY_Y);
          // const long double pressure  = geometry(i)->FastGetSolutionStepValue(PRESSURE);

          // long double expectedVelocityX =  pow(posX,2) * (1.0-posX)*(1.0-posX) * ( 2.0*posY - 6.0*pow(posY,2) + 4.0*pow(posY,3) );
          // long double expectedVelocityY = -pow(posY,2) * (1.0-posY)*(1.0-posY) * ( 2.0*posX - 6.0*pow(posX,2) + 4.0*pow(posX,3) );
          // long double expectedPressure  = -posX * (1.0-posX);

          // long double nodalErrorVelocityX = velX - expectedVelocityX;
          // long double nodalErrorVelocityY = velY - expectedVelocityY;
          // long double nodalErrorPressure  = pressure - expectedPressure;

          // sumErrorL2Velocity  +=  (pow(nodalErrorVelocityX,2) + pow(nodalErrorVelocityY,2)) * nodalArea;
          // sumErrorL2VelocityX +=  pow(nodalErrorVelocityX,2) * nodalArea;
          // sumErrorL2VelocityY +=  pow(nodalErrorVelocityY,2) * nodalArea;
          // sumErrorL2Pressure  +=  pow(nodalErrorPressure,2)  * nodalArea;
          // eleErrorL2Velocity  +=  pow(nodalErrorVelocityX,2) + pow(nodalErrorVelocityY,2);
          // eleErrorL2VelocityX +=  pow(nodalErrorVelocityX,2);
          // eleErrorL2VelocityY +=  pow(nodalErrorVelocityY,2);
          // eleErrorL2Pressure  +=  pow(nodalErrorPressure,2);

          bariPosX += nodalPosX / 3.0;
          bariPosY += nodalPosY / 3.0;
        }

        const long double posX = bariPosX;
        const long double posY = bariPosY;
        long double expectedVelocityX = pow(posX, 2) * (1.0 - posX) * (1.0 - posX) * (2.0 * posY - 6.0 * pow(posY, 2) + 4.0 * pow(posY, 3));
        long double expectedVelocityY = -pow(posY, 2) * (1.0 - posY) * (1.0 - posY) * (2.0 * posX - 6.0 * pow(posX, 2) + 4.0 * pow(posX, 3));
        long double expectedPressure = -posX * (1.0 - posX);

        eleErrorL2VelocityX = elementalVelocityX - expectedVelocityX;
        eleErrorL2VelocityY = elementalVelocityY - expectedVelocityY;
        eleErrorL2Pressure = elementalPressure - expectedPressure;

        sumErrorL2VelocityX += pow(eleErrorL2VelocityX, 2) * geometry.Area();
        sumErrorL2VelocityY += pow(eleErrorL2VelocityY, 2) * geometry.Area();
        sumErrorL2Pressure += pow(eleErrorL2Pressure, 2) * geometry.Area();

        // sumErrorL2Velocity  +=  eleErrorL2Velocity * geometry.Area();
        // sumErrorL2VelocityX +=  eleErrorL2VelocityX * geometry.Area();
        // sumErrorL2VelocityY +=  eleErrorL2VelocityY * geometry.Area();
        // sumErrorL2Pressure  +=  eleErrorL2Pressure * geometry.Area();

        const long double tauXX = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_XX);
        const long double tauYY = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_YY);
        const long double tauXY = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_XY);

        long double expectedTauXX = 2.0 * (-4.0 * (1.0 - bariPosX) * bariPosX * (-1.0 + 2.0 * bariPosX) * bariPosY * (1.0 - 3.0 * bariPosY + 2.0 * pow(bariPosY, 2)));
        long double expectedTauYY = 2.0 * (4.0 * bariPosX * (1.0 - 3.0 * bariPosX + 2.0 * pow(bariPosX, 2)) * (1.0 - bariPosY) * bariPosY * (-1.0 + 2.0 * bariPosY));
        long double expectedTauXY = (2.0 * (1.0 - 6.0 * bariPosY + 6.0 * pow(bariPosY, 2)) * (1.0 - bariPosX) * (1.0 - bariPosX) * pow(bariPosX, 2) - 2.0 * (1.0 - 6.0 * bariPosX + 6.0 * pow(bariPosX, 2)) * (1.0 - bariPosY) * (1 - bariPosY) * pow(bariPosY, 2));

        long double nodalErrorTauXX = tauXX - expectedTauXX;
        long double nodalErrorTauYY = tauYY - expectedTauYY;
        long double nodalErrorTauXY = tauXY - expectedTauXY;

        // std::cout<<"tauXX "<<tauXX<<"     expectedtauXX "<<expectedTauXX<<"     nodalErrorTauXX "<<nodalErrorTauXX<<std::endl;
        // std::cout<<"tauyy "<<tauYY<<"     expectedtauYY "<<expectedTauYY<<"     nodalErrorTauYY "<<nodalErrorTauYY<<std::endl;
        // std::cout<<"tauXY "<<tauXY<<"     expectedtauXY "<<expectedTauXY<<"     nodalErrorTauXY "<<nodalErrorTauXY<<std::endl;

        sumErrorL2TauXX += pow(nodalErrorTauXX, 2) * geometry.Area();
        sumErrorL2TauYY += pow(nodalErrorTauYY, 2) * geometry.Area();
        sumErrorL2TauXY += pow(nodalErrorTauXY, 2) * geometry.Area();
      }
    }

    // long double errorL2Velocity  = sumErrorL2Velocity;
    // long double errorL2VelocityX = sumErrorL2VelocityX;
    // long double errorL2VelocityY = sumErrorL2VelocityY;
    // long double errorL2Pressure  = sumErrorL2Pressure;
    long double errorL2Velocity = sqrt(sumErrorL2Velocity);
    long double errorL2VelocityX = sqrt(sumErrorL2VelocityX);
    long double errorL2VelocityY = sqrt(sumErrorL2VelocityY);
    long double errorL2Pressure = sqrt(sumErrorL2Pressure);
    long double errorL2TauXX = sqrt(sumErrorL2TauXX);
    long double errorL2TauYY = sqrt(sumErrorL2TauYY);
    long double errorL2TauXY = sqrt(sumErrorL2TauXY);

    std::ofstream myfileVelocity;
    myfileVelocity.open("errorL2VelocityFile.txt", std::ios::app);
    myfileVelocity << currentTime << "\t" << errorL2Velocity << "\n";
    myfileVelocity.close();

    std::ofstream myfileVelocityX;
    myfileVelocityX.open("errorL2VelocityXFile.txt", std::ios::app);
    myfileVelocityX << currentTime << "\t" << errorL2VelocityX << "\n";
    myfileVelocityX.close();

    std::ofstream myfileVelocityY;
    myfileVelocityY.open("errorL2VelocityYFile.txt", std::ios::app);
    myfileVelocityY << currentTime << "\t" << errorL2VelocityY << "\n";
    myfileVelocityY.close();

    std::ofstream myfilePressure;
    myfilePressure.open("errorL2PressureFile.txt", std::ios::app);
    myfilePressure << currentTime << "\t" << errorL2Pressure << "\n";
    myfilePressure.close();

    std::ofstream myfileTauXX;
    myfileTauXX.open("errorL2TauXXFile.txt", std::ios::app);
    myfileTauXX << currentTime << "\t" << errorL2TauXX << "\n";
    myfileTauXX.close();

    std::ofstream myfileTauYY;
    myfileTauYY.open("errorL2TauYYFile.txt", std::ios::app);
    myfileTauYY << currentTime << "\t" << errorL2TauYY << "\n";
    myfileTauYY.close();

    std::ofstream myfileTauXY;
    myfileTauXY.open("errorL2TauXYFile.txt", std::ios::app);
    myfileTauXY << currentTime << "\t" << errorL2TauXY << "\n";
    myfileTauXY.close();
  }

  void ComputeErrorL2NormCasePoiseuille()
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    const double currentTime = rCurrentProcessInfo[TIME];
    const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

    double sumErrorL2VelocityTheta = 0;
    double sumErrorL2TauTheta = 0;

    double r_in = 0.2;
    double R_out = 0.5;
    double kappa = r_in / R_out;
    double omega = 0.5;
    double viscosity = 100.0;

#pragma omp parallel
    {
      ModelPart::ElementIterator ElemBegin;
      ModelPart::ElementIterator ElemEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElemBegin, ElemEnd);
      for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
      {

        Element::GeometryType &geometry = itElem->GetGeometry();
        long double nodalArea = 0;

        if (dimension == 2)
        {
          nodalArea = geometry.Area() / 3.0;
        }
        else if (dimension == 3)
        {
          nodalArea = geometry.Volume() * 0.25;
        }

        long double bariPosX = 0;
        long double bariPosY = 0;

        long double eleErrorL2Velocity = 0;
        long double eleErrorL2VelocityX = 0;
        long double eleErrorL2VelocityY = 0;
        long double eleErrorL2Pressure = 0;

        //ShapeFunctionDerivativesArrayType DN_DX;
        Matrix NContainer;
        NContainer = geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
        //this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);

        const Vector &N = row(NContainer, 0);
        //  itElem->EvaluateInPoint(elementalPressure,PRESSURE,N);
        const unsigned int NumNodes = geometry.size();

        double elementalPressure = N[0] * geometry(0)->FastGetSolutionStepValue(PRESSURE);
        double elementalVelocityX = N[0] * geometry(0)->FastGetSolutionStepValue(VELOCITY_X);
        double elementalVelocityY = N[0] * geometry(0)->FastGetSolutionStepValue(VELOCITY_Y);
        ;

        for (unsigned int i = 1; i < NumNodes; i++)
        {
          elementalPressure += N[i] * geometry(i)->FastGetSolutionStepValue(PRESSURE);
          elementalVelocityX += N[i] * geometry(i)->FastGetSolutionStepValue(VELOCITY_X);
          elementalVelocityY += N[i] * geometry(i)->FastGetSolutionStepValue(VELOCITY_Y);
        }

        for (unsigned int i = 0; i < geometry.size(); i++)
        {

          // index = i*dimension;
          const long double nodalPosX = geometry(i)->X();
          const long double nodalPosY = geometry(i)->Y();

          bariPosX += nodalPosX / 3.0;
          bariPosY += nodalPosY / 3.0;
        }

        const long double posX = bariPosX;
        const long double posY = bariPosY;
        const double rPos = sqrt(pow(posX, 2) + pow(posY, 2));
        const double cosalfa = posX / rPos;
        const double sinalfa = posY / rPos;
        const double sin2alfa = 2.0 * cosalfa * sinalfa;
        const double cos2alfa = 1.0 - 2.0 * pow(sinalfa, 2);

        double expectedVelocityTheta = pow(kappa, 2) * omega * R_out / (1.0 - pow(kappa, 2)) * (R_out / rPos - rPos / R_out);
        double computedVelocityTheta = sqrt(pow(elementalVelocityX, 2) + pow(elementalVelocityY, 2));
        double nodalErrorVelocityTheta = computedVelocityTheta - expectedVelocityTheta;

        const long double tauXX = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_XX);
        const long double tauYY = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_YY);
        const long double tauXY = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_XY);

        double expectedTauTheta = (2.0 * viscosity * pow(kappa, 2) * omega * pow(R_out, 2)) / (1.0 - pow(kappa, 2)) / pow(rPos, 2);
        double computedTauTheta = (tauXX - tauYY) * sin2alfa / 2.0 - tauXY * cos2alfa;
        double nodalErrorTauTheta = computedTauTheta - expectedTauTheta;

        sumErrorL2VelocityTheta += pow(nodalErrorVelocityTheta, 2) * geometry.Area();
        sumErrorL2TauTheta += pow(nodalErrorTauTheta, 2) * geometry.Area();
      }
    }

    double errorL2VelocityTheta = sqrt(sumErrorL2VelocityTheta);
    double errorL2TauTheta = sqrt(sumErrorL2TauTheta);

    std::ofstream myfileVelocity;
    myfileVelocity.open("errorL2Poiseuille.txt", std::ios::app);
    myfileVelocity << currentTime << "\t" << errorL2VelocityTheta << "\t" << errorL2TauTheta << "\n";
    myfileVelocity.close();
  }

  bool CheckVelocityConvergence(const double NormDv, double &errorNormDv)
  {
    ModelPart &rModelPart = BaseType::GetModelPart();

    double NormV = 0.00;
    errorNormDv = 0;

#pragma omp parallel reduction(+ \
                               : NormV)
    {
      ModelPart::NodeIterator NodeBegin;
      ModelPart::NodeIterator NodeEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
      for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
      {
        const array_1d<double, 3> &Vel = itNode->FastGetSolutionStepValue(VELOCITY);

        double NormVelNode = 0;

        for (unsigned int d = 0; d < 3; ++d)
        {
          NormVelNode += Vel[d] * Vel[d];
          NormV += Vel[d] * Vel[d];
        }
      }
    }

    BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(NormV);

    NormV = sqrt(NormV);

    if (NormV == 0.0)
      NormV = 1.00;

    errorNormDv = NormDv / NormV;

    if (BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
    {
      std::cout << "The norm of velocity increment is: " << NormDv << std::endl;
      std::cout << "The norm of velocity is: " << NormV << std::endl;
      std::cout << "Velocity error: " << errorNormDv << "mVelocityTolerance: " << mVelocityTolerance << std::endl;
    }
    /* else{ */
    /*   std::cout<<"Velocity error: "<< errorNormDv <<" velTol: " << mVelocityTolerance<< std::endl; */
    /* } */

    if (errorNormDv < mVelocityTolerance)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  bool CheckPressureConvergence(const double NormDp, double &errorNormDp)
  {
    ModelPart &rModelPart = BaseType::GetModelPart();

    double NormP = 0.00;
    errorNormDp = 0;

#pragma omp parallel reduction(+ \
                               : NormP)
    {
      ModelPart::NodeIterator NodeBegin;
      ModelPart::NodeIterator NodeEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
      for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
      {
        const double Pr = itNode->FastGetSolutionStepValue(PRESSURE);
        NormP += Pr * Pr;
      }
    }

    BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(NormP);

    NormP = sqrt(NormP);

    if (NormP == 0.0)
      NormP = 1.00;

    errorNormDp = NormDp / NormP;

    if (BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
    {
      std::cout << "         The norm of pressure increment is: " << NormDp << std::endl;
      std::cout << "         The norm of pressure is: " << NormP << std::endl;
      std::cout << "         Pressure error: " << errorNormDp << std::endl;
    }
    /* else{ */
    /*     std::cout<<"         Pressure error: "<<errorNormDp <<" presTol: "<<mPressureTolerance << std::endl; */
    /* } */

    if (errorNormDp < mPressureTolerance)
    {
      return true;
    }
    else
      return false;
  }

  bool FixTimeStepMomentum(const double DvErrorNorm)
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    double currentTime = rCurrentProcessInfo[TIME];
    double timeInterval = rCurrentProcessInfo[DELTA_TIME];
    double minTolerance = 0.005;
    bool fixedTimeStep = false;
    if (currentTime < 10 * timeInterval)
    {
      minTolerance = 10;
    }

    if ((DvErrorNorm > minTolerance || (DvErrorNorm < 0 && DvErrorNorm > 0) || (DvErrorNorm != DvErrorNorm)) &&
        DvErrorNorm != 0 &&
        (DvErrorNorm != 1 || currentTime > timeInterval))
    {
      rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, true);
      std::cout << "NOT GOOD CONVERGENCE!!! I'll reduce the next time interval" << DvErrorNorm << std::endl;
      minTolerance = 0.05;
      if (DvErrorNorm > minTolerance)
      {
        std::cout << "BAD CONVERGENCE!!! I GO AHEAD WITH THE PREVIOUS VELOCITY AND PRESSURE FIELDS" << DvErrorNorm << std::endl;
        fixedTimeStep = true;
#pragma omp parallel
        {
          ModelPart::NodeIterator NodeBegin;
          ModelPart::NodeIterator NodeEnd;
          OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
          for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
          {
            itNode->FastGetSolutionStepValue(VELOCITY, 0) = itNode->FastGetSolutionStepValue(VELOCITY, 1);
            itNode->FastGetSolutionStepValue(PRESSURE, 0) = itNode->FastGetSolutionStepValue(PRESSURE, 1);
            itNode->FastGetSolutionStepValue(ACCELERATION, 0) = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
          }
        }
      }
    }
    else
    {
      rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, false);
    }
    return fixedTimeStep;
  }

  bool CheckMomentumConvergence(const double DvErrorNorm)
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    double currentTime = rCurrentProcessInfo[TIME];
    double timeInterval = rCurrentProcessInfo[DELTA_TIME];
    double minTolerance = 0.99999;
    bool fixedTimeStep = false;

    if ((DvErrorNorm > minTolerance || (DvErrorNorm < 0 && DvErrorNorm > 0) || (DvErrorNorm != DvErrorNorm)) &&
        DvErrorNorm != 0 &&
        (DvErrorNorm != 1 || currentTime > timeInterval))
    {
      rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, true);
      std::cout << "           BAD CONVERGENCE DETECTED DURING THE ITERATIVE LOOP!!! error: " << DvErrorNorm << " higher than 0.9999" << std::endl;
      std::cout << "      I GO AHEAD WITH THE PREVIOUS VELOCITY AND PRESSURE FIELDS" << std::endl;
      fixedTimeStep = true;
#pragma omp parallel
      {
        ModelPart::NodeIterator NodeBegin;
        ModelPart::NodeIterator NodeEnd;
        OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
        for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
        {
          itNode->FastGetSolutionStepValue(VELOCITY, 0) = itNode->FastGetSolutionStepValue(VELOCITY, 1);
          itNode->FastGetSolutionStepValue(PRESSURE, 0) = itNode->FastGetSolutionStepValue(PRESSURE, 1);
          itNode->FastGetSolutionStepValue(ACCELERATION, 0) = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
        }
      }
    }
    else
    {
      rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, false);
    }
    return fixedTimeStep;
  }

  bool FixTimeStepContinuity(const double DvErrorNorm)
  {
    ModelPart &rModelPart = BaseType::GetModelPart();
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    double currentTime = rCurrentProcessInfo[TIME];
    double timeInterval = rCurrentProcessInfo[DELTA_TIME];
    double minTolerance = 0.01;
    bool fixedTimeStep = false;
    if (currentTime < 10 * timeInterval)
    {
      minTolerance = 10;
    }

    if ((DvErrorNorm > minTolerance || (DvErrorNorm < 0 && DvErrorNorm > 0) || (DvErrorNorm != DvErrorNorm)) &&
        DvErrorNorm != 0 &&
        (DvErrorNorm != 1 || currentTime > timeInterval))
    {
      fixedTimeStep = true;
      rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, true);
      if (DvErrorNorm > 10 * minTolerance)
      {
        rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, true);
        std::cout << "           BAD CONVERGENCE DETECTED DURING THE ITERATIVE LOOP!!! error: " << DvErrorNorm << " higher than 0.9999" << std::endl;
        std::cout << "      I GO AHEAD WITH THE PREVIOUS VELOCITY AND PRESSURE FIELDS" << std::endl;
        fixedTimeStep = true;
#pragma omp parallel
        {
          ModelPart::NodeIterator NodeBegin;
          ModelPart::NodeIterator NodeEnd;
          OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
          for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
          {
            itNode->FastGetSolutionStepValue(VELOCITY, 0) = itNode->FastGetSolutionStepValue(VELOCITY, 1);
            itNode->FastGetSolutionStepValue(PRESSURE, 0) = itNode->FastGetSolutionStepValue(PRESSURE, 1);
            itNode->FastGetSolutionStepValue(ACCELERATION, 0) = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
          }
        }
      }
    }
    else
    {
      rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, false);
    }
    return fixedTimeStep;
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

  ///@name Static Member Variables
  ///@{

  ///@}
  ///@name Member Variables
  ///@{

  double mVelocityTolerance;

  double mPressureTolerance;

  unsigned int mMaxPressureIter;

  unsigned int mDomainSize;

  unsigned int mTimeOrder;

  bool mReformDofSet;

  // Fractional step index.
  /*  1 : Momentum step (calculate fractional step velocity)
      * 2-3 : Unused (reserved for componentwise calculation of frac step velocity)
      * 4 : Pressure step
      * 5 : Computation of projections
      * 6 : End of step velocity
      */
  //    unsigned int mStepId;

  /// Scheme for the solution of the momentum equation
  StrategyPointerType mpMomentumStrategy;

  /// Scheme for the solution of the mass equation
  StrategyPointerType mpPressureStrategy;

  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  virtual void InitializeStrategy(SolverSettingsType &rSolverConfig)
  {
    KRATOS_TRY;

    mTimeOrder = rSolverConfig.GetTimeOrder();

    // Check that input parameters are reasonable and sufficient.
    this->Check();

    //ModelPart& rModelPart = this->GetModelPart();

    mDomainSize = rSolverConfig.GetDomainSize();

    mReformDofSet = rSolverConfig.GetReformDofSet();

    BaseType::SetEchoLevel(rSolverConfig.GetEchoLevel());

    // Initialize strategies for each step
    bool HaveVelStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Velocity, mpMomentumStrategy);

    if (HaveVelStrategy)
    {
      rSolverConfig.FindTolerance(SolverSettingsType::Velocity, mVelocityTolerance);
      /* rSolverConfig.FindMaxIter(SolverSettingsType::Velocity,mMaxVelocityIter); */
    }
    else
    {
      KRATOS_THROW_ERROR(std::runtime_error, "TwoStepVPStrategy error: No Velocity strategy defined in FractionalStepSettings", "");
    }

    bool HavePressStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Pressure, mpPressureStrategy);

    if (HavePressStrategy)
    {
      rSolverConfig.FindTolerance(SolverSettingsType::Pressure, mPressureTolerance);
      rSolverConfig.FindMaxIter(SolverSettingsType::Pressure, mMaxPressureIter);
    }
    else
    {
      KRATOS_THROW_ERROR(std::runtime_error, "TwoStepVPStrategy error: No Pressure strategy defined in FractionalStepSettings", "");
    }

    // Check input parameters
    this->Check();

    KRATOS_CATCH("");
  }

  ///@}
  ///@name Private  Access
  ///@{

  ///@}
  ///@name Private Inquiry
  ///@{

  ///@}
  ///@name Un accessible methods
  ///@{

  /// Assignment operator.
  TwoStepVPStrategy &operator=(TwoStepVPStrategy const &rOther) {}

  /// Copy constructor.
  TwoStepVPStrategy(TwoStepVPStrategy const &rOther) {}

  ///@}

}; /// Class TwoStepVPStrategy

///@}
///@name Type Definitions
///@{

///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_V_P_STRATEGY_H
