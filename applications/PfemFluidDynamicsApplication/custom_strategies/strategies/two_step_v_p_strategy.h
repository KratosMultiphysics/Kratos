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
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

#include "custom_utilities/solver_settings.h"

#include "custom_strategies/strategies/gauss_seidel_linear_strategy.h"

#include "pfem_fluid_dynamics_application_variables.h"

#include "processes/integration_values_extrapolation_to_nodes_process.h"

#include "v_p_strategy.h"

#include <stdio.h>
#include <cmath>

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
  class TwoStepVPStrategy : public VPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
  {
  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(TwoStepVPStrategy);

    typedef VPStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

    typedef TwoStepVPSolverSettings<TSparseSpace, TDenseSpace, TLinearSolver> SolverSettingsType;

    ///@}
    ///@name Life Cycle
    ///@{

    TwoStepVPStrategy(ModelPart &rModelPart,
                      typename TLinearSolver::Pointer pVelocityLinearSolver,
                      typename TLinearSolver::Pointer pPressureLinearSolver,
                      bool ReformDofSet = true,
                      double VelTol = 0.0001,
                      double PresTol = 0.0001,
                      int MaxPressureIterations = 1, // Only for predictor-corrector
                      unsigned int TimeOrder = 2,
                      unsigned int DomainSize = 2) : BaseType(rModelPart, pVelocityLinearSolver, pPressureLinearSolver, ReformDofSet, DomainSize),
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
      typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
      typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

      // initializing fractional velocity solution step
      typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;
      typename SchemeType::Pointer pScheme;

      typename SchemeType::Pointer Temp = typename SchemeType::Pointer(new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>());
      pScheme.swap(Temp);

      // CONSTRUCTION OF VELOCITY
      // BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pVelocityLinearSolver));
      BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pVelocityLinearSolver));

      this->mpMomentumStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pVelocityLinearSolver, vel_build, ReformDofAtEachIteration, CalculateNormDxFlag));

      vel_build->SetCalculateReactionsFlag(false);

      /* BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pPressureLinearSolver, PRESSURE)); */
      BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pPressureLinearSolver));

      this->mpPressureStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pPressureLinearSolver, pressure_build, ReformDofAtEachIteration, CalculateNormDxFlag));

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

      ModelPart &rModelPart = BaseType::GetModelPart();

      const auto &r_current_process_info = rModelPart.GetProcessInfo();
      for (const auto &r_element : rModelPart.Elements())
      {
        ierr = r_element.Check(r_current_process_info);
        if (ierr != 0)
        {
          break;
        }
      }

      return ierr;

      KRATOS_CATCH("");
    }

    void SetTimeCoefficients(ProcessInfo &rCurrentProcessInfo)
    {
      KRATOS_TRY;

      if (mTimeOrder == 2)
      {
        // calculate the BDF coefficients
        double Dt = rCurrentProcessInfo[DELTA_TIME];
        double OldDt = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

        double Rho = OldDt / Dt;
        double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

        Vector &BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
        BDFcoeffs.resize(3, false);

        BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho);        // coefficient for step n+1 (3/2Dt if Dt is constant)
        BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0); // coefficient for step n (-4/2Dt if Dt is constant)
        BDFcoeffs[2] = TimeCoeff;                                  // coefficient for step n-1 (1/2Dt if Dt is constant)
      }
      else if (mTimeOrder == 1)
      {
        double Dt = rCurrentProcessInfo[DELTA_TIME];
        double TimeCoeff = 1.0 / Dt;

        Vector &BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
        BDFcoeffs.resize(2, false);

        BDFcoeffs[0] = TimeCoeff;  // coefficient for step n+1 (1/Dt)
        BDFcoeffs[1] = -TimeCoeff; // coefficient for step n (-1/Dt)
      }

      KRATOS_CATCH("");
    }

    bool SolveSolutionStep() override
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      double currentTime = rCurrentProcessInfo[TIME];
      double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      bool timeIntervalChanged = rCurrentProcessInfo[TIME_INTERVAL_CHANGED];
      unsigned int stepsWithChangedDt = rCurrentProcessInfo[STEPS_WITH_CHANGED_DT];
      bool converged = false;

      unsigned int maxNonLinearIterations = mMaxPressureIter;

      KRATOS_INFO("\nSolution with two_step_vp_strategy at t=") << currentTime << "s" << std::endl;

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

      double pressureNorm = 0;
      double velocityNorm = 0;

      this->SetBlockedAndIsolatedFlags();

      for (unsigned int it = 0; it < maxNonLinearIterations; ++it)
      {
        momentumConverged = this->SolveMomentumIteration(it, maxNonLinearIterations, fixedTimeStep, velocityNorm);

        this->UpdateTopology(rModelPart, BaseType::GetEchoLevel());

        if (fixedTimeStep == false)
        {
          continuityConverged = this->SolveContinuityIteration(it, maxNonLinearIterations, pressureNorm);
        }
        if (it == maxNonLinearIterations - 1 || ((continuityConverged && momentumConverged) && it > 2))
        {

          // double tensilStressSign = 1.0;
          //  ComputeErrorL2Norm(tensilStressSign);
          this->UpdateStressStrain();
        }

        if ((continuityConverged && momentumConverged) && it > 2)
        {
          rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, false);
          rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, false);
          converged = true;

          KRATOS_INFO("TwoStepVPStrategy") << "V-P strategy converged in " << it + 1 << " iterations." << std::endl;

          break;
        }
        if (fixedTimeStep == true)
        {
          break;
        }
      }

      if (!continuityConverged && !momentumConverged && BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
        std::cout << "Convergence tolerance not reached." << std::endl;

      if (mReformDofSet)
        this->Clear();

      return converged;
    }

    void FinalizeSolutionStep() override
    {
    }

    void InitializeSolutionStep() override
    {
        auto& r_model_part = BaseType::GetModelPart();

        // InitializeSolutionStep of momentum strategy
        r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP, 1);
        mpMomentumStrategy->InitializeSolutionStep();

        // InitializeSolutionStep of continuity strategy
        r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP, 5);
        mpPressureStrategy->InitializeSolutionStep();
    }

    void UpdateStressStrain() override
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();

#pragma omp parallel
      {
        ModelPart::ElementIterator ElemBegin;
        ModelPart::ElementIterator ElemEnd;
        OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElemBegin, ElemEnd);
        for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
        {
          itElem->InitializeNonLinearIteration(rCurrentProcessInfo);
        }
      }

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
     * @param rCurrentProcessInfo ProcessInfo instance from the fluid ModelPart. Must contain DELTA_TIME variables.
     */

    bool SolveMomentumIteration(unsigned int it, unsigned int maxIt, bool &fixedTimeStep, double &velocityNorm) override
    {

      ModelPart &rModelPart = BaseType::GetModelPart();
      int Rank = rModelPart.GetCommunicator().MyPID();
      bool ConvergedMomentum = false;
      double NormDv = 0;
      fixedTimeStep = false;
      // build momentum system and solve for fractional step velocity increment
      rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP, 1);

      // if (it == 0)
      // {
      //   mpMomentumStrategy->InitializeSolutionStep();
      // }

      NormDv = mpMomentumStrategy->Solve();

      if (BaseType::GetEchoLevel() > 1 && Rank == 0)
        std::cout << "-------------- s o l v e d ! ------------------" << std::endl;

      if (it == 0)
      {
        velocityNorm = this->ComputeVelocityNorm();
      }

      double DvErrorNorm = NormDv / velocityNorm;
      unsigned int iterationForCheck = 2;
      // Check convergence
      if (it == maxIt - 1)
      {
        KRATOS_INFO("Iteration") << it << "  Final Velocity error: " << DvErrorNorm << std::endl;
        ConvergedMomentum = this->FixTimeStepMomentum(DvErrorNorm, fixedTimeStep);
      }
      else if (it > iterationForCheck)
      {
        KRATOS_INFO("Iteration") << it << "  Velocity error: " << DvErrorNorm << std::endl;
        ConvergedMomentum = this->CheckMomentumConvergence(DvErrorNorm, fixedTimeStep);
      }
      else
      {
        KRATOS_INFO("Iteration") << it << "  Velocity error: " << DvErrorNorm << std::endl;
      }

      if (!ConvergedMomentum && BaseType::GetEchoLevel() > 0 && Rank == 0)
        std::cout << "Momentum equations did not reach the convergence tolerance." << std::endl;

      return ConvergedMomentum;
    }

    bool SolveContinuityIteration(unsigned int it, unsigned int maxIt, double &NormP) override
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      int Rank = rModelPart.GetCommunicator().MyPID();
      bool ConvergedContinuity = false;
      bool fixedTimeStep = false;
      double NormDp = 0;

      // 2. Pressure solution
      rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP, 5);

    //   if (it == 0)
    //   {
    //     mpPressureStrategy->InitializeSolutionStep();
    //   }
      NormDp = mpPressureStrategy->Solve();

      if (BaseType::GetEchoLevel() > 0 && Rank == 0)
        std::cout << "The norm of pressure is: " << NormDp << std::endl;

      if (it == 0)
      {
        NormP = this->ComputePressureNorm();
      }

      double DpErrorNorm = NormDp / (NormP);

      // Check convergence
      if (it == (maxIt - 1))
      {
        KRATOS_INFO("Iteration") << it << "  Final Pressure error: " << DpErrorNorm << std::endl;
        ConvergedContinuity = this->FixTimeStepContinuity(DpErrorNorm, fixedTimeStep);
      }
      else
      {
        KRATOS_INFO("Iteration") << it << "  Pressure error: " << DpErrorNorm << std::endl;
        ConvergedContinuity = this->CheckContinuityConvergence(DpErrorNorm, fixedTimeStep);
      }

      if (!ConvergedContinuity && BaseType::GetEchoLevel() > 0 && Rank == 0)
        std::cout << "Continuity equation did not reach the convergence tolerance." << std::endl;

      return ConvergedContinuity;
    }

    bool CheckVelocityConvergence(const double NormDv, double &errorNormDv) override
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      const int n_nodes = rModelPart.NumberOfNodes();

      double NormV = 0.00;
      errorNormDv = 0;

#pragma omp parallel for reduction(+ \
                                   : NormV)
      for (int i_node = 0; i_node < n_nodes; ++i_node)
      {
        const auto it_node = rModelPart.NodesBegin() + i_node;
        const auto &r_vel = it_node->FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < 3; ++d)
        {
          NormV += r_vel[d] * r_vel[d];
        }
      }
      NormV = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(NormV);
      NormV = sqrt(NormV);

      const double zero_tol = 1.0e-12;
      errorNormDv = (NormV < zero_tol) ? NormDv : NormDv / NormV;

      if (BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
      {
        std::cout << "The norm of velocity increment is: " << NormDv << std::endl;
        std::cout << "The norm of velocity is: " << NormV << std::endl;
        std::cout << "Velocity error: " << errorNormDv << "mVelocityTolerance: " << mVelocityTolerance << std::endl;
      }

      if (errorNormDv < mVelocityTolerance)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    bool CheckPressureConvergence(const double NormDp, double &errorNormDp, double &NormP) override
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      const int n_nodes = rModelPart.NumberOfNodes();

      NormP = 0.00;
      errorNormDp = 0;
      double tmp_NormP = 0.0;

#pragma omp parallel for reduction(+ : tmp_NormP)
      for (int i_node = 0; i_node < n_nodes; ++i_node)
      {
        const auto it_node = rModelPart.NodesBegin() + i_node;
        const double Pr = it_node->FastGetSolutionStepValue(PRESSURE);
        tmp_NormP += Pr * Pr;
      }
      NormP = tmp_NormP;
      NormP = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(NormP);
      NormP = sqrt(NormP);

      const double zero_tol = 1.0e-12;
      errorNormDp = (NormP < zero_tol) ? NormDp : NormDp / NormP;

      if (BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
      {
        std::cout << "         The norm of pressure increment is: " << NormDp << std::endl;
        std::cout << "         The norm of pressure is: " << NormP << std::endl;
        std::cout << "         Pressure error: " << errorNormDp << std::endl;
      }

      if (errorNormDp < mPressureTolerance)
      {
        return true;
      }
      else
        return false;
    }

    bool FixTimeStepMomentum(const double DvErrorNorm, bool &fixedTimeStep) override
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      double currentTime = rCurrentProcessInfo[TIME];
      double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      double minTolerance = 0.005;
      bool converged = false;

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
        if (DvErrorNorm < mVelocityTolerance)
        {
          converged = true;
        }
      }
      return converged;
    }

    bool CheckMomentumConvergence(const double DvErrorNorm, bool &fixedTimeStep) override
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      double currentTime = rCurrentProcessInfo[TIME];
      double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      double minTolerance = 0.99999;
      bool converged = false;

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
        if (DvErrorNorm < mVelocityTolerance)
        {
          converged = true;
        }
      }
      return converged;
    }

    bool FixTimeStepContinuity(const double DvErrorNorm, bool &fixedTimeStep) override
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      double currentTime = rCurrentProcessInfo[TIME];
      double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      double minTolerance = 0.01;
      bool converged = false;
      if (currentTime < 10 * timeInterval)
      {
        minTolerance = 10;
      }

      if ((DvErrorNorm > minTolerance || (DvErrorNorm < 0 && DvErrorNorm > 0) || (DvErrorNorm != DvErrorNorm)) &&
          DvErrorNorm != 0 &&
          (DvErrorNorm != 1 || currentTime > timeInterval))
      {
        fixedTimeStep = true;
        // rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, true);
        if (DvErrorNorm > 0.9999)
        {
          rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, true);
          std::cout << "           BAD PRESSURE CONVERGENCE DETECTED DURING THE ITERATIVE LOOP!!! error: " << DvErrorNorm << " higher than 0.1" << std::endl;
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
      else if (DvErrorNorm < mPressureTolerance)
      {
        converged = true;
        fixedTimeStep = false;
      }
      rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, false);
      return converged;
    }

    bool CheckContinuityConvergence(const double DvErrorNorm, bool &fixedTimeStep) override
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      bool converged = false;

      if (DvErrorNorm < mPressureTolerance)
      {
        converged = true;
        fixedTimeStep = false;
      }
      rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, false);
      return converged;
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
