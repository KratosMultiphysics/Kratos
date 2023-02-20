//
//   Project Name:        KratosPFEMFluidDynamicsApplication $
//   Last modified by:    $Author:                   AFranci $
//   Date:                $Date:                January 2023 $
//   Revision:            $Revision:                     0.0 $
//
//

#ifndef KRATOS_TWO_STEP_V_P_THERMAL_STRATEGY_H
#define KRATOS_TWO_STEP_V_P_THERMAL_STRATEGY_H

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

#include "two_step_v_p_thermal_strategy.h"

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
  class TwoStepVPThermalStrategy : public TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
  {
  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(TwoStepVPThermalStrategy);

    typedef TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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
    using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mVelocityTolerance;
    using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mPressureTolerance;
    using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mMaxPressureIter;
    using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mDomainSize;
    using TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mReformDofSet;
    ///@}
    ///@name Life Cycle
    ///@{

    TwoStepVPThermalStrategy(ModelPart &rModelPart,
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
    virtual ~TwoStepVPThermalStrategy() {}

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
      buffer << "TwoStepVPThermalStrategy";
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "TwoStepVPThermalStrategy";
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
#pragma omp parallel
      {
        ModelPart::NodeIterator NodeBegin;
        ModelPart::NodeIterator NodeEnd;
        OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
        for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
        {
          if (itNode->SolutionStepsDataHas(HEAT_FLUX))
          {
            itNode->FastGetSolutionStepValue(HEAT_FLUX) = 0;
          }
        }
      }
      if (it == 0)
      {
        mpMomentumStrategy->InitializeSolutionStep();
      }

      NormDv = mpMomentumStrategy->Solve();

      Parameters extrapolation_parameters(R"(
        {
        	"list_of_variables": ["MECHANICAL_DISSIPATION"]
        })");
      auto extrapolation_process = IntegrationValuesExtrapolationToNodesProcess(rModelPart, extrapolation_parameters);
      extrapolation_process.Execute();

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

    void SetEchoLevel(int Level) override
    {
      BaseType::SetEchoLevel(Level);
      int StrategyLevel = Level > 0 ? Level - 1 : 0;
      mpMomentumStrategy->SetEchoLevel(StrategyLevel);
      mpPressureStrategy->SetEchoLevel(StrategyLevel);
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

    // double mVelocityTolerance;

    // double mPressureTolerance;

    // unsigned int mMaxPressureIter;

    // unsigned int mDomainSize;

    // unsigned int mTimeOrder;

    // bool mReformDofSet;

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
    TwoStepVPThermalStrategy &operator=(TwoStepVPThermalStrategy const &rOther) {}

    /// Copy constructor.
    TwoStepVPThermalStrategy(TwoStepVPThermalStrategy const &rOther) {}

    ///@}

  }; /// Class TwoStepVPThermalStrategy

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}

  ///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_V_P_STRATEGY_H
