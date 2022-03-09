//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Reza Najian Asl
//

#if !defined(KRATOS_HELMHOLTZ_STRATEGY)
#define KRATOS_HELMHOLTZ_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "shape_optimization_application_variables.h"
#include "processes/find_nodal_neighbours_process.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
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

/// Short class definition.
/**   Detail class definition.
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class HelmholtzStrategy
    : public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {

public:
  /**@name Type Definitions */
  /*@{ */

  /** Counted pointer of ClassName */
  KRATOS_CLASS_POINTER_DEFINITION(HelmholtzStrategy);

  typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
  typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
  typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;

  /*@} */
  /**@name Life Cycle
   */
  /*@{ */

  HelmholtzStrategy(ModelPart &ModelPart,
                              typename TLinearSolver::Pointer pNewLinearSolver,
                              bool ReformDofSetAtEachStep = false,
                              bool ComputeReactions = false,
                              int EchoLevel = 0)
      : ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(ModelPart) {

    KRATOS_TRY;

    mreform_dof_set_at_each_step = ReformDofSetAtEachStep;
    mecho_level = EchoLevel;
    mcompute_reactions = ComputeReactions;
    bool calculate_norm_dx_flag = false;

    typename SchemeType::Pointer pscheme = typename SchemeType::Pointer(
        new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,
                                                       TDenseSpace>());


    typedef Variable<double> VarComponent;

    mpbuilder_and_solver_x = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedEliminationBuilderAndSolverComponentwise<
            TSparseSpace, TDenseSpace, TLinearSolver, VarComponent>(
            pNewLinearSolver, HELMHOLTZ_VARS_X));

    mpbuilder_and_solver_y = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedEliminationBuilderAndSolverComponentwise<
            TSparseSpace, TDenseSpace, TLinearSolver, VarComponent>(
            pNewLinearSolver, HELMHOLTZ_VARS_Y));

    mpbuilder_and_solver_z = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedEliminationBuilderAndSolverComponentwise<
            TSparseSpace, TDenseSpace, TLinearSolver, VarComponent>(
            pNewLinearSolver, HELMHOLTZ_VARS_Z));

    mstrategy_x = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,TLinearSolver>(
        BaseType::GetModelPart(),
        pscheme,
        mpbuilder_and_solver_x,
        ComputeReactions,
        mreform_dof_set_at_each_step,
        calculate_norm_dx_flag));

    mstrategy_x->SetEchoLevel(mecho_level);

    mstrategy_y = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,TLinearSolver>(
        BaseType::GetModelPart(),
        pscheme,
        mpbuilder_and_solver_y,
        ComputeReactions,
        mreform_dof_set_at_each_step,
        calculate_norm_dx_flag));

    mstrategy_y->SetEchoLevel(mecho_level);

    mstrategy_z = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,TLinearSolver>(
        BaseType::GetModelPart(),
        pscheme,
        mpbuilder_and_solver_z,
        ComputeReactions,
        mreform_dof_set_at_each_step,
        calculate_norm_dx_flag));

    mstrategy_z->SetEchoLevel(mecho_level);

    KRATOS_CATCH("")
  }

  virtual ~HelmholtzStrategy(){

  };

  void Initialize() override {
    FindNodalNeighboursProcess find_nodal_neighbours_process(BaseType::GetModelPart());
    find_nodal_neighbours_process.Execute();
  }

  double Solve() override {
    KRATOS_TRY;

    ProcessInfo &rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

    VariableUtils().UpdateCurrentToInitialConfiguration(
        BaseType::GetModelPart().GetCommunicator().LocalMesh().Nodes());

    unsigned int dimension =
        BaseType::GetModelPart().GetProcessInfo()[DOMAIN_SIZE];

    if (dimension == 2) {
      // X DIRECTION
      rCurrentProcessInfo[HELMHOLTZ_DIRECTION] = 1;
      mstrategy_x->Solve();
      // Y DIRECTION
      rCurrentProcessInfo[HELMHOLTZ_DIRECTION] = 2;
      mstrategy_y->Solve();
    } else {
      // X DIRECTION
      rCurrentProcessInfo[HELMHOLTZ_DIRECTION] = 1;
      mstrategy_x->Solve();
      // Y DIRECTION
      rCurrentProcessInfo[HELMHOLTZ_DIRECTION] = 2;
      mstrategy_y->Solve();
      // Z DIRECTION
      rCurrentProcessInfo[HELMHOLTZ_DIRECTION] = 3;
      mstrategy_z->Solve();
    }

    if (mreform_dof_set_at_each_step == true)
    {
        mstrategy_x->Clear();
        mstrategy_y->Clear();
        mstrategy_z->Clear();
    }

    return 0.0;

    KRATOS_CATCH("");
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

protected:
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

private:
  /**@name Static Member Variables */
  /*@{ */

  /*@} */
  /**@name Member Variables */
  /*@{ */

  typename BaseType::Pointer mstrategy_x;
  typename BaseType::Pointer mstrategy_y;
  typename BaseType::Pointer mstrategy_z;
  typename TBuilderAndSolverType::Pointer mpbuilder_and_solver_x;
  typename TBuilderAndSolverType::Pointer mpbuilder_and_solver_y;
  typename TBuilderAndSolverType::Pointer mpbuilder_and_solver_z;

  bool mreform_dof_set_at_each_step;
  bool mcompute_reactions;
  int mecho_level;

  /*@} */
  /**@name Private Operators*/
  /*@{ */

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
  /**@name Unaccessible methods */
  /*@{ */

  /** Copy constructor.
   */
  HelmholtzStrategy(const HelmholtzStrategy &Other);

  /*@} */

}; /* Class HelmholtzStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */
}
/* namespace Kratos.*/

#endif /* KRATOS_HELMHOLTZ_STRATEGY  defined */
