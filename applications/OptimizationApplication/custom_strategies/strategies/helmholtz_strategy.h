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

#if !defined(KRATOS_HELMHOLTZ_VEC_STRATEGY)
#define KRATOS_HELMHOLTZ_VEC_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "containers/model.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/variable_utils.h"
#include "utilities/condition_number_utility.h"

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
  typedef typename BaseType::TSystemMatrixType                          TSystemMatrixType;
  typedef typename BaseType::TSystemVectorType                          TSystemVectorType;  
  typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;

  /*@} */
  /**@name Life Cycle
   */
  /*@{ */

  /** Constructor.
   */
  HelmholtzStrategy(ModelPart &model_part,
                               typename TLinearSolver::Pointer pNewLinearSolver,
                               bool ReformDofSetAtEachStep = false,
                               bool ComputeReactions = false,
                               int EchoLevel = 0,
                               const double PoissonRatio = 0.3)
      : ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part) {
    KRATOS_TRY

    mreform_dof_set_at_each_step = ReformDofSetAtEachStep;
    mcompute_reactions = ComputeReactions;
    mecho_level = EchoLevel;
    bool calculate_norm_dx_flag = false;

    mpscheme = typename SchemeType::Pointer(
        new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,
                                                       TDenseSpace>());

    mpbulider_and_solver = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace,
                                               TLinearSolver>(
            pNewLinearSolver));

    mpstrategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,TLinearSolver>(
        BaseType::GetModelPart(),
        mpscheme,
        mpbulider_and_solver,
        mcompute_reactions,
        mreform_dof_set_at_each_step,
        calculate_norm_dx_flag));

    mpstrategy->SetEchoLevel(mecho_level);

    KRATOS_CATCH("")
  }

  virtual ~HelmholtzStrategy()
  {

  }

  void Initialize() override {}

  double Solve() override {
    KRATOS_TRY;

    VariableUtils().UpdateCurrentToInitialConfiguration(
        BaseType::GetModelPart().GetCommunicator().LocalMesh().Nodes());

    // Solve for the mesh movement
    mpstrategy->Solve();

    // Clearing the system if needed
    if (mreform_dof_set_at_each_step == true)
      mpstrategy->Clear();

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
    TSystemMatrixType &GetSystemMatrix() override
    {
        return mpstrategy->GetSystemMatrix();
    }

    TSystemVectorType &GetSystemVector() override
    {
        return mpstrategy->GetSystemVector();
    }

    TSystemVectorType &GetSolutionVector() override
    {
        return mpstrategy->GetSolutionVector();
    }   

    void ExportSystem() 
    { 
        mpstrategy->Clear();
        mpstrategy->Initialize();
        mpstrategy->InitializeSolutionStep();
        mpstrategy->Predict();         
        mpbulider_and_solver->Build(mpscheme,BaseType::GetModelPart(),mpstrategy->GetSystemMatrix(),mpstrategy->GetSystemVector());

        std::stringstream matrix_market_name;
        matrix_market_name << "A_wo_D_BC.mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), mpstrategy->GetSystemMatrix(), false);

        std::stringstream matrix_market_vectname;
        matrix_market_vectname << "b_wo_D_BC.mm";
        TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), mpstrategy->GetSystemVector()); 

        mpstrategy->SolveSolutionStep();

        matrix_market_name.str("");
        matrix_market_name << "A_wi_D_BC.mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), mpstrategy->GetSystemMatrix(), false);

        matrix_market_vectname.str("");
        matrix_market_vectname << "b_wi_D_BC.mm";
        TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), mpstrategy->GetSystemVector()); 

        mpstrategy->FinalizeSolutionStep();

        mpstrategy->Clear();

    }       

    typename BaseType::Pointer GetStrategy() 
    {
        return mpstrategy;
    }     

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

  typename SchemeType::Pointer mpscheme;
  typename BaseType::Pointer mpstrategy;
  typename TBuilderAndSolverType::Pointer mpbulider_and_solver;

  int mecho_level;
  bool mreform_dof_set_at_each_step;
  bool mcompute_reactions;

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
  /**@name Un accessible methods */
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

#endif /* KRATOS_HELMHOLTZ_VEC_STRATEGY  defined */
