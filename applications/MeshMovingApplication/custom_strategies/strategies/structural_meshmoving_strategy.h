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
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

#if !defined(KRATOS_NEW_STRUCTURAL_MESHMOVING_STRATEGY)
#define KRATOS_NEW_STRUCTURAL_MESHMOVING_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "custom_elements/structural_meshmoving_element.h"
#include "custom_utilities/move_mesh_utilities.h"
#include "containers/model.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/variable_utils.h"

#include "includes/mesh_moving_variables.h"

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
class StructuralMeshMovingStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
  /**@name Type Definitions */
  /*@{ */

  /** Counted pointer of ClassName */
  KRATOS_CLASS_POINTER_DEFINITION(StructuralMeshMovingStrategy);

  typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
  typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
  typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;

  /*@} */
  /**@name Life Cycle
   */
  /*@{ */

  /** Constructor.
   */
  StructuralMeshMovingStrategy(ModelPart &model_part,
                               typename TLinearSolver::Pointer pNewLinearSolver,
                               int TimeOrder = 2,
                               bool ReformDofSetAtEachStep = false,
                               bool ComputeReactions = false,
                               bool CalculateMeshVelocities = true,
                               int EchoLevel = 0,
                               const double PoissonRatio = 0.3,
                               const bool ReInitializeModelPartEachStep = false)
      : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part) {
    KRATOS_TRY

    mreform_dof_set_at_each_step = ReformDofSetAtEachStep || ReInitializeModelPartEachStep;
    mcompute_reactions = ComputeReactions;
    mcalculate_mesh_velocities = CalculateMeshVelocities;
    mecho_level = EchoLevel;
    mtime_order = TimeOrder;
    mreinitialize_model_part_at_each_step = ReInitializeModelPartEachStep;
    bool calculate_norm_dx_flag = false;

    typename SchemeType::Pointer pscheme = typename SchemeType::Pointer(
        new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,
                                                       TDenseSpace>());

    const std::string element_type = "StructuralMeshMovingElement";
    mpmesh_model_part = MoveMeshUtilities::GenerateMeshPart(
        BaseType::GetModelPart(), element_type);

    mpmesh_model_part->pGetProperties(0)->GetValue(MESH_POISSON_RATIO) = PoissonRatio;

    mpbulider_and_solver = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace,
                                               TLinearSolver>(
            pNewLinearSolver));

    mstrategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,TLinearSolver>(
        *mpmesh_model_part,
        pscheme,
        mpbulider_and_solver,
        mcompute_reactions,
        mreform_dof_set_at_each_step,
        calculate_norm_dx_flag));

    mstrategy->SetEchoLevel(mecho_level);

    KRATOS_CATCH("")
  }

  virtual ~StructuralMeshMovingStrategy()
  {
    Model& owner_model = mpmesh_model_part->GetModel();
    std::string name = mpmesh_model_part->Name();
    owner_model.DeleteModelPart(name);
  }

  void Initialize() override {}

  double Solve() override {
    KRATOS_TRY;

    if (mreinitialize_model_part_at_each_step) {
        MoveMeshUtilities::InitializeMeshPartWithElements(
            *mpmesh_model_part, BaseType::GetModelPart(),
            mpmesh_model_part->pGetProperties(0), "StructuralMeshMovingElement");
    }

    VariableUtils().UpdateCurrentToInitialConfiguration(
        mpmesh_model_part->GetCommunicator().LocalMesh().Nodes());

    // Solve for the mesh movement
    mstrategy->Solve();

    // Clearing the system if needed
    if (mreform_dof_set_at_each_step == true)
      mstrategy->Clear();

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
  ModelPart* mpmesh_model_part;

  typename BaseType::Pointer mstrategy;
  typename TBuilderAndSolverType::Pointer mpbulider_and_solver;

  int mecho_level;
  int mtime_order;
  bool mreform_dof_set_at_each_step;
  bool mcompute_reactions;
  bool mcalculate_mesh_velocities;
  bool mreinitialize_model_part_at_each_step;

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
  StructuralMeshMovingStrategy(const StructuralMeshMovingStrategy &Other);

  /*@} */

}; /* Class StructuralMeshMovingStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */
}
/* namespace Kratos.*/

#endif /* KRATOS_NEW_STRUCTURAL_MESHMOVING_STRATEGY  defined */
