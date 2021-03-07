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

#if !defined(KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY)
#define KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/mesh_moving_variables.h"
#include "custom_elements/laplacian_meshmoving_element.h"
#include "custom_utilities/move_mesh_utilities.h"
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
class LaplacianMeshMovingStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {

public:
  /**@name Type Definitions */
  /*@{ */

  /** Counted pointer of ClassName */
  KRATOS_CLASS_POINTER_DEFINITION(LaplacianMeshMovingStrategy);

  typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
  typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
  typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;

  /*@} */
  /**@name Life Cycle
   */
  /*@{ */

  LaplacianMeshMovingStrategy(ModelPart &ModelPart,
                              typename TLinearSolver::Pointer pNewLinearSolver,
                              int TimeOrder = 1,
                              bool ReformDofSetAtEachStep = false,
                              bool ComputeReactions = false,
                              bool CalculateMeshVelocities = true,
                              int EchoLevel = 0,
                              const bool ReInitializeModelPartEachStep = false)
      : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(ModelPart) {

    KRATOS_TRY;

    mreform_dof_set_at_each_step = ReformDofSetAtEachStep || ReInitializeModelPartEachStep;
    mecho_level = EchoLevel;
    mcompute_reactions = ComputeReactions;
    mcalculate_mesh_velocities = CalculateMeshVelocities;
    mtime_order = TimeOrder;
    bool calculate_norm_dx_flag = false;
    mreinitialize_model_part_at_each_step = ReInitializeModelPartEachStep;

    typename SchemeType::Pointer pscheme = typename SchemeType::Pointer(
        new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,
                                                       TDenseSpace>());

    const std::string element_type = "LaplacianMeshMovingElement";
    mpmesh_model_part = MoveMeshUtilities::GenerateMeshPart(
        BaseType::GetModelPart(), element_type);

    typedef Variable<double> VarComponent;

    mpbuilder_and_solver_x = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedEliminationBuilderAndSolverComponentwise<
            TSparseSpace, TDenseSpace, TLinearSolver, VarComponent>(
            pNewLinearSolver, MESH_DISPLACEMENT_X));

    mpbuilder_and_solver_y = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedEliminationBuilderAndSolverComponentwise<
            TSparseSpace, TDenseSpace, TLinearSolver, VarComponent>(
            pNewLinearSolver, MESH_DISPLACEMENT_Y));

    mpbuilder_and_solver_z = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedEliminationBuilderAndSolverComponentwise<
            TSparseSpace, TDenseSpace, TLinearSolver, VarComponent>(
            pNewLinearSolver, MESH_DISPLACEMENT_Z));

    mstrategy_x = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,TLinearSolver>(
        *mpmesh_model_part,
        pscheme,
        mpbuilder_and_solver_x,
        ComputeReactions,
        mreform_dof_set_at_each_step,
        calculate_norm_dx_flag));

    mstrategy_x->SetEchoLevel(mecho_level);

    mstrategy_y = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,TLinearSolver>(
        *mpmesh_model_part,
        pscheme,
        mpbuilder_and_solver_y,
        ComputeReactions,
        mreform_dof_set_at_each_step,
        calculate_norm_dx_flag));

    mstrategy_y->SetEchoLevel(mecho_level);

    mstrategy_z = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,TLinearSolver>(
        *mpmesh_model_part,
        pscheme,
        mpbuilder_and_solver_z,
        ComputeReactions,
        mreform_dof_set_at_each_step,
        calculate_norm_dx_flag));

    mstrategy_z->SetEchoLevel(mecho_level);

    KRATOS_CATCH("")
  }

  virtual ~LaplacianMeshMovingStrategy(){

  };

  void Initialize() override {
    FindNodalNeighboursProcess find_nodal_neighbours_process(*mpmesh_model_part);
    find_nodal_neighbours_process.Execute();
  }

  double Solve() override {
    KRATOS_TRY;

    if (mreinitialize_model_part_at_each_step) {
        MoveMeshUtilities::InitializeMeshPartWithElements(
            *mpmesh_model_part, BaseType::GetModelPart(),
            mpmesh_model_part->pGetProperties(0), "LaplacianMeshMovingElement");
    }

    ProcessInfo &rCurrentProcessInfo = (mpmesh_model_part)->GetProcessInfo();

    VariableUtils().UpdateCurrentToInitialConfiguration(
        mpmesh_model_part->GetCommunicator().LocalMesh().Nodes());

    unsigned int dimension =
        BaseType::GetModelPart().GetProcessInfo()[DOMAIN_SIZE];

    if (dimension == 2) {
      // X DIRECTION
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] = 1;
      mstrategy_x->Solve();
      // Y DIRECTION
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] = 2;
      mstrategy_y->Solve();
    } else {
      // X DIRECTION
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] = 1;
      mstrategy_x->Solve();
      // Y DIRECTION
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] = 2;
      mstrategy_y->Solve();
      // Z DIRECTION
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] = 3;
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
  ModelPart* mpmesh_model_part;

  typename BaseType::Pointer mstrategy_x;
  typename BaseType::Pointer mstrategy_y;
  typename BaseType::Pointer mstrategy_z;
  typename TBuilderAndSolverType::Pointer mpbuilder_and_solver_x;
  typename TBuilderAndSolverType::Pointer mpbuilder_and_solver_y;
  typename TBuilderAndSolverType::Pointer mpbuilder_and_solver_z;

  bool mreform_dof_set_at_each_step;
  bool mcompute_reactions;
  int mtime_order;
  int mecho_level;
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
  /**@name Unaccessible methods */
  /*@{ */

  /** Copy constructor.
   */
  LaplacianMeshMovingStrategy(const LaplacianMeshMovingStrategy &Other);

  /*@} */

}; /* Class LaplacianMeshMovingStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */
}
/* namespace Kratos.*/

#endif /* KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY  defined */
