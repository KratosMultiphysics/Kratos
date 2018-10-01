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
#include "mesh_moving_application.h"
#include "custom_elements/laplacian_meshmoving_element.h"
#include "custom_utilities/move_mesh_utilities.h"
#include "processes/find_nodal_neighbours_process.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/solving_strategy.h"

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
                              int EchoLevel = 0)
      : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(ModelPart) {

    KRATOS_TRY;

    mreform_dof_set_at_each_step = ReformDofSetAtEachStep;
    mecho_level = EchoLevel;
    mcompute_reactions = ComputeReactions;
    mcalculate_mesh_velocities = CalculateMeshVelocities;
    mtime_order = TimeOrder;
    bool calculate_norm_dx_flag = false;

    typename SchemeType::Pointer pscheme = typename SchemeType::Pointer(
        new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,
                                                       TDenseSpace>());

    const std::string element_type = "LaplacianMeshMovingElement";
    mpmesh_model_part = MoveMeshUtilities::GenerateMeshPart(
        BaseType::GetModelPart(), element_type);

    typedef typename Kratos::VariableComponent<
        Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>>
        VarComponent;

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

    mstrategy_x = typename BaseType::Pointer(
        new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,
                                        TLinearSolver>(
            *mpmesh_model_part, pscheme, pNewLinearSolver,
            mpbuilder_and_solver_x, ComputeReactions,
            mreform_dof_set_at_each_step, calculate_norm_dx_flag));

    mstrategy_x->SetEchoLevel(mecho_level);

    mstrategy_y = typename BaseType::Pointer(
        new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,
                                        TLinearSolver>(
            *mpmesh_model_part, pscheme, pNewLinearSolver,
            mpbuilder_and_solver_y, ComputeReactions,
            mreform_dof_set_at_each_step, calculate_norm_dx_flag));

    mstrategy_y->SetEchoLevel(mecho_level);

    mstrategy_z = typename BaseType::Pointer(
        new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,
                                        TLinearSolver>(
            *mpmesh_model_part, pscheme, pNewLinearSolver,
            mpbuilder_and_solver_z, ComputeReactions,
            mreform_dof_set_at_each_step, calculate_norm_dx_flag));

    mstrategy_z->SetEchoLevel(mecho_level);

    KRATOS_CATCH("")
  }

  virtual ~LaplacianMeshMovingStrategy(){

  };

  void Initialize() override {
    unsigned int n_average_elements = 10;
    unsigned int n_average_nodes = 10;
    FindNodalNeighboursProcess find_nodal_neighbours_process =
        FindNodalNeighboursProcess(*mpmesh_model_part, n_average_elements,
                                   n_average_nodes);
    find_nodal_neighbours_process.Execute();
  }

  double Solve() override {
    KRATOS_TRY;

    ProcessInfo &rCurrentProcessInfo = (mpmesh_model_part)->GetProcessInfo();

    MoveMeshUtilities::SetMeshToInitialConfiguration(
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
    // Update FEM-base
    const double delta_time =
        BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];

    if (mcalculate_mesh_velocities == true)
        MoveMeshUtilities::CalculateMeshVelocities(mpmesh_model_part.get(), mtime_order,
                                                   delta_time);
    MoveMeshUtilities::MoveMesh(
        mpmesh_model_part->GetCommunicator().LocalMesh().Nodes());

    if (mreform_dof_set_at_each_step == true)
    {
        mstrategy_x->Clear();
        mstrategy_y->Clear();
        mstrategy_z->Clear();
    }

    return 0.0;

    KRATOS_CATCH("");
  }

  void UpdateReferenceMesh()
  {
  MoveMeshUtilities::UpdateReferenceMesh(BaseType::GetModelPart());
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
  std::unique_ptr<ModelPart> mpmesh_model_part;

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
