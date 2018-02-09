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
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

#if !defined(KRATOS_NEW_STRUCTURAL_MESHMOVING_STRATEGY)
#define KRATOS_NEW_STRUCTURAL_MESHMOVING_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "custom_elements/structural_meshmoving_element.h"
#include "includes/model_part.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_utilities/move_mesh_utilities.h"

#include "ale_application.h"

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
                               bool ComputeReactions = false, int EchoLevel = 0)
      : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part) {
    KRATOS_TRY

    m_reform_dof_set_at_each_step = ReformDofSetAtEachStep;
    m_compute_reactions = ComputeReactions;
    m_echo_level = EchoLevel;
    m_time_order = TimeOrder;
    bool calculate_norm_dx_flag = false;

    typename SchemeType::Pointer pscheme = typename SchemeType::Pointer(
        new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,
                                                       TDenseSpace>());

    std::string element_type = "StructuralMeshMovingElement";
    mp_mesh_model_part = MoveMeshUtilities::GenerateMeshPart(BaseType::GetModelPart(), element_type);

    mp_bulider_and_solver = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace,
                                               TLinearSolver>(
            pNewLinearSolver));

    m_strategy = typename BaseType::Pointer(
        new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,
                                        TLinearSolver>(
            *mp_mesh_model_part, pscheme, pNewLinearSolver,
            mp_bulider_and_solver, m_compute_reactions,
            m_reform_dof_set_at_each_step, calculate_norm_dx_flag));

    m_strategy->SetEchoLevel(m_echo_level);

    KRATOS_CATCH("")
  }

  virtual ~StructuralMeshMovingStrategy() {}

  void Initialize() override {}

  double Solve() override {
    KRATOS_TRY;

    MoveMeshUtilities::SetMeshToInitialConfiguration(mp_mesh_model_part->GetCommunicator().LocalMesh().Nodes());

    // Solve for the mesh movement
    m_strategy->Solve();

    // Update FEM-base
    const double delta_time = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];
    MoveMeshUtilities::CalculateMeshVelocities(mp_mesh_model_part, m_time_order, delta_time);
    MoveMeshUtilities::MoveMesh(mp_mesh_model_part->GetCommunicator().LocalMesh().Nodes());

    // Clearing the system if needed
    if (m_reform_dof_set_at_each_step == true)
      m_strategy->Clear();

    return 0.0;

    KRATOS_CATCH("");
  }

  void UpdateReferenceMesh() {
    for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
         i != BaseType::GetModelPart().NodesEnd(); ++i) {
      (i)->X0() = (i)->X();
      (i)->Y0() = (i)->Y();
      (i)->Z0() = (i)->Z();
    }
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
  ModelPart::Pointer mp_mesh_model_part;

  typename BaseType::Pointer m_strategy;
  typename TBuilderAndSolverType::Pointer mp_bulider_and_solver;

  int m_echo_level;
  int m_time_order;
  bool m_reform_dof_set_at_each_step;
  bool m_compute_reactions;

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
