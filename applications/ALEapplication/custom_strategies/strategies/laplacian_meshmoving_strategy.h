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

#if !defined(KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY)
#define KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "ale_application.h"
#include "custom_elements/laplacian_meshmoving_element.h"
#include "includes/model_part.h"
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

  LaplacianMeshMovingStrategy(ModelPart &model_part,
                              typename TLinearSolver::Pointer pNewLinearSolver,
                              int TimeOrder = 1,
                              bool ReformDofSetAtEachStep = false,
                              bool ComputeReactions = false, int EchoLevel = 0)
      : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part) {

    KRATOS_TRY;

    m_reform_dof_set_at_each_step = ReformDofSetAtEachStep;
    m_echo_level = EchoLevel;
    m_compute_reactions = ComputeReactions;
    m_time_order = TimeOrder;
    bool calculate_norm_dx_flag = false;

    
    typename SchemeType::Pointer pscheme = typename SchemeType::Pointer(
        new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,
                                                       TDenseSpace>());

    GenerateMeshPart();

    typedef typename Kratos::VariableComponent<
        Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>>
        VarComponent;

    mp_builder_and_solver_x = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedEliminationBuilderAndSolverComponentwise<
            TSparseSpace, TDenseSpace, TLinearSolver, VarComponent>(
            pNewLinearSolver, MESH_DISPLACEMENT_X));

    mp_builder_and_solver_y = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedEliminationBuilderAndSolverComponentwise<
            TSparseSpace, TDenseSpace, TLinearSolver, VarComponent>(
            pNewLinearSolver, MESH_DISPLACEMENT_Y));

    mp_builder_and_solver_z = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedEliminationBuilderAndSolverComponentwise<
            TSparseSpace, TDenseSpace, TLinearSolver, VarComponent>(
            pNewLinearSolver, MESH_DISPLACEMENT_Z));

    m_strategy_x = typename BaseType::Pointer(
        new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,
                                        TLinearSolver>(
            *mp_mesh_model_part, pscheme, pNewLinearSolver,
            mp_builder_and_solver_x, ComputeReactions,
            m_reform_dof_set_at_each_step, calculate_norm_dx_flag));

    m_strategy_x->SetEchoLevel(m_echo_level);

    m_strategy_y = typename BaseType::Pointer(
        new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,
                                        TLinearSolver>(
            *mp_mesh_model_part, pscheme, pNewLinearSolver,
            mp_builder_and_solver_y, ComputeReactions,
            m_reform_dof_set_at_each_step, calculate_norm_dx_flag));

    m_strategy_y->SetEchoLevel(m_echo_level);

    m_strategy_z = typename BaseType::Pointer(
        new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,
                                        TLinearSolver>(
            *mp_mesh_model_part, pscheme, pNewLinearSolver,
            mp_builder_and_solver_z, ComputeReactions,
            m_reform_dof_set_at_each_step, calculate_norm_dx_flag));

    m_strategy_z->SetEchoLevel(m_echo_level);

    KRATOS_CATCH("")
  }

  virtual ~LaplacianMeshMovingStrategy(){

  };

  void Initialize() override {
    unsigned int n_average_elements = 10;
    unsigned int n_average_nodes = 10;
    FindNodalNeighboursProcess find_nodal_neighbours_process =
        FindNodalNeighboursProcess(*mp_mesh_model_part, n_average_elements,
                                   n_average_nodes);
    find_nodal_neighbours_process.Execute();
  }

  double Solve() override {
    KRATOS_TRY;

    ProcessInfo &rCurrentProcessInfo = (mp_mesh_model_part)->GetProcessInfo();

    // Setting mesh to initial configuration
    SetMeshToInitialConfiguration();

    unsigned int dimension =
        BaseType::GetModelPart().GetProcessInfo()[DOMAIN_SIZE];

    if (dimension == 2) {
      // X DIRECTION
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] = 1;
      m_strategy_x->Solve();
      // Y DIRECTION
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] = 2;
      m_strategy_y->Solve();
    } else {
      // X DIRECTION
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] = 1;
      m_strategy_x->Solve();
      // Y DIRECTION
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] = 2;
      m_strategy_y->Solve();
      // Z DIRECTION
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] = 3;
      m_strategy_z->Solve();
    }
    // Update FEM-base
    CalculateMeshVelocities();
    MoveMesh();

    if (m_reform_dof_set_at_each_step == true)
      m_strategy_x->Clear();
      m_strategy_y->Clear();
      m_strategy_z->Clear();

    return 0.0;

    KRATOS_CATCH("");
  }

  void CalculateMeshVelocities() {
    KRATOS_TRY;

    const double delta_time = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];

    KRATOS_ERROR_IF(delta_time <= 0.0)<< "Invalid DELTA_TIME." << std::endl;

    const double coeff = 1 / delta_time;

    if (m_time_order == 1) {
      for (ModelPart::NodeIterator i = (*mp_mesh_model_part).GetCommunicator().LocalMesh().NodesBegin();
           i != (*mp_mesh_model_part).GetCommunicator().LocalMesh().NodesEnd(); ++i) {

        array_1d<double, 3> &mesh_v =
            (i)->FastGetSolutionStepValue(MESH_VELOCITY);
        array_1d<double, 3> &disp =
            (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        array_1d<double, 3> &dispold =
            (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
        noalias(mesh_v) = disp - dispold;
        mesh_v *= coeff;
      }
    } else if (m_time_order == 2) {
      const double c1 = 1.50 * coeff;
      const double c2 = -2.0 * coeff;
      const double c3 = 0.50 * coeff;

      for (ModelPart::NodeIterator i = (*mp_mesh_model_part).GetCommunicator().LocalMesh().NodesBegin();
           i != (*mp_mesh_model_part).GetCommunicator().LocalMesh().NodesEnd(); ++i) {

        array_1d<double, 3> &mesh_v =
            (i)->FastGetSolutionStepValue(MESH_VELOCITY);
        noalias(mesh_v) = c1 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        noalias(mesh_v) +=
            c2 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
        noalias(mesh_v) +=
            c3 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 2);
      }
    } else {
      KRATOS_ERROR << "Wrong TimeOrder: Acceptable values are: 1 and 2"
                   << std::endl;
    }

    KRATOS_CATCH("");
  }

  void MoveMesh() override {
    for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
         i != BaseType::GetModelPart().NodesEnd(); ++i) {
      (i)->X() = (i)->X0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_X);
      (i)->Y() = (i)->Y0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_Y);
      (i)->Z() = (i)->Z0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_Z);
    }
  }

  void SetMeshToInitialConfiguration() {
    for (ModelPart::NodeIterator i = (*mp_mesh_model_part).NodesBegin();
         i != (*mp_mesh_model_part).NodesEnd(); ++i) {

      (i)->X() = (i)->X0();
      (i)->Y() = (i)->Y0();
      (i)->Z() = (i)->Z0();
    }
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

  typename BaseType::Pointer m_strategy_x;
  typename BaseType::Pointer m_strategy_y;
  typename BaseType::Pointer m_strategy_z;
  typename TBuilderAndSolverType::Pointer mp_builder_and_solver_x;
  typename TBuilderAndSolverType::Pointer mp_builder_and_solver_y;
  typename TBuilderAndSolverType::Pointer mp_builder_and_solver_z;

  bool m_reform_dof_set_at_each_step;
  bool m_compute_reactions;
  int m_time_order;
  int m_echo_level;
  /*@} */
  /**@name Private Operators*/
  /*@{ */

  /*@} */
  /**@name Private Operations*/
  /*@{ */

  void GenerateMeshPart() {
    mp_mesh_model_part = ModelPart::Pointer(new ModelPart("MeshPart", 1));

    // initializing mesh nodes
    mp_mesh_model_part->Nodes() = BaseType::GetModelPart().Nodes();

    // creating mesh elements
    ModelPart::ElementsContainerType &MeshElems =
        mp_mesh_model_part->Elements();
    Element::Pointer pElem;

    for (ModelPart::ElementsContainerType::iterator it =
             BaseType::GetModelPart().ElementsBegin();

         it != BaseType::GetModelPart().ElementsEnd(); ++it) {

      pElem = Kratos::make_shared<LaplacianMeshMovingElement>(
          (*it).Id(), (*it).pGetGeometry(), (*it).pGetProperties());
      MeshElems.push_back(pElem);
    }
  }

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
