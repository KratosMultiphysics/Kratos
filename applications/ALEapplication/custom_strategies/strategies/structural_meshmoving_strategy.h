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

    GenerateMeshPart();

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

    // Setting mesh to initial configuration
    for (ModelPart::NodeIterator i = (*mp_mesh_model_part).NodesBegin();
         i != (*mp_mesh_model_part).NodesEnd(); ++i) {

      (i)->X() = (i)->X0();
      (i)->Y() = (i)->Y0();
      (i)->Z() = (i)->Z0();
    }

    // Solve for the mesh movement
    m_strategy->Solve();

    // Update FEM-base
    CalculateMeshVelocities();
    MoveMesh();

    // Clearing the system if needed
    if (m_reform_dof_set_at_each_step == true)
      m_strategy->Clear();

    return 0.0;

    KRATOS_CATCH("");
  }

  void CalculateMeshVelocities() {
    KRATOS_TRY;

    double delta_time = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];

    KRATOS_ERROR_IF(delta_time <= 0.0)<< "Invalid DELTA_TIME." << std::endl;

    double coeff = 1 / delta_time;
    if (m_time_order == 1) // mesh velocity calculated as (x(n+1)-x(n))/Dt
    {
      for (ModelPart::NodeIterator i = (*mp_mesh_model_part).NodesBegin();
           i != (*mp_mesh_model_part).NodesEnd(); ++i) {

        array_1d<double, 3> &mesh_v =
            (i)->FastGetSolutionStepValue(MESH_VELOCITY);
        array_1d<double, 3> &disp =
            (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        array_1d<double, 3> &dispold =
            (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
        noalias(mesh_v) = disp - dispold;
        mesh_v *= coeff;
      }
    } else if (m_time_order ==
               2) // mesh velocity calculated as (3*x(n+1)-4*x(n)+x(n-1))/(2*Dt)
    {
      double c1 = 1.50 * coeff;
      double c2 = -2.0 * coeff;
      double c3 = 0.50 * coeff;

      for (ModelPart::NodeIterator i = (*mp_mesh_model_part).NodesBegin();
           i != (*mp_mesh_model_part).NodesEnd(); ++i) {

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

  void SetEchoLevel(int Level) override { m_strategy->SetEchoLevel(Level); }

  void MoveMesh() override {
    for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
         i != BaseType::GetModelPart().NodesEnd(); ++i) {
      (i)->X() = (i)->X0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_X);
      (i)->Y() = (i)->Y0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_Y);
      (i)->Z() = (i)->Z0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_Z);
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

  void GenerateMeshPart() {
    mp_mesh_model_part = ModelPart::Pointer(new ModelPart("MeshPart", 1));

    mp_mesh_model_part->Nodes() = BaseType::GetModelPart().Nodes();

    mp_mesh_model_part->GetNodalSolutionStepVariablesList() =
        BaseType::GetModelPart().GetNodalSolutionStepVariablesList();
    mp_mesh_model_part->SetBufferSize(BaseType::GetModelPart().GetBufferSize());

    // Creating mesh elements
    ModelPart::ElementsContainerType &MeshElems =
        mp_mesh_model_part->Elements();
    Element::Pointer pElem;

    for (ModelPart::ElementsContainerType::iterator it =
             BaseType::GetModelPart().ElementsBegin();
         it != BaseType::GetModelPart().ElementsEnd(); ++it) {

      pElem = Kratos::make_shared<StructuralMeshMovingElement>(
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
