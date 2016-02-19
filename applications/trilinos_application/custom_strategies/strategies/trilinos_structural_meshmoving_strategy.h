/* ****************************************************************************
 *  Projectname:         $KratosATrilinosApplication
 *  Last Modified by:    $Author: Andreas.Mini@tum.de $
 *  Date:                $Date: November 2015 $
 *  Revision:            $Revision: 1.4 $
 * ***************************************************************************/

#if !defined(KRATOS_TRILINOS_STRUCTURAL_MESHMOVING_STRATEGY )
#define  KRATOS_TRILINOS_STRUCTURAL_MESHMOVING_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "../ALEapplication/custom_elements/structural_meshmoving_element.h"

/* Trilinos includes */
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
#include "custom_utilities/parallel_fill_communicator.h"

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
template<class TSparseSpace, class TDenseSpace,  //= DenseSpace<double>,
    class TLinearSolver  //= LinearSolver<TSparseSpace,TDenseSpace>
>
class TrilinosStructuralMeshMovingStrategy : public SolvingStrategy<
    TSparseSpace, TDenseSpace, TLinearSolver> {
 public:
  /**@name Type Definitions */
  /*@{ */

  /** Counted pointer of ClassName */
  KRATOS_CLASS_POINTER_DEFINITION( TrilinosStructuralMeshMovingStrategy );

  typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

  /*@} */
  /**@name Life Cycle
   */
  /*@{ */

  /** Constructor.
   */
  TrilinosStructuralMeshMovingStrategy(
      Epetra_MpiComm& Comm,
      ModelPart& model_part,
      typename TLinearSolver::Pointer pNewLinearSolver,
      int velocity_order = 1,
      bool reform_dof_at_every_step = false
  )
  : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part)
  {
    KRATOS_TRY

    // Passed variables
    mvel_order = velocity_order;
    mreform_dof_at_every_step = reform_dof_at_every_step;

    // Definitions for trilinos
    int guess_row_size;
    guess_row_size = 15;

    bool CalculateReactions = false;
    bool ReformDofAtEachIteration = false;
    bool CalculateNormDxFlag = false;

    //Generating Mesh Part
    GenerateMeshPart();

    typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
    typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
    ( new TrilinosResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace >() );

    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;
    BuilderSolverTypePointer builderSolver = BuilderSolverTypePointer
    (new TrilinosBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(Comm, guess_row_size, pNewLinearSolver) );

    mstrategy = typename BaseType::Pointer
    ( new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver >(*mpMeshModelPart,pscheme,pNewLinearSolver,builderSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag) );

    KRATOS_CATCH("")
  }

  /** Destructor.
   */
  virtual ~TrilinosStructuralMeshMovingStrategy() {}

  /** Destructor.
   */

  //*********************************************************************************
  //*********************************************************************************
  double Solve()
  {
    KRATOS_TRY;

    // Setting mesh to initial configuration
    for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin(); i != (*mpMeshModelPart).NodesEnd(); ++i)
    {
      (i)->X() = (i)->X0();
      (i)->Y() = (i)->Y0();
      (i)->Z() = (i)->Z0();
    }

    // Solve for mesh movement
    mstrategy->Solve();

    // Update FEM database
    CalculateMeshVelocities();
    BaseType::MoveMesh();

    // Clearing the system if needed
    if(mreform_dof_at_every_step == true)
    mstrategy->Clear();

    return 0.0;

    KRATOS_CATCH("")
  }

  //*********************************************************************************
  //*********************************************************************************
  void CalculateMeshVelocities()
  {
    KRATOS_TRY;

    double DeltaTime = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];

    if (DeltaTime <= 0.0)
    KRATOS_THROW_ERROR(std::logic_error, "Invalid DELTA_TIME.","");

    double coeff = 1/DeltaTime;
    if( mvel_order == 1)  //mesh velocity calculated as (x(n+1)-x(n))/Dt
    {
      for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin();
          i != (*mpMeshModelPart).NodesEnd(); ++i)
      {
        array_1d<double,3>& mesh_v = (i)->FastGetSolutionStepValue(MESH_VELOCITY);
        array_1d<double,3>& disp = (i)->FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double,3>& dispold = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
        noalias(mesh_v) = disp - dispold;
        mesh_v *= coeff;
      }
    }
    else  //mesh velocity calculated as (3*x(n+1)-4*x(n)+x(n-1))/(2*Dt)
    {
      double c1 = 1.50*coeff;
      double c2 = -2.0*coeff;
      double c3 = 0.50*coeff;

      for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin();
          i != (*mpMeshModelPart).NodesEnd(); ++i)
      {
        array_1d<double,3>& mesh_v = (i)->FastGetSolutionStepValue(MESH_VELOCITY);
        noalias(mesh_v) = c1 * (i)->FastGetSolutionStepValue(DISPLACEMENT);
        noalias(mesh_v) += c2 * (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
        noalias(mesh_v) += c3 * (i)->FastGetSolutionStepValue(DISPLACEMENT,2);
      }
    }

    KRATOS_CATCH("")
  }

  virtual void SetEchoLevel(int Level)
  {
    mstrategy->SetEchoLevel(Level);
  }

  void MoveNodes()
  {
    CalculateMeshVelocities();
    BaseType::MoveMesh();
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
  ModelPart::Pointer mpMeshModelPart;

  typename BaseType::Pointer mstrategy;

  int mvel_order;
  bool mreform_dof_at_every_step;

  /*@} */
  /**@name Private Operators*/
  /*@{ */

  /*@} */
  /**@name Private Operations*/
  /*@{ */

  void GenerateMeshPart()
  {
    // Initialize auxiliary model part storing the mesh elements
    mpMeshModelPart = ModelPart::Pointer( new ModelPart("MeshPart",1) );

    // Initializing mesh nodes
    mpMeshModelPart->Nodes() = BaseType::GetModelPart().Nodes();

    // Setup new auxiliary model part
    mpMeshModelPart->GetNodalSolutionStepVariablesList() = BaseType::GetModelPart().GetNodalSolutionStepVariablesList();
    mpMeshModelPart->SetBufferSize(BaseType::GetModelPart().GetBufferSize());

    // Create a communicator for the new model part and copy the partition information about nodes.
    typename Communicator::Pointer pMeshMPIComm = typename Communicator::Pointer( new MPICommunicator( &(BaseType::GetModelPart().GetNodalSolutionStepVariablesList()) ) );

    mpMeshModelPart->SetCommunicator( pMeshMPIComm );

    // Creating mesh elements
    ModelPart::ElementsContainerType& MeshElems = mpMeshModelPart->Elements();
    Element::Pointer pElem;

    for(ModelPart::ElementsContainerType::iterator it = BaseType::GetModelPart().ElementsBegin();
        it != BaseType::GetModelPart().ElementsEnd(); ++it)
    {
      pElem = Element::Pointer(new StructuralMeshMovingElement(
              (*it).Id(),
              (*it).pGetGeometry(),
              (*it).pGetProperties() ) );
      MeshElems.push_back(pElem);
    }

    // Optimize communicaton plan
    ParallelFillCommunicator CommunicatorGeneration( *mpMeshModelPart );
    CommunicatorGeneration.Execute();
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
  TrilinosStructuralMeshMovingStrategy(const TrilinosStructuralMeshMovingStrategy& Other);

  /*@} */

}; /* Class TrilinosStructuralMeshMovingStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */

}
/* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_MESHMOVING_STRATEGY  defined */

