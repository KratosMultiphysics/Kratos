/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: Andreas.Mini@tum.de $
 *  Date:                $Date: November 2015 $
 *  Revision:            $Revision: 1.4 $
 * ***************************************************************************/

#if !defined(KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY )
#define  KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "custom_elements/laplacian_meshmoving_element.h"

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
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class LaplacianMeshMovingStrategy : public SolvingStrategy<TSparseSpace,
    TDenseSpace, TLinearSolver> {

 public:
  /**@name Type Definitions */
  /*@{ */

  /** Counted pointer of ClassName */
  KRATOS_CLASS_POINTER_DEFINITION( LaplacianMeshMovingStrategy );

  typedef SolvingStrategy< TSparseSpace,TDenseSpace,TLinearSolver > BaseType;

  /*@} */
  /**@name Life Cycle
   */
  /*@{ */

  /** Constructor.
   */
  LaplacianMeshMovingStrategy(
      ModelPart& model_part,
      typename TLinearSolver::Pointer pNewLinearSolver,
      int dimension = 3,
      int velocity_order = 1,
      bool reform_dof_at_every_step = true
  )
  :SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part) {

    KRATOS_TRY

    GenerateMeshPart(dimension);

    mdimension = dimension;
    mvel_order = velocity_order;
    mreform_dof_at_every_step = reform_dof_at_every_step;

    typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;
    typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
    ( new ResidualBasedIncrementalUpdateStaticScheme
        <TSparseSpace, TDenseSpace>());

    bool CalculateReactions = false;
    bool ReformDofAtEachIteration = true;
    bool CalculateNormDxFlag = false;

    typedef typename BuilderAndSolver
    <TSparseSpace,TDenseSpace,TLinearSolver>
    ::Pointer BuilderSolverTypePointer;

    typedef typename Kratos::VariableComponent
    <Kratos::VectorComponentAdaptor<Kratos
    ::array_1d<double, 3> > > VarComponent;

    BuilderSolverTypePointer dir_x_build = BuilderSolverTypePointer
    (new ResidualBasedEliminationBuilderAndSolverComponentwise
        <TSparseSpace,TDenseSpace,TLinearSolver, VarComponent>
        ( pNewLinearSolver,DISPLACEMENT_X));

    BuilderSolverTypePointer dir_y_build = BuilderSolverTypePointer
    (new ResidualBasedEliminationBuilderAndSolverComponentwise
        <TSparseSpace,TDenseSpace,TLinearSolver, VarComponent>
        (pNewLinearSolver,DISPLACEMENT_Y));

    BuilderSolverTypePointer dir_z_build = BuilderSolverTypePointer
    (new ResidualBasedEliminationBuilderAndSolverComponentwise
        <TSparseSpace,TDenseSpace,TLinearSolver, VarComponent>
        (pNewLinearSolver,DISPLACEMENT_Z));

    mstrategy_x = typename BaseType::Pointer
    (new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
        (*mpMeshModelPart,pscheme,pNewLinearSolver,dir_x_build,
            CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag));

    mstrategy_y = typename BaseType::Pointer
    (new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
        (*mpMeshModelPart,pscheme,pNewLinearSolver,dir_y_build,
            CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag));

    mstrategy_z = typename BaseType::Pointer
    (new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
        (*mpMeshModelPart,pscheme,pNewLinearSolver,dir_z_build,
            CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag));

    KRATOS_CATCH("")
  }

  /** Destructor.
   */
  virtual ~LaplacianMeshMovingStrategy() {

  }

  //***************************************************************************
  //***************************************************************************
  double Solve()
  {
    KRATOS_TRY;

    ProcessInfo& rCurrentProcessInfo = (mpMeshModelPart)->GetProcessInfo();

    // Setting mesh to initial configuration
    for(ModelPart::NodeIterator i = (*mpMeshModelPart ).NodesBegin();
        i != (*mpMeshModelPart ).NodesEnd(); ++i) {

      (i)->X() = (i)->X0();
      (i)->Y() = (i)->Y0();
      (i)->Z() = (i)->Z0();
    }

    //X DIRECTION
    rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
    mstrategy_x->Solve();

    //Y DIRECTION
    if(mdimension > 1) {
      rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
      mstrategy_y->Solve();
    }

    //Z DIRECTION
    if(mdimension > 2) {
      rCurrentProcessInfo[FRACTIONAL_STEP] = 3;
      mstrategy_z->Solve();
    }

    MoveNodes();

    //clearing the system if needed
    if(mreform_dof_at_every_step == true) {
      mstrategy_x->Clear();
      mstrategy_y->Clear();
      mstrategy_z->Clear();
    }

    return 0.0;

    KRATOS_CATCH("");
  }

  //***************************************************************************
  //***************************************************************************
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
          i != (*mpMeshModelPart).NodesEnd(); ++i) {

        array_1d<double,3>& mesh_v = (i)
        ->FastGetSolutionStepValue(MESH_VELOCITY);
        array_1d<double,3>& disp = (i)
        ->FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double,3>& dispold = (i)
        ->FastGetSolutionStepValue(DISPLACEMENT,1);
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
          i != (*mpMeshModelPart).NodesEnd(); ++i) {

        array_1d<double,3>& mesh_v = (i)
        ->FastGetSolutionStepValue(MESH_VELOCITY);
        noalias(mesh_v) = c1 * (i)
        ->FastGetSolutionStepValue(DISPLACEMENT);
        noalias(mesh_v) += c2 * (i)
        ->FastGetSolutionStepValue(DISPLACEMENT,1);
        noalias(mesh_v) += c3 * (i)
        ->FastGetSolutionStepValue(DISPLACEMENT,2);
      }
    }

    KRATOS_CATCH("");
  }

  //***************************************************************************
  //***************************************************************************
  virtual void SetEchoLevel(int Level)
  {
    mstrategy_x->SetEchoLevel(Level);
    mstrategy_y->SetEchoLevel(Level);
    mstrategy_z->SetEchoLevel(Level);
  }

  //***************************************************************************
  //***************************************************************************
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

  typename BaseType::Pointer mstrategy_x;
  typename BaseType::Pointer mstrategy_y;
  typename BaseType::Pointer mstrategy_z;

  int mdimension;
  int mvel_order;
  bool mreform_dof_at_every_step;

  /*@} */
  /**@name Private Operators*/
  /*@{ */

  /*@} */
  /**@name Private Operations*/
  /*@{ */

  //***************************************************************************
  //***************************************************************************
  void GenerateMeshPart(int dimension)
  {
    mpMeshModelPart = ModelPart::Pointer( new ModelPart("MeshPart",1) );

    //initializing mesh nodes
    mpMeshModelPart->Nodes() = BaseType::GetModelPart().Nodes();

    //creating mesh elements
    ModelPart::ElementsContainerType& MeshElems = mpMeshModelPart->Elements();
    Element::Pointer pElem;

    if(dimension == 2)
    for(ModelPart::ElementsContainerType::iterator it
        = BaseType::GetModelPart().ElementsBegin();

        it != BaseType::GetModelPart().ElementsEnd(); ++it) {

      pElem = Element::Pointer(new LaplacianMeshMovingElement<2>(
              (*it).Id(),
              (*it).pGetGeometry(),
              (*it).pGetProperties() ) );
      MeshElems.push_back(pElem);
    }

    if(dimension == 3)
    for(ModelPart::ElementsContainerType::iterator it
        = BaseType::GetModelPart().ElementsBegin();

        it != BaseType::GetModelPart().ElementsEnd(); ++it) {

      pElem = Element::Pointer(new LaplacianMeshMovingElement<3>(
              (*it).Id(),
              (*it).pGetGeometry(),
              (*it).pGetProperties() ) );
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
  LaplacianMeshMovingStrategy(const LaplacianMeshMovingStrategy& Other);

  /*@} */

}; /* Class LaplacianMeshMovingStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */

}
/* namespace Kratos.*/

#endif /* KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY  defined */

