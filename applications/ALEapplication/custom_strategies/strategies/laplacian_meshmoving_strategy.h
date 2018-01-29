// ==============================================================================
/*
 KratosALEApllication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Pooyan Dadvand, Riccardo Rossi, Andreas Winterstein
                     pooyan@cimne.upc.edu
                     rrossi@cimne.upc.edu
                     a.winterstein@tum.de
- CIMNE (International Center for Numerical Methods in Engineering),
  Gran Capita' s/n, 08034 Barcelona, Spain
- Chair of Structural Analysis, Technical University of Munich
  Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================


/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: A.Winterstein@tum.de $
 *  Date:                $Date: June 2016 $
 *  Revision:            $Revision: 1.5 $
 * ***************************************************************************/

#if !defined(KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY )
#define  KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY

/* System includes */
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp"

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "ale_application.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
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
template<class TSparseSpace, 
         class TDenseSpace, 
         class TLinearSolver
         >
class LaplacianMeshMovingStrategy 
    : public SolvingStrategy<TSparseSpace,TDenseSpace, TLinearSolver> 
{

 public:
  /**@name Type Definitions */
  /*@{ */

  /** Counted pointer of ClassName */
  KRATOS_CLASS_POINTER_DEFINITION( LaplacianMeshMovingStrategy );

  typedef SolvingStrategy< TSparseSpace,TDenseSpace,TLinearSolver > BaseType;
  typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
  typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;

  /*@} */
  /**@name Life Cycle
   */
  /*@{ */

  /************************************* CONSTRUCTOR *********************************/
  /***********************************************************************************/
  LaplacianMeshMovingStrategy(
      ModelPart& model_part,
      typename TLinearSolver::Pointer pNewLinearSolver,
      int TimeOrder = 1,
      bool ReformDofSetAtEachStep = false,
      int EchoLevel = 0
    )
    :SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part) {

    KRATOS_TRY;

    mReformDofSetAtEachStep = ReformDofSetAtEachStep;
    mEchoLevel = EchoLevel;

    mTimeOrder = TimeOrder;
    bool compute_reactions = false;
    bool calculate_norm_dx_flag = false;
  
    typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
    ( new ResidualBasedIncrementalUpdateStaticScheme
        <TSparseSpace, TDenseSpace>());

    GenerateMeshPart();

    
    mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
    (new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(pNewLinearSolver) );

    mStrategy = typename BaseType::Pointer
    ( new ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver >
    (*mpMeshModelPart,pscheme,pNewLinearSolver,mpBuilderAndSolver,compute_reactions,mReformDofSetAtEachStep,calculate_norm_dx_flag) );
    
    mStrategy->SetEchoLevel(mEchoLevel);
    /*
    
    typedef typename BuilderAndSolver
    <TSparseSpace,TDenseSpace,TLinearSolver>
    ::Pointer BuilderSolverTypePointer;

    typedef typename Kratos::VariableComponent
    <Kratos::VectorComponentAdaptor<Kratos
    ::array_1d<double, 3> > > VarComponent;

    mpBuilderAndSolverX = typename TBuilderAndSolverType::Pointer
    (new ResidualBasedEliminationBuilderAndSolverComponentwise
        <TSparseSpace,TDenseSpace,TLinearSolver, VarComponent>
        ( pNewLinearSolver, MESH_DISPLACEMENT_X));

    mpBuilderAndSolverY = typename TBuilderAndSolverType::Pointer
    (new ResidualBasedEliminationBuilderAndSolverComponentwise
        <TSparseSpace,TDenseSpace,TLinearSolver, VarComponent>
        ( pNewLinearSolver, MESH_DISPLACEMENT_Y));

    mpBuilderAndSolverZ = typename TBuilderAndSolverType::Pointer
    (new ResidualBasedEliminationBuilderAndSolverComponentwise
        <TSparseSpace,TDenseSpace,TLinearSolver, VarComponent>
        ( pNewLinearSolver, MESH_DISPLACEMENT_Z));


    mStrategy_X = typename BaseType::Pointer
    (new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
        (*mpMeshModelPart,pscheme,pNewLinearSolver,mpBuilderAndSolverX,
            compute_reactions,mReformDofSetAtEachStep,calculate_norm_dx_flag));

    mStrategy_Y = typename BaseType::Pointer
    (new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
        (*mpMeshModelPart,pscheme,pNewLinearSolver,mpBuilderAndSolverY,
            compute_reactions,mReformDofSetAtEachStep,calculate_norm_dx_flag));

    mStrategy_Z = typename BaseType::Pointer
    (new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
        (*mpMeshModelPart,pscheme,pNewLinearSolver,mpBuilderAndSolverZ,
            compute_reactions,mReformDofSetAtEachStep,calculate_norm_dx_flag));
    */
    KRATOS_CATCH("")
  }

  virtual ~LaplacianMeshMovingStrategy() {

  };
   

  double Solve() override
  {
    KRATOS_TRY;

    ProcessInfo& rCurrentProcessInfo = (mpMeshModelPart)->GetProcessInfo();

    rCurrentProcessInfo[TIME] = BaseType::GetModelPart().GetProcessInfo()[TIME];
    rCurrentProcessInfo[DELTA_TIME] = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];

    // Setting mesh to initial configuration
    SetMeshToInitialConfiguration();

    //X DIRECTION
    //rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
    //mStrategy_X->Solve();

    //Y DIRECTION
    //rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
    //mStrategy_Y->Solve();

    //Z DIRECTION
    //rCurrentProcessInfo[FRACTIONAL_STEP] = 3;
    //mStrategy_Z->Solve();

    rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
    // Solve for the mesh movement
    mStrategy->Solve();

    rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
    // Solve for the mesh movement
    mStrategy->Solve();

    rCurrentProcessInfo[FRACTIONAL_STEP] = 3;
    // Solve for the mesh movement
    mStrategy->Solve();

    // Update FEM-base
    CalculateMeshVelocities();
    MoveMesh();

    return 0.0;

    KRATOS_CATCH("");
  }

  void CalculateMeshVelocities()
  {
    KRATOS_TRY;

    double DeltaTime = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];

    if (DeltaTime <= 0.0)
        KRATOS_ERROR<<"Invalid DELTA_TIME."<<std::endl;

    double coeff = 1/DeltaTime;

    if( mTimeOrder == 1)
    {
      for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin();
          i != (*mpMeshModelPart).NodesEnd(); ++i) {

        array_1d<double,3>& mesh_v = (i)
        ->FastGetSolutionStepValue(MESH_VELOCITY);
        array_1d<double,3>& disp = (i)
        ->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        array_1d<double,3>& dispold = (i)
        ->FastGetSolutionStepValue(MESH_DISPLACEMENT,1);
        noalias(mesh_v) = disp - dispold;
        mesh_v *= coeff;
      }
    }
    else if (mTimeOrder == 2)
    {
      double c1 = 1.50*coeff;
      double c2 = -2.0*coeff;
      double c3 = 0.50*coeff;

      for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin();
          i != (*mpMeshModelPart).NodesEnd(); ++i) {

        array_1d<double,3>& mesh_v = (i)
        ->FastGetSolutionStepValue(MESH_VELOCITY);
        noalias(mesh_v) = c1 * (i)
        ->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        noalias(mesh_v) += c2 * (i)
        ->FastGetSolutionStepValue(MESH_DISPLACEMENT,1);
        noalias(mesh_v) += c3 * (i)
        ->FastGetSolutionStepValue(MESH_DISPLACEMENT,2);
      }
    }
    else 
    {
      KRATOS_ERROR<<"Wrong TimeOrder: Acceptable values are: 1 and 2"<<std::endl;
    }
    
    KRATOS_CATCH("");
  }

  void MoveMesh() override
  {
    for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
            i != BaseType::GetModelPart().NodesEnd(); ++i)
    {
        (i)->X() = (i)->X0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_X);
        (i)->Y() = (i)->Y0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_Y);
        (i)->Z() = (i)->Z0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_Z);
    }
  }

  void SetMeshToInitialConfiguration()
  {
  for(ModelPart::NodeIterator i = (*mpMeshModelPart ).NodesBegin();
      i != (*mpMeshModelPart ).NodesEnd(); ++i) {

      (i)->X() = (i)->X0();
      (i)->Y() = (i)->Y0();
      (i)->Z() = (i)->Z0();
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
  ModelPart::Pointer mpMeshModelPart;

  typename BaseType::Pointer mStrategy_X;
  typename BaseType::Pointer mStrategy_Y;
  typename BaseType::Pointer mStrategy_Z;
  typename BaseType::Pointer mStrategy;
  typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;
  typename TBuilderAndSolverType::Pointer mpBuilderAndSolverX;
  typename TBuilderAndSolverType::Pointer mpBuilderAndSolverY;
  typename TBuilderAndSolverType::Pointer mpBuilderAndSolverZ;
  


  bool mReformDofSetAtEachStep;
  int mTimeOrder;
  int mEchoLevel;
  /*@} */
  /**@name Private Operators*/
  /*@{ */

  /*@} */
  /**@name Private Operations*/
  /*@{ */

  void GenerateMeshPart()
  {
    mpMeshModelPart = ModelPart::Pointer( new ModelPart("MeshPart",1) );

    //initializing mesh nodes
    mpMeshModelPart->Nodes() = BaseType::GetModelPart().Nodes();

    //creating mesh elements
    ModelPart::ElementsContainerType& MeshElems = mpMeshModelPart->Elements();
    Element::Pointer pElem;

    for(ModelPart::ElementsContainerType::iterator it
        = BaseType::GetModelPart().ElementsBegin();

        it != BaseType::GetModelPart().ElementsEnd(); ++it) {

      pElem = Element::Pointer(new LaplacianMeshMovingElement(
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
