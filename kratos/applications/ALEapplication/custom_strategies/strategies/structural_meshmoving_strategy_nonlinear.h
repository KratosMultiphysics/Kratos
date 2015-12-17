/* *********************************************************
*
*   Last Modified by:    $Author: AMini $
*   Date:                $Date: Mai 2015 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/


#if !defined(KRATOS_NEW_STRUCTURAL_MESHMOVING_STRATEGY_NONLINEAR )
#define  KRATOS_NEW_STRUCTURAL_MESHMOVING_STRATEGY_NONLINEAR


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_elements/structural_meshmoving_element_nonlinear.h"
#include "ale_application.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"


//includes for non-linear solver
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/schemes/residualbased_incremental_aitken_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "spaces/ublas_space.h"


namespace Kratos
{

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

/** Short class definition.

Detail class definition.

Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

The non-linear system of equations is solved by a Newton-Raphson iteration.

 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class StructuralMeshMovingStrategyNonlinear
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( StructuralMeshMovingStrategyNonlinear );

    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;



    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    StructuralMeshMovingStrategyNonlinear(
        ModelPart& model_part,
        typename TLinearSolver::Pointer pNewLinearSolver,
        int dimension = 3,
        int velocity_order = 1,
        bool reform_dof_at_every_step = true,
        double NonLinearTol = 10e-5,    //default value
        int MaxIter = 50                //default value

    )
        : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part)
    {
        KRATOS_TRY

        //Generating Mesh Part
        GenerateMeshPart(dimension);

        mdimension = dimension;
        mvel_order = velocity_order;
        mreform_dof_at_every_step = reform_dof_at_every_step;
        mtol = NonLinearTol;
        mmax_it = MaxIter;


        typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
        typename SchemeType::Pointer pscheme = typename SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace >() );
        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;


        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        //bool CalculateNormDxFlag = false;


        BuilderSolverTypePointer pBuilderSolver = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(pNewLinearSolver) );

        // Convergence criteria
        typedef typename ConvergenceCriteria< TSparseSpace, TDenseSpace >::Pointer ConvergenceCriteriaPointerType;
        typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;


        double NearlyZero = 1.0e-20;
        ConvergenceCriteriaPointerType pConvCriteria = ConvergenceCriteriaPointerType( new ResidualCriteria<TSparseSpace,TDenseSpace>(NonLinearTol,NearlyZero) );

        bool MoveMesh = false;

        mstrategy = StrategyPointerType( new ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(*mpMeshModelPart,pscheme,pNewLinearSolver,pConvCriteria,pBuilderSolver,MaxIter,CalculateReactions,ReformDofAtEachIteration,MoveMesh));
        mstrategy->SetEchoLevel(2);

        //====================================================

        KRATOS_CATCH("")
    }



    /** Destructor.
    */
    virtual ~StructuralMeshMovingStrategyNonlinear() {}

    /** Destructor.
    */


    //*********************************************************************************
    //*********************************************************************************
    double Solve()
    {
        KRATOS_TRY;

        // Setting mesh to initial configuration
        for(ModelPart::NodeIterator i = (*mpMeshModelPart ).NodesBegin();
            i != (*mpMeshModelPart ).NodesEnd(); ++i) {

          (i)->X() = (i)->X0();
          (i)->Y() = (i)->Y0();
          (i)->Z() = (i)->Z0();
        }


        // Solve for the mesh movement
        mstrategy->Solve();

        //copy back

        // Update FEM-base
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
        if( mvel_order == 1) //mesh velocity calculated as (x(n+1)-x(n))/Dt
        {
            for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ;
                    i != (*mpMeshModelPart).NodesEnd() ; ++i)
            {
                array_1d<double,3>& mesh_v = (i)->FastGetSolutionStepValue(MESH_VELOCITY);
                array_1d<double,3>& disp = (i)->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double,3>& dispold = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
                noalias(mesh_v) =  disp - dispold;
                mesh_v *= coeff;
            }
        }
        else //mesh velocity calculated as (3*x(n+1)-4*x(n)+x(n-1))/(2*Dt)
        {
            double c1 = 1.50*coeff;
            double c2 = -2.0*coeff;
            double c3 = 0.50*coeff;

            for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ;
                    i != (*mpMeshModelPart).NodesEnd() ; ++i)
            {
                array_1d<double,3>& mesh_v = (i)->FastGetSolutionStepValue(MESH_VELOCITY);
                noalias(mesh_v) =  c1 * (i)->FastGetSolutionStepValue(DISPLACEMENT);
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

    int mdimension;
    int mvel_order;
    bool mreform_dof_at_every_step;
    double mtol;
    int mmax_it;

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */

    void GenerateMeshPart(int dimension)
    {
        mpMeshModelPart = ModelPart::Pointer( new ModelPart("MeshPart",1) );

        mpMeshModelPart->SetProcessInfo( BaseType::GetModelPart().pGetProcessInfo() );
        mpMeshModelPart->SetBufferSize( BaseType::GetModelPart().GetBufferSize() );
        mpMeshModelPart->SetProperties( BaseType::GetModelPart().pProperties() );

        // Initializing mesh nodes
        mpMeshModelPart->Nodes() = BaseType::GetModelPart().Nodes();

        // Removing existing mesh conditions
        mpMeshModelPart->Conditions().clear();


        // Creating mesh elements
        ModelPart::ElementsContainerType& MeshElems = mpMeshModelPart->Elements();
        Element::Pointer pElem;

        if(dimension == 2)
            for(ModelPart::ElementsContainerType::iterator it =  BaseType::GetModelPart().ElementsBegin();
                    it != BaseType::GetModelPart().ElementsEnd(); it++)
            {
                pElem = Element::Pointer(new StructuralMeshMovingElementNonlinear<2>(
                                             (*it).Id(),
                                             (*it).pGetGeometry(),
                                             (*it).pGetProperties() ) );
                MeshElems.push_back(pElem);
            }
        else if(dimension == 3)
            for(ModelPart::ElementsContainerType::iterator it =  BaseType::GetModelPart().ElementsBegin();
                    it != BaseType::GetModelPart().ElementsEnd(); it++)
            {
                pElem = Element::Pointer(new StructuralMeshMovingElementNonlinear<3>(
                                             (*it).Id(),
                                             (*it).pGetGeometry(),
                                             (*it).pGetProperties() ) );
                MeshElems.push_back(pElem);
            }
    }

//    void ReGenerateMeshPart()
//    {
//        std::cout << "regenerating elements for the mesh motion scheme" << std::endl;

//        mpMeshModelPart->SetProcessInfo( BaseType::GetModelPart().pGetProcessInfo() );
//        mpMeshModelPart->SetBufferSize( BaseType::GetModelPart().GetBufferSize() );
//        mpMeshModelPart->SetProperties( BaseType::GetModelPart().pProperties() );

//        // Initializing mesh nodes
//        mpMeshModelPart->Nodes().clear();
//        mpMeshModelPart->Nodes() = BaseType::GetModelPart().Nodes();

//        //creating mesh elements
//        ModelPart::ElementsContainerType& MeshElems = mpMeshModelPart->Elements();
//        Element::Pointer pElem;

//        MeshElems.clear();
//        MeshElems.reserve( MeshElems.size() );

//        // Removing existing mesh conditions
//        mpMeshModelPart->Conditions().clear();

//        if(mdimension == 2)
//            for(ModelPart::ElementsContainerType::iterator it =  BaseType::GetModelPart().ElementsBegin();
//                    it != BaseType::GetModelPart().ElementsEnd(); it++)
//            {
//                pElem = Element::Pointer(new StructuralMeshMovingElementNonlinear<2>(
//                                             (*it).Id(),
//                                             (*it).pGetGeometry(),
//                                             (*it).pGetProperties() ) );
//                MeshElems.push_back(pElem);
//            }
//        else if(mdimension == 3)
//            for(ModelPart::ElementsContainerType::iterator it =  BaseType::GetModelPart().ElementsBegin();
//                    it != BaseType::GetModelPart().ElementsEnd(); it++)
//            {
//                pElem = Element::Pointer(new StructuralMeshMovingElementNonlinear<3>(
//                                             (*it).Id(),
//                                             (*it).pGetGeometry(),
//                                             (*it).pGetProperties() ) );
//                MeshElems.push_back(pElem);
//            }
//    }

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
    StructuralMeshMovingStrategyNonlinear(const StructuralMeshMovingStrategyNonlinear& Other);


    /*@} */

}; /* Class StructuralMeshMovingStrategyNonlin */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_NEW_STRUCTURAL_MESHMOVING_STRATEGY  defined */

