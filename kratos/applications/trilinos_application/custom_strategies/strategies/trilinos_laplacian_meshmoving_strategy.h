/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: A.Winterstein@tum.de $
 *  Date:                $Date: June 2016 $
 *  Revision:            $Revision: 1.5 $
 * ***************************************************************************/


#if !defined(KRATOS_NEW_TRILINOS_LAPLACIAN_MESHMOVING_STRATEGY )
#define  KRATOS_NEW_TRILINOS_LAPLACIAN_MESHMOVING_STRATEGY


/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/communicator.h"
#include "includes/mpi_communicator.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "../ALEapplication/ale_application.h"
#include "../ALEapplication/custom_elements/laplacian_meshmoving_element.h"


/* Trilinos includes */
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
#include "custom_utilities/parallel_fill_communicator.h"

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

/// Short class definition.
/**   Detail class definition.



*/
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class TrilinosLaplacianMeshMovingStrategy
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosLaplacianMeshMovingStrategy );

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
    TrilinosLaplacianMeshMovingStrategy(
        Epetra_MpiComm& Comm,
        ModelPart& model_part,
        typename TLinearSolver::Pointer pNewLinearSolver,
        int dimension = 3,
        int velocity_order = 1,
        bool reform_dof_at_every_step = true
    )
        : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part)
    {
        KRATOS_TRY

        // Passed variables
        mdimension = dimension;
        mvel_order = velocity_order;
        mreform_dof_at_every_step = reform_dof_at_every_step;

        // Definitions for trilinos
        int guess_row_size;
        guess_row_size = 15;

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        // Generating Mesh Part
        GenerateMeshPart();

        typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
        typename SchemeType::Pointer pscheme = typename SchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

        BuilderSolverTypePointer aux_var_build = BuilderSolverTypePointer(new TrilinosBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(Comm, guess_row_size, pNewLinearSolver) );

        mstrategy = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >(*mpMeshModelPart,pscheme,pNewLinearSolver,aux_var_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        mstrategy->SetEchoLevel(2);

        for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ;
                i != (*mpMeshModelPart).NodesEnd() ; ++i)
        {
            if(!i->IsFixed(DISPLACEMENT_X))
                (i)->GetSolutionStepValue(DISPLACEMENT_X) = 0.00;
            if(!i->IsFixed(DISPLACEMENT_Y))
                (i)->GetSolutionStepValue(DISPLACEMENT_Y) = 0.00;
            if(!i->IsFixed(DISPLACEMENT_Z))
                (i)->GetSolutionStepValue(DISPLACEMENT_Z) = 0.00;
        }

        KRATOS_CATCH("");
    }



    /** Destructor.
    */
    virtual ~TrilinosLaplacianMeshMovingStrategy() {}

    /** Destructor.
    */


    //*********************************************************************************
    //*********************************************************************************
    double Solve()
    {
        KRATOS_TRY;

        // Mesh has to be regenerated in each solving step


        ProcessInfo& rCurrentProcessInfo = (mpMeshModelPart)->GetProcessInfo();

        // Updating the time
        rCurrentProcessInfo[TIME] = BaseType::GetModelPart().GetProcessInfo()[TIME];
        rCurrentProcessInfo[DELTA_TIME] = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];

        // Fixing Dofs As Needed
        for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ;
                i != (*mpMeshModelPart).NodesEnd() ; ++i)
        {
            if(i->IsFixed(DISPLACEMENT_X))
                (i)->Fix(AUX_MESH_VAR);
            if(i->IsFixed(DISPLACEMENT_Y))
                (i)->Fix(AUX_MESH_VAR);
            if(i->IsFixed(DISPLACEMENT_Z))
                (i)->Fix(AUX_MESH_VAR);
        }
        // X DIRECTION
        rCurrentProcessInfo[FRACTIONAL_STEP] = 1; //laplacian mesh moving type corresponds to -1
        for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ;
                i != (*mpMeshModelPart).NodesEnd() ; ++i)
        {
            (i)->FastGetSolutionStepValue(AUX_MESH_VAR) = (i)->GetSolutionStepValue(DISPLACEMENT_X);
        }
        mstrategy->Solve();
        for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ;
                i != (*mpMeshModelPart).NodesEnd() ; ++i)
        {
            (i)->GetSolutionStepValue(DISPLACEMENT_X) = (i)->FastGetSolutionStepValue(AUX_MESH_VAR);
        }

        //Y DIRECTION
        if(mdimension > 1)
        {
            rCurrentProcessInfo[FRACTIONAL_STEP] = 2; //laplacian mesh moving type corresponds to -1
            for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ;
                    i != (*mpMeshModelPart).NodesEnd() ; ++i)
            {
                (i)->FastGetSolutionStepValue(AUX_MESH_VAR) = (i)->GetSolutionStepValue(DISPLACEMENT_Y);
            }
            mstrategy->Solve();
            for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ;
                    i != (*mpMeshModelPart).NodesEnd() ; ++i)
            {
                (i)->GetSolutionStepValue(DISPLACEMENT_Y) = (i)->FastGetSolutionStepValue(AUX_MESH_VAR);
            }
        }

        //Z DIRECTION
        if(mdimension > 2)
        {
            rCurrentProcessInfo[FRACTIONAL_STEP] = 3; //laplacian mesh moving type corresponds to -1
            for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ;
                    i != (*mpMeshModelPart).NodesEnd() ; ++i)
            {
                (i)->FastGetSolutionStepValue(AUX_MESH_VAR) = (i)->GetSolutionStepValue(DISPLACEMENT_Z);
            }
            mstrategy->Solve();
            for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ;
                    i != (*mpMeshModelPart).NodesEnd() ; ++i)
            {
                (i)->GetSolutionStepValue(DISPLACEMENT_Z) = (i)->FastGetSolutionStepValue(AUX_MESH_VAR);
            }
        }

		// Update FEM database
        CalculateMeshVelocities();
        BaseType::MoveMesh();

        //clearing the system if needed
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
        mpMeshModelPart->SetProcessInfo(BaseType::GetModelPart().pGetProcessInfo());
        mpMeshModelPart->SetProperties(BaseType::GetModelPart().pProperties());

        // Create a communicator for the new model part and copy the partition information about nodes.
        Communicator& rReferenceComm = BaseType::GetModelPart().GetCommunicator();
        typename Communicator::Pointer pMeshMPIComm = typename Communicator::Pointer( new MPICommunicator( &(BaseType::GetModelPart().GetNodalSolutionStepVariablesList()) ) );
        pMeshMPIComm->SetNumberOfColors( rReferenceComm.GetNumberOfColors() ) ;
        pMeshMPIComm->NeighbourIndices() = rReferenceComm.NeighbourIndices();
        pMeshMPIComm->LocalMesh().SetNodes( rReferenceComm.LocalMesh().pNodes() );
        pMeshMPIComm->InterfaceMesh().SetNodes( rReferenceComm.InterfaceMesh().pNodes() );
        pMeshMPIComm->GhostMesh().SetNodes( rReferenceComm.GhostMesh().pNodes() );

        // Setup communication plan
        for (unsigned int i = 0; i < rReferenceComm.GetNumberOfColors(); i++)
        {
            pMeshMPIComm->pInterfaceMesh(i)->SetNodes( rReferenceComm.pInterfaceMesh(i)->pNodes() );
            pMeshMPIComm->pLocalMesh(i)->SetNodes( rReferenceComm.pLocalMesh(i)->pNodes() );
            pMeshMPIComm->pGhostMesh(i)->SetNodes( rReferenceComm.pGhostMesh(i)->pNodes() );
        }
        mpMeshModelPart->SetCommunicator( pMeshMPIComm );

        // Creating mesh elements
        ModelPart::ElementsContainerType& MeshElems = mpMeshModelPart->Elements();
        Element::Pointer pElem;

            for(ModelPart::ElementsContainerType::iterator it =  BaseType::GetModelPart().ElementsBegin();
                    it != BaseType::GetModelPart().ElementsEnd(); it++)
            {
                pElem = Element::Pointer(new LaplacianMeshMovingElement(
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
    TrilinosLaplacianMeshMovingStrategy(const TrilinosLaplacianMeshMovingStrategy& Other);


    /*@} */

}; /* Class TrilinosLaplacianMeshMovingStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_NEW_TRILINOS_LAPLACIAN_MESHMOVING_STRATEGY  defined */

