#ifndef KRATOS_TRILINOS_BUILDER_AND_SOLVER_ML_PERIODIC_H
#define KRATOS_TRILINOS_BUILDER_AND_SOLVER_ML_PERIODIC_H

/* System includes */
#include <set>
/* #include <omp.h> */

/* External includes */
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp"


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_mixed.h"

//trilinos includes
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

//aztec solver includes
#include "AztecOO.h"

#include "Amesos.h"
// #include "AmesosClassType.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelPreconditioner.h"

namespace Kratos
{

///@addtogroup TrilinosApplication
///@{

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

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class TrilinosBuilderAndSolverMLPeriodic
        : public TrilinosBuilderAndSolverMLmixed< TSparseSpace,TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosBuilderAndSolverMLPeriodic );


    typedef TrilinosBuilderAndSolverMLmixed<TSparseSpace,TDenseSpace, TLinearSolver > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    TrilinosBuilderAndSolverMLPeriodic(Epetra_MpiComm& Comm,
                                       int guess_row_size,
                                       int dim,
                                       int ndofs,
                                       typename TLinearSolver::Pointer pNewLinearSystemSolver,
                                       const Kratos::Variable<int>& PeriodicIdVar):
        TrilinosBuilderAndSolverMLmixed< TSparseSpace,TDenseSpace,TLinearSolver >(Comm,guess_row_size,dim,pNewLinearSystemSolver),
        mNumDof(ndofs),
        mPeriodicIdVar(PeriodicIdVar)
    {}



    /** Destructor.
    */
    virtual ~TrilinosBuilderAndSolverMLPeriodic() {}

    /*@} */
    /**@name Operations */
    /*@{ */


    /// Assign an Equation Id to all degrees of freedom in the system.
    /**
     * To properly set up the periodic conditions, it is assumed that all nodes on
     * a periodic boundary have the Id of their periodic "image" stored as Node->GetSolutionStepValue(mPeriodicIdVar,0).
     * @note It is assumed that the problem's partition takes the presence of periodic conditions into account, that is,
     * all processes that own periodic conditions know the nodes on both ends (either as local or ghost nodes).
     * @param rModelPart The problem's ModelPart
     */
    virtual void SetUpSystem(ModelPart &rModelPart)
    {
        KRATOS_TRY;

        // This sort helps us preserve consistency of dof ordering across processors (as searching for a dof in
        // the nodal list of dofs will sort the list first).
        // If/when the DofList becomes a static array, this step can be skipped.
        for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)
            itNode->GetDofs().Sort();

        unsigned int Rank = this->mrComm.MyPID();

        // Count the Dofs on this partition (on periodic node pairs, only the dofs on the node with higher Id are counted)
        int DofCount = 0;

        for (typename DofsArrayType::iterator itDof = this->mDofSet.begin(); itDof != this->mDofSet.end(); ++itDof)
            if ( (itDof->GetSolutionStepValue(PARTITION_INDEX) == Rank) && ( itDof->GetSolutionStepValue(mPeriodicIdVar) < static_cast<int>(itDof->Id())) )
                DofCount++;

        // Comunicate the Dof counts to set a global numbering
        int DofOffset;

        int ierr = this->mrComm.ScanSum(&DofCount,&DofOffset,1);
        if (ierr != 0)
            KRATOS_THROW_ERROR(std::runtime_error,"In TrilinosBuilderAndSolverMLPeriodic::SetUpSystem: Found Epetra_MpiComm::ScanSum failure with error code ",ierr);

        DofOffset -= DofCount;
        this->mFirstMyId = DofOffset;

        // Assing Id to all local nodes except for the second nodes on a periodic pair
        for (typename DofsArrayType::iterator itDof = this->mDofSet.begin(); itDof != this->mDofSet.end(); ++itDof)
            if ( (itDof->GetSolutionStepValue(PARTITION_INDEX) == Rank) && ( itDof->GetSolutionStepValue(mPeriodicIdVar) < static_cast<int>(itDof->Id()) ) )
                itDof->SetEquationId(DofOffset++);

        // Synchronize Dof Ids across processes
        rModelPart.GetCommunicator().SynchronizeDofs();

        // Once the non-repeated EquationId are assigned and synchronized, copy the EquationId to the second node on each periodic pair
        for (ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond)
        {
            // PeriodicCondition always have exactly 2 nodes
            ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
            if (rGeom.PointsNumber() == 2)
            {
                int Node0 = rGeom[0].Id();
                int Node0Pair = rGeom[0].FastGetSolutionStepValue(mPeriodicIdVar);

                int Node1 = rGeom[1].Id();
                int Node1Pair = rGeom[1].FastGetSolutionStepValue(mPeriodicIdVar);

                // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                if ( ( Node0 == Node1Pair ) && ( Node1 == Node0Pair ) )
                {
                    if ( Node0 < Node0Pair ) // If Node0 is the one with lower Id (the one that does not have an EquationId yet)
                        CopyEquationId(rGeom[1],rGeom[0]);
                    else
                        CopyEquationId(rGeom[0],rGeom[1]);
                }
            }
        }

        // Synchronize Dof Ids across processes a second time to transfer the Ids of periodic nodes.
        // Note that a custom communication routine is used, as the one provided by the communicator synchronizes all ghost values to the
        // value on the owner process, while here the EquationId value is given by the owner of the PeriodicCondition,
        // which is not necessarily the same as the owner of the node.
        CommunicateEquationId(rModelPart.GetCommunicator());

        // Store local and gloabal system size
        int TotalDofNum;
        ierr = this->mrComm.SumAll(&DofCount,&TotalDofNum,1);
        if (ierr != 0)
            KRATOS_THROW_ERROR(std::runtime_error,"In TrilinosBuilderAndSolverMLPeriodic::SetUpSystem: Found Epetra_MpiComm::SumAll failure with error code ",ierr);

        this->mLocalSystemSize = DofCount;
        this->mEquationSystemSize = TotalDofNum;
        this->mLastMyId = DofOffset;

        /* the following prints can be helpful for debug
        for (typename DofsArrayType::iterator itDof = this->mDofSet.begin(); itDof != this->mDofSet.end(); ++itDof)
            if (itDof->GetSolutionStepValue(PARTITION_INDEX) == Rank)
            {
                if (itDof->EquationId() < this->mFirstMyId || itDof->EquationId() >= this->mLastMyId)
                {
                    KRATOS_WATCH(Rank);
                    KRATOS_WATCH(itDof->EquationId());
                    KRATOS_WATCH(this->mFirstMyId);
                    KRATOS_WATCH(this->mLastMyId);
                }
            }
        */

        KRATOS_CATCH("");
    }

    virtual void SystemSolveML(
            TSystemMatrixType& A,
            TSystemVectorType& Dx,
            TSystemVectorType& b,
            ModelPart& r_model_part
            )
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            //***********************************************
            //new attempt
            Epetra_LinearProblem AztecProblem(&A, &Dx, &b);

            Epetra_Vector scaling_vect(A.RowMap());
            A.InvColSums(scaling_vect);
            AztecProblem.LeftScale(scaling_vect);

            AztecOO solver(AztecProblem);

            // create an empty parameter list for ML options
            Teuchos::ParameterList MLList;

            // set defaults for classic smoothed aggregation with heavy smoothers
            // (of domain decomposition type, i.e. one-level Schwarz with incomplete
            // factorizations on each subdomain/process)
            // We need to define the solvers on each subdomain (== processor).
            // Here we use an incomplete LU factorization, with no fill-in
            // and no overlap. To that aim, we use Aztec's preconditioning function.
            // Aztec requires two more vectors. Note: the following options and params
            // will be used ONLY for the smoother, and will NOT affect the Aztec solver
            // NOTE: to use exact solvers change to AZ_lu (requires AztecOO configured
            // with option--enable-aztecoo-azlu), of use IFPACK smoothers
            // (requires Trilinos to be built with options--enable-ifpack --enable-amesos)

            int options[AZ_OPTIONS_SIZE];
            double params[AZ_PARAMS_SIZE];
            AZ_defaults(options, params);
            options[AZ_precond] = AZ_dom_decomp;
            options[AZ_subdomain_solve] = AZ_ilu;
            options[AZ_graph_fill] = 0;
            options[AZ_overlap] = 0;

            // SetDefaults() will call AZ_defaults(options,params), and will also set the
            // preconditioner as `AZ_dom_decomp'.
            // NOTE THAT THE `options' AND `params' VECTORS ARE NOT COPIED into
            // the list, only the pointers is stored, so do not delete options
            // and params before the end of the linear system solution!
            // Alternatively, you can also call SetDefaults() without passing
            // `options' and `params.' This way, the code will allocate a int
            // and a double vector, that must be freed by the user.
            // `DD' means to set default values for domain decomposition
            // preconditioners

            ML_Epetra::SetDefaults("NSSA", MLList, options, params);
            //                ML_Epetra::SetDefaults("DD", MLList, options, params);

            // Overwrite some parameters. Please refer to the user's guide
            // for more information
            // Some parameters are reported here to better explain the process
            // even if they are as defaults.
            // NOTE: To use `METIS' as aggregation scheme, you need to configure
            // ML with the option --with-ml_metis. Otherwise, the code will
            // creates aggregates containing all the local nodes (that is,
            // the dimension of the coarse problem will be equal to the
            // number of processors)

            //                MLList.set("aggregation: type", "METIS");
            //                MLList.set("smoother: type", "Aztec");

            //                  MLList.set("aggregation: nodes per aggregate", 128);
            //                  MLList.set("smoother: pre or post", "pre");
            //                  MLList.set("coarse: type","Amesos-KLU");



            MLList.set("PDE equations", mNumDof);
            MLList.set("null space: add default vectors", true);
            MLList.set("aggregation: type","Uncoupled");

            // Create the preconditioning object. We suggest to use `new' and
            // `delete' because the destructor contains some calls to MPI (as
            // required by ML and possibly Amesos). This is an issue only if the
            // destructor is called **after** MPI_Finalize().

            ML_Epetra::MultiLevelPreconditioner* MLPrec =
                new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);

            // tell AztecOO to use this preconditioner, then solve
            solver.SetPrecOperator(MLPrec);

            // =========================== end of ML part =============================

            // Instruct AztecOO to use GMRES with estimation of the condition
            // number. Also, requires output every 32 iterations
            // Then, solve with 500 iterations and 1e-8 as tolerance on the
            // relative residual

            solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
            solver.SetAztecOption(AZ_output, AZ_none);
            solver.SetAztecOption(AZ_kspace, 200);
            solver.Iterate(500, 1e-8);

            // delete the preconditioner. Do it BEFORE MPI_Finalize
            delete MLPrec;


            //***********************************************


            //
            //                Teuchos::ParameterList MLList;
            //                MLList.set("energy minimization: enable", true);
            //                MLList.set("energy minimization: type", 3); // 1,2,3 cheap -> expensive
            //                MLList.set("aggregation: block scaling", false);
            //                MLList.set("aggregation: type", "Uncoupled");
            //                //				MLList.set("smoother: type (level 0)","symmetric Gauss-Seidel");
            //                MLList.set("smoother: type (level 0)", "Jacobi");
            //                MLList.set("smoother: sweeps (level 0)", 2);
            //                MLList.set("smoother: damping factor (level 0)", 0.89);
            //
            //                MLList.set("PDE equations", numdf);
            //                MLList.set("null space: dimension", dimns);
            //                MLList.set("null space: type", "pre-computed");
            //                MLList.set("null space: add default vectors", false);
            //                MLList.set("null space: vectors", nullsp);
            //
            //                // create the preconditioner
            //                ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);
            //
            //                // create an AztecOO solver
            //                AztecOO Solver(AztecProblem);
            //
            //                // set preconditioner and solve
            //                Solver.SetPrecOperator(MLPrec);
            //                Solver.SetAztecOption(AZ_solver, AZ_gmres);
            //                Solver.SetAztecOption(AZ_kspace, 200);
            //                Solver.SetAztecOption(AZ_output, 15); //SetAztecOption(AZ_output, AZ_none);
            //
            //                int mmax_iter = 300;
            //                Solver.Iterate(mmax_iter, 1e-9);
            //                delete MLPrec;

            //		BaseType::mpLinearSystemSolver->Solve(A,Dx,b);

        }
        else
        {
            TSparseSpace::SetToZero(Dx);
        }

        //prints informations about the current time
        if (this->GetEchoLevel() > 1)
        {
            if (r_model_part.GetCommunicator().MyPID() == 0)
                std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
        }

        KRATOS_CATCH("")

    }
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


    /*@} */
    /**@name Operators
    */
    /*@{ */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    const int mNumDof;

    const Kratos::Variable<int>& mPeriodicIdVar;

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */

    /// Duplicate EquationIds to the second node on each periodic pair
    void CopyEquationId(ModelPart::NodeType& rOrigin,
                        ModelPart::NodeType& rDest)
    {
        for (typename DofsArrayType::iterator itDof = rOrigin.GetDofs().begin(); itDof != rOrigin.GetDofs().end(); itDof++)
            rDest.pGetDof( itDof->GetVariable() )->SetEquationId( itDof->EquationId() );
    }

    /// Send the Equation Id of periodic nodes to the owner of the node, so it can be sync'ed across processes.
    void CommunicateEquationId(Communicator& rComm)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int destination = 0;

        Communicator::NeighbourIndicesContainerType& neighbours_indices = rComm.NeighbourIndices();

        for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
            if ((destination = neighbours_indices[i_color]) >= 0)
            {
                Communicator::NodesContainerType& r_origin_nodes = rComm.GhostMesh(i_color).Nodes();
                Communicator::NodesContainerType& r_dest_nodes = rComm.LocalMesh(i_color).Nodes();

                // Calculating send and received buffer size
                unsigned int send_buffer_size = 0;
                unsigned int receive_buffer_size = 0;

                for (Communicator::NodesContainerType::iterator i_node = r_origin_nodes.begin(); i_node != r_origin_nodes.end(); ++i_node)
                    send_buffer_size += i_node->GetDofs().size();

                for (Communicator::NodesContainerType::iterator i_node = r_dest_nodes.begin(); i_node != r_dest_nodes.end(); ++i_node)
                    receive_buffer_size += i_node->GetDofs().size();

                unsigned int position = 0;
                int* send_buffer = new int[send_buffer_size];
                int* receive_buffer = new int[receive_buffer_size];


                // Filling the buffer
                for (ModelPart::NodeIterator i_node = r_origin_nodes.begin(); i_node != r_origin_nodes.end(); ++i_node)
                    for (ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin(); i_dof != i_node->GetDofs().end(); i_dof++)
                        send_buffer[position++] = i_dof->EquationId();


                MPI_Status status;

                if (position > send_buffer_size)
                    std::cout << rank << " Error in estimating send buffer size...." << std::endl;


                int send_tag = i_color;
                int receive_tag = i_color;

                MPI_Sendrecv(send_buffer, send_buffer_size, MPI_INT, destination, send_tag, receive_buffer, receive_buffer_size, MPI_INT, destination, receive_tag,
                             MPI_COMM_WORLD, &status);

                // Updating nodes
                position = 0;
                for (ModelPart::NodeIterator i_node = r_dest_nodes.begin();
                     i_node != r_dest_nodes.end(); i_node++)
                    for (ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin(); i_dof != i_node->GetDofs().end(); i_dof++)
                    {
                        unsigned int NewId = static_cast<unsigned int>(receive_buffer[position++]);
                        if (NewId > i_dof->EquationId())
                            i_dof->SetEquationId(NewId);
                    }

                if (position > receive_buffer_size)
                    std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                delete [] send_buffer;
                delete [] receive_buffer;
            }

        // Once we are sure that the owner process of all nodes knows its equation id, replicate its value to all ghost images
        rComm.SynchronizeDofs();
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


    /*@} */

}; /* Class TrilinosBuilderAndSolverMLPeriodic */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

///@} //addtogroup

}  /* namespace Kratos.*/

#endif // KRATOS_TRILINOS_BUILDER_AND_SOLVER_ML_PERIODIC_H
